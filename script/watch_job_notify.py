#!/usr/bin/env python3
"""Send one aggregate email after every subjob in a SUMS run has terminated.

No analysis changes, resubmission, cancellation, or merge. Run on the submit host.
"""
import argparse
import collections
import datetime
import email.message
import email.policy
import email.utils
import fcntl
import getpass
import json
import os
from pathlib import Path
import re
import socket
import subprocess
import time
import xml.etree.ElementTree as ET


def now():
    return datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()


def save_json(path, data):
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w") as out:
        json.dump(data, out, ensure_ascii=False, indent=2)
        out.write("\n")
        out.flush()
        os.fsync(out.fileno())
    os.replace(str(tmp), str(path))


def parse_events(text):
    submitted, terminal = set(), {}
    # Ignore an incomplete last event while Condor is appending to the log.
    for block in re.split(r"(?m)^\.\.\.[ \t]*\r?$", text)[:-1]:
        match = re.match(r"\s*(\d{3}) \((\d+)\.(\d+)\.\d+\)", block)
        if not match:
            continue
        code, cluster, process = match.groups()
        key = "{}.{}".format(int(cluster), int(process))
        if code == "000":
            submitted.add(key)
            terminal.pop(key, None)
        elif code in ("001", "004", "012", "013"):
            terminal.pop(key, None)
        elif code == "005":
            result = re.search(r"Normal termination \(return value (\d+)\)", block)
            terminal[key] = "exit_" + result.group(1) if result else "signal_or_unknown"
        elif code == "009":
            terminal[key] = "removed"
    return submitted, terminal


def queue_states(jobid, submit_host):
    result = subprocess.run(
        ["condor_q", "-name", submit_host, "-constraint",
         'regexp("{}",Out)'.format(jobid), "-af", "ClusterId", "ProcId", "JobStatus"],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=45,
        check=True)
    states = {}
    for line in result.stdout.splitlines():
        fields = line.split()
        if len(fields) != 3 or not all(x.isdigit() for x in fields):
            raise RuntimeError("Unexpected condor_q output: " + line)
        states["{}.{}".format(int(fields[0]), int(fields[1]))] = int(fields[2])
    return states


def ready(expected, submitted, terminal, queue):
    return (len(submitted) == expected and set(terminal) == submitted
            and not queue)


def output_summary(xml, jobid, indices):
    output = xml.find("output").attrib
    if not output["toURL"].startswith("file:"):
        raise ValueError("Only local file: output is supported")
    root_dir = Path(output["toURL"][5:])
    template = Path(output["fromScratch"]).name
    log_template = xml.find("stdout").attrib["URL"]
    if not log_template.startswith("file:") or "$JOBID" not in template:
        raise ValueError("Unrecognized output/log template")
    missing, bad_header, incomplete = [], [], []
    roots = 0
    for i in indices:
        subjob = "{}_{}".format(jobid, i)
        path = root_dir / template.replace("$JOBID", subjob)
        if not path.is_file() or path.stat().st_size < 4:
            missing.append(i)
        else:
            with path.open("rb") as source:
                if source.read(4) == b"root":
                    roots += 1
                else:
                    bad_header.append(i)
        log_path = Path(log_template[5:].replace("$JOBID", subjob))
        log = log_path.read_text(errors="replace") if log_path.is_file() else ""
        matches = re.findall(r"requested=(\d+)\s+completed=(\d+)[^\n]*?status=(\d+)", log)
        if not matches or matches[-1][0] != matches[-1][1] or matches[-1][2] != "0":
            incomplete.append(i)
    return {"rootDirectory": str(root_dir), "rootsWithValidHeader": roots,
            "missingRoots": missing, "badRootHeaders": bad_header,
            "analysisCompletionUnconfirmed": incomplete}


def make_message(recipient, meta, status):
    expected = status["expectedSubjobs"]
    good = (status["outcomes"] == {"exit_0": expected}
            and status["output"]["rootsWithValidHeader"] == expected
            and not status["output"]["analysisCompletionUnconfirmed"])
    label = "完了" if good else "終了・要確認"
    message = email.message.EmailMessage(policy=email.policy.SMTP)
    message["From"] = "{}@{}".format(getpass.getuser(), socket.getfqdn())
    message["To"] = recipient
    message["Subject"] = "[STAR {}] {} ({} jobs)".format(label, meta["anaName"], expected)
    message["Date"] = email.utils.format_datetime(datetime.datetime.now().astimezone())
    message["Message-ID"] = "<star-{}-completion@{}>".format(meta["jobid"], socket.getfqdn())
    message.set_content(
        ("STAR解析の正常完了を確認しました。\n\n" if good else
         "バッチjobは全件終了しましたが、解析の正常完了は確認できません。結果を確認してください。\n\n") +
        "解析: {ana}\nSUMS ID: {jobid}\n確認時刻（日本時間）: {checked}\n"
        "終了job数: {n}\nCondor終了内訳: {outcomes}\n"
        "回収ROOT（ヘッダー確認済み）: {roots}/{n}\n"
        "出力先: {directory}\n\n"
        "未回収ROOTのsubjob番号: {missing}\n"
        "ROOTヘッダー異常: {bad}\n"
        "解析完了ログを確認できないsubjob番号: {incomplete}\n\n"
        "物理QA・ヒストグラム全体の破損検査・mergeは未実施です。\n"
        "これは指定した投入1回分の通知です。自動再投入は行いません。\n"
        "監視記録: {state}\n".format(
            ana=meta["anaName"], jobid=meta["jobid"], checked=status["checkedAt"],
            n=expected, outcomes=status["outcomes"],
            roots=status["output"]["rootsWithValidHeader"],
            directory=status["output"]["rootDirectory"],
            missing=status["output"]["missingRoots"], bad=status["output"]["badRootHeaders"],
            incomplete=status["output"]["analysisCompletionUnconfirmed"],
            state=status["stateDirectory"]))
    return message


def send_once(state_dir, message, status, recipient):
    attempt = state_dir / "delivery_attempt.json"
    if attempt.exists() or (state_dir / "sent.json").exists():
        raise SystemExit("Delivery already attempted; verify mail records before any manual retry")
    # Persist BEFORE invoking the MTA. A timeout may still mean mail was accepted.
    save_json(attempt, dict(email=recipient, attemptedAt=now(),
                            messageId=str(message["Message-ID"])))
    try:
        subprocess.run(["/usr/sbin/sendmail", "-i", "-t"],
                       input=message.as_bytes(), check=True, timeout=45)
        status.update(mailStatus="accepted_by_local_sendmail", email=recipient,
                      sentAt=now())
        save_json(state_dir / "sent.json", status)
    except Exception as error:
        raise SystemExit("Delivery outcome uncertain; automatic resend disabled: {}".format(error))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runmeta", type=Path, required=True)
    parser.add_argument("--event-log", type=Path, required=True)
    parser.add_argument("--email", required=True)
    parser.add_argument("--poll-sec", type=int, default=300)
    parser.add_argument("--settle-sec", type=int, default=300)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check", action="store_true", help="One read-only check; never send mail")
    mode.add_argument("--watch", action="store_true", help="Watch and send a single completion email")
    args = parser.parse_args()
    if not re.fullmatch(r"[A-Za-z0-9.!#$%&'*+/=?^_`{|}~-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}", args.email):
        parser.error("Invalid email address")
    if args.poll_sec < 30 or args.settle_sec < 0:
        parser.error("poll-sec must be >=30 and settle-sec >=0")
    meta = json.loads(args.runmeta.read_text())
    if not re.fullmatch(r"[0-9A-Fa-f]{32}", meta["jobid"]):
        parser.error("Invalid SUMS id")
    if socket.getfqdn() != meta["submit"]["host"]:
        parser.error("Run this monitor on the original submit host")
    root = Path(meta["projectRoot"])
    pattern = re.compile(re.escape(meta["jobPrefix"]) + r"_(\d+)\.list$")
    indices = sorted(int(pattern.fullmatch(p.name).group(1))
                     for p in (root / "job/run").iterdir() if pattern.fullmatch(p.name))
    if not indices or indices != list(range(len(indices))):
        parser.error("Missing/non-contiguous SUMS subjob lists")
    xml = ET.parse(meta["submit"]["joblistSnapshot"]["path"]).getroot()
    state_dir = root / "job/run/notify" / (meta["anaName"] + "_" + meta["jobid"])
    lock = None
    if args.watch:
        state_dir.mkdir(parents=True, exist_ok=True)
        lock = (state_dir / "monitor.lock").open("a")
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if (state_dir / "sent.json").exists():
            print("Already accepted by local mail service; no duplicate notification", flush=True)
            return
        if (state_dir / "delivery_attempt.json").exists():
            raise SystemExit("Previous delivery outcome uncertain; verify mail records before restarting")
        save_json(state_dir / "settings.json", dict(
            runmeta=str(args.runmeta.resolve()), eventLog=str(args.event_log.resolve()),
            email=args.email, pollSeconds=args.poll_sec, settleSeconds=args.settle_sec,
            expectedSubjobs=len(indices), pid=os.getpid(), host=socket.getfqdn(), startedAt=now()))
    settled_since = None
    while True:
        try:
            submitted, terminal = parse_events(args.event_log.read_text())
            queue = queue_states(meta["jobid"], meta["submit"]["host"])
            status = dict(checkedAt=now(), expectedSubjobs=len(indices),
                          submitted=len(submitted), terminal=len(terminal),
                          queued=len(queue), queueStates=dict(collections.Counter(queue.values())),
                          outcomes=dict(collections.Counter(terminal.values())),
                          stateDirectory=str(state_dir))
            complete = ready(len(indices), submitted, terminal, queue)
            if args.check:
                status["allTerminal"] = complete
                print(json.dumps(status, ensure_ascii=False, indent=2))
                return
            if complete:
                if settled_since is None:
                    settled_since = time.monotonic()
                if time.monotonic() - settled_since >= args.settle_sec:
                    status["output"] = output_summary(xml, meta["jobid"], indices)
                    message = make_message(args.email, meta, status)
                    (state_dir / "message.eml").write_bytes(message.as_bytes())
                    save_json(state_dir / "status.json", status)
                    send_once(state_dir, message, status, args.email)
                    print("{} Mail accepted by local sendmail for {}".format(now(), args.email), flush=True)
                    return
            else:
                settled_since = None
            save_json(state_dir / "status.json", status)
            print(json.dumps(status, ensure_ascii=False), flush=True)
        except (OSError, ValueError, subprocess.SubprocessError) as error:
            settled_since = None
            if args.check:
                raise
            print("{} Retry without sending: {}".format(now(), error), flush=True)
        time.sleep(args.poll_sec)


if __name__ == "__main__":
    main()
