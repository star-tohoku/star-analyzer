#!/bin/bash
# Retry one local reference without rerunning successful species outputs.
set -euo pipefail
[[ $# == 2 ]] || { echo 'Usage: bash tests/run_femto_lambda_reference.sh 13p5|3p85 legacy|standalone' >&2; exit 2; }
ref_energy="$1"
ref_kind="$2"
case "$ref_energy" in 13p5|3p85) ;; *) exit 2;; esac
case "$ref_kind" in legacy|standalone) ;; *) exit 2;; esac
ref_project=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)
cd "$ref_project"
ref_base=rootfile/femto_lambda_kf_validation_20260930
ref_output="$ref_base/$ref_energy/$ref_kind.root"
ref_log="$ref_base/provenance/${ref_energy}_${ref_kind}_retry1.log"
ref_status="$ref_base/provenance/${ref_energy}_${ref_kind}_retry1.exitcode"
[[ ! -e "$ref_output" && ! -e "$ref_log" && ! -e "$ref_status" ]] || { echo 'Existing retry output/log; refusing overwrite' >&2; exit 2; }
trap 'ref_exit=$?; printf "%s\n" "$ref_exit" > "$ref_status"' EXIT
ref_input="config/picoDstList/auau${ref_energy}_femtoLambda_validation_20260930.list"
if [[ "$ref_kind" == legacy ]]; then
  ref_config="config/mainconf/main_auau${ref_energy}_anaLambdaNuclearId.yaml"
  ref_macro=run_anaLambdaNuclearId
else
  ref_config="config/mainconf/main_auau${ref_energy}_KFParticle_femto_control_20260930.yaml"
  ref_macro=run_anaLambda_KFParticle
fi
bash tests/run_femto_lambda_root5.sh "$ref_config" \
  "analysis/${ref_macro}.C(\"$ref_input\",\"$ref_output\",\"reference_${ref_energy}_${ref_kind}\",10000,\"$ref_config\")" \
  > "$ref_log" 2>&1
if [[ "$ref_kind" == legacy ]]; then
  rg -q 'StLambdaMaker::Finish\(\) processed 10000 events' "$ref_log"
fi
[[ -s "$ref_output" ]] || { echo 'Reference output is missing' >&2; exit 1; }
