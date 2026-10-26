#!/usr/bin/env python3
"""Exercise the generated worker input contract with csh; no farm or analysis run."""
from pathlib import Path
import subprocess,tempfile,unittest,xml.etree.ElementTree as ET
ROOT=Path(__file__).resolve().parents[1]
class BatchInput(unittest.TestCase):
 def invoke(self,contents):
  with tempfile.TemporaryDirectory() as temporary:
   d=Path(temporary);p=d/'one.picoDst.root';p.write_text('readable fixture');listing=d/'input.list';listing.write_text(contents.replace('INPUT',str(p)));cmd=ET.parse(ROOT/'job/joblist/job_template_from_conf.xml').getroot().findtext('command');start=cmd.index('    if ( ! $?SCRATCH )');end=cmd.index('    set runtime_dir=');block=cmd[start:end].replace('__FILES_PER_JOB__','3');script=d/'check.csh';script.write_text('setenv SCRATCH "'+str(d)+'"\nsetenv FILELIST "'+str(listing)+'"\n'+block)
   r=subprocess.run(['csh','-f',str(script)],stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True);normalized=d/'production_input.list';return r.returncode,normalized.read_text() if normalized.exists() else '',str(p)
 def test_path_only(self):
  rc,out,path=self.invoke('INPUT\n');self.assertEqual(rc,0);self.assertEqual(out,path+'\n')
 def test_sums_path_count(self):
  rc,out,path=self.invoke('INPUT 4545\n');self.assertEqual(rc,0);self.assertEqual(out,path+'\n')
 def test_mixed_lines(self):
  rc,out,path=self.invoke('# comment\nINPUT 4545\n\nINPUT\n');self.assertEqual(rc,0);self.assertEqual(out,(path+'\n')*2)
 def test_invalid_count_fails(self):self.assertEqual(self.invoke('INPUT bad\n')[0],83)
 def test_extra_column_fails(self):self.assertEqual(self.invoke('INPUT 1 extra\n')[0],83)
 def test_unstaged_fails(self):self.assertEqual(self.invoke('/missing/pico.root 123\n')[0],84)
 def test_empty_fails(self):self.assertEqual(self.invoke('')[0],83)
if __name__=='__main__':unittest.main()
