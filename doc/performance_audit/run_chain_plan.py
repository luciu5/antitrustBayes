import csv,subprocess,os
from pathlib import Path
root=Path('doc/statistical_toolkit_runs/20260916/performance-audit')
env=os.environ.copy();env['R_LIBS']='/tmp/toolkit-publish-Rlib:/tmp/cga-Rlib'
for row in csv.DictReader((root/'chain_plan.csv').open()):
 out=root/row['case_id']
 if Path(str(out)+'.csv').exists():
  print('Preserved completed',row['case_id'],flush=True);continue
 print('Running',row['case_id'],flush=True)
 with Path(str(out)+'.log').open('w') as log:
  p=subprocess.run(['python3','doc/performance_audit/process_memory.py',str(out)+'.memory.csv','Rscript','doc/performance_audit/chain_benchmark.R',row['fixture'],row['cores'],row['seed'],str(out)],env=env,stdout=log,stderr=subprocess.STDOUT)
 print('Finished',row['case_id'],'exit',p.returncode,flush=True)
 if p.returncode:
  print('Failure preserved; continuing frozen plan without tuning.',flush=True)
