import tarfile, pandas as pd
from pathlib import Path
raw=Path('D:/cs2/data/raw/ds_geo_gse157007/GSE157007_RAW.tar')
with tarfile.open(raw) as t:
    names=t.getnames()
print('\n'.join(names[:120]))
print('n',len(names))
lab=Path('D:/cs1/runs/cs1-paper-factory/run_20260529_095727/agents/02-knowledge/workspace/gse157007_confirmed_group_assignments.tsv')
print('labels exists', lab.exists())
print(pd.read_csv(lab, sep='\t').to_string())
