import tarfile,gzip,io
from pathlib import Path
raw=Path('D:/cs2/data/raw/ds_geo_gse157007/GSE157007_RAW.tar')
with tarfile.open(raw) as t:
  for target in ['GSM4750298_F002_matrix.tsv.gz','GSM5684306_OH14_matrix.mtx.gz','GSM4750298_F002_features.tsv.gz']:
    f=t.extractfile(target)
    gz=gzip.GzipFile(fileobj=f)
    data=gz.read(500)
    print('\n--',target,'--')
    print(data.decode(errors='replace'))
