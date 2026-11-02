"""Small fixture: build outputs must not invalidate package-source identity."""
from pathlib import Path
from tempfile import TemporaryDirectory
from package_source_hash import source_fingerprint


with TemporaryDirectory() as directory:
    root = Path(directory)
    for name in ('src', 'R', 'inst/stan', 'data', 'man'):
        (root / name).mkdir(parents=True)
    inputs = ('DESCRIPTION', 'NAMESPACE', 'configure', 'src/Makevars',
              'src/model.h', 'src/model.cc', 'R/api.R', 'inst/stan/model.stan',
              'data/fixture.rda', 'man/api.Rd')
    for name in inputs:
        (root / name).write_text('input')
    before = source_fingerprint(root)
    for name in ('src/model.o', 'src/package.so', 'src/compiler.log',
                 'src/.cache/item', '00LOCK/package/item', '.git/index'):
        path = root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text('generated')
        assert source_fingerprint(root) == before, name
        path.write_text('changed generated')
        assert source_fingerprint(root) == before, name
    for name in inputs:
        path = root / name
        path.write_text('changed input')
        assert source_fingerprint(root) != before, name
        path.write_text('input')
        assert source_fingerprint(root) == before, name
    print('32 source-fingerprint assertions PASS')
