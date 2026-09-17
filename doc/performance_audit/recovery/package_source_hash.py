"""Fingerprint package inputs without fingerprinting native build outputs."""
from pathlib import Path
import hashlib
import sys


def source_fingerprint(root):
    root = Path(root)
    files = []
    for name in ('DESCRIPTION', 'NAMESPACE', 'configure', 'configure.win',
                 '.Rbuildignore'):
        path = root / name
        if path.is_file():
            files.append(path)
    for name in ('R', 'src', 'inst', 'data', 'man'):
        directory = root / name
        if not directory.is_dir():
            continue
        for path in directory.rglob('*'):
            if not path.is_file():
                continue
            relative = path.relative_to(directory)
            if any(part.startswith('.') or part.startswith('00LOCK')
                   for part in relative.parts):
                continue
            if path.suffix.lower() in ('.o', '.so', '.dll', '.dylib', '.a',
                                       '.log', '.tmp', '.pyc'):
                continue
            files.append(path)
    fingerprint = hashlib.sha256()
    for path in sorted(files, key=lambda p: p.relative_to(root).as_posix()):
        content = hashlib.sha256(path.read_bytes()).hexdigest()
        fingerprint.update((path.relative_to(root).as_posix() + '\0' + content
                            + '\n').encode())
    return fingerprint.hexdigest()


if __name__ == '__main__':
    print(source_fingerprint(sys.argv[1]))
