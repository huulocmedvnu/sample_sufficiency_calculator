import os
import sys

_HERE = os.path.dirname(__file__)
# make `calculator` (src/) and `verify_theory` (tests/) importable in tests
for _p in ("src", "tests"):
    _path = os.path.join(_HERE, _p)
    if _path not in sys.path:
        sys.path.insert(0, _path)
