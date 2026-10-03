# -*- coding: utf-8 -*-
# Defines data paths from the same environment variables as _config.R (README section 4).
# Imported by the Python scripts.

import os
ROOT = os.environ.get("RICE_ROOT", os.getcwd())
DATA_ROOT = os.environ.get("RICE_DATA", ROOT)
GRID_ROOT = os.environ.get("RICE_GRID_ROOT")
os.makedirs(ROOT, exist_ok=True)
if not os.path.isdir(DATA_ROOT):
    raise SystemExit("RICE_DATA does not exist: %s" % DATA_ROOT)
CODE = os.path.dirname(os.path.abspath(__file__))
while not os.path.exists(os.path.join(CODE, "_config.py")) and os.path.dirname(CODE) != CODE:
    CODE = os.path.dirname(CODE)
_meta_default = os.path.join(DATA_ROOT, "00_META")
META_DIR = os.environ.get("RICE_META", _meta_default if os.path.isdir(_meta_default) else os.path.join(CODE, "meta"))

