#!/usr/bin/env python
"""Check that DEPROJ_DATA_ROOT holds the data products used by the examples.

Usage:  python scripts/check_data.py
"""
import os
import pathlib
import sys

DATA = pathlib.Path(os.environ.get(
    'DEPROJ_DATA_ROOT', '/sharefs/alicpt/users/doujiazheng/data/deproj'))

REQUIRED = [
    'mask/AliCPT_20uKcut150_C_1024.fits',
    'params/deproj_beam_params.pkl',
    'params/beamsys/dg_150G.txt',
    'params/beamsys/de_conti_150G.txt',
    'params/fp/theta.npy',
    'params/fp/chi.npy',
    'maps/LensWithBeamSys/000/LensInfo/CMB_LEN.fits',
    'maps/LensWithBeamSys/000/Ali_150.fits',
    'maps/LensWithBeamSys/000/DEPROJ_Ali_150.fits',
    'maps/Noise/SIGMA_ALI_150.fits',
    'FG/group1_map_150GHz.fits',
]

missing = []
for rel in REQUIRED:
    ok = (DATA / rel).exists()
    missing += [] if ok else [rel]
    print('%s  %s' % ('  OK  ' if ok else 'MISSING', rel))

print()
if missing:
    print('Data root %s is incomplete: %d file(s) missing.' % (DATA, len(missing)))
    sys.exit(1)
print('Data root %s looks complete.' % DATA)
