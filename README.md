# AliCPT beam systematics & deprojection — example notebooks

Example code for the dataset **“Simulated CMB dataset for beam systematics and
deprojection studies in AliCPT-1”** (Dou & Li). This repository
(<https://github.com/djzinchina/alicpt-deproj>) contains *code only*: it shows
how to read the beam-mismatch parameters and the simulated maps from the data
server, and how to reproduce the three Technical Validation figures of the data
paper.

> **No data are stored in this repository.** All FITS / TXT / PKL products live
> on the IHEP data server and are read from there; see [Data](#data).

<!-- TODO(owner): replace <DATA_DOI> once the dataset DOI is assigned. -->

## Contents

| Path | Description |
| --- | --- |
| `deproj_utils.py` | Small helpers: map plotting, IQU → T/E/B, E-mode recycling |
| `examples/plot_beam_params.ipynb` | Figure 2 — recovered beam-mismatch parameters |
| `examples/plot_maps_and_spectra.ipynb` | Figures 3 and 4 — B-mode maps and TT/EE/TE/BB power spectra |
| `examples/generate_noise_maps.ipynb` | Gaussian white-noise realizations from the noise r.m.s. maps |
| `scripts/check_data.py` | Check that `DEPROJ_DATA_ROOT` holds the expected products |
| `environment.yml` | Conda environment used by the examples |
| `pyproject.toml` | Optional `pip install -e .` metadata for `deproj_utils` |
| `CITATION.cff` | Citation metadata |

## Data

The dataset is hosted on the AliCPT shared filesystem at IHEP; see the data paper for access details. The examples read the
data through the `DEPROJ_DATA_ROOT` environment variable and never write back:

```bash
export DEPROJ_DATA_ROOT=/sharefs/alicpt/users/doujiazheng/data/deproj
```

The examples expect this layout (only the parts they use are listed):

```
$DEPROJ_DATA_ROOT/
├── mask/
│   └── AliCPT_20uKcut150_C_1024.fits      analysis mask, Nside = 1024
├── maps/
│   ├── LensWithBeamSys/<sim>/             <sim> = 000 … 299
│   │   ├── LensInfo/
│   │   │   ├── CMB_LEN.fits               lensed CMB, IQU, beam-free
│   │   │   └── CMB_UNLEN.fits             unlensed CMB, IQU, beam-free
│   │   ├── Ali_95.fits                    AliCPT 95 GHz, with beam systematics
│   │   ├── Ali_150.fits                   AliCPT 150 GHz, with beam systematics
│   │   ├── DEPROJ_Ali_95.fits             AliCPT 95 GHz, deprojected
│   │   ├── DEPROJ_Ali_150.fits            AliCPT 150 GHz, deprojected
│   │   └── HFI_*.fits, WMAP_K.fits        ancillary channels
│   └── Noise/
│       └── SIGMA_<channel>.fits           IQU noise r.m.s. maps
├── FG/
│   └── group1_map_<freq>GHz.fits          Galactic foregrounds, IQU
└── params/
    ├── deproj_beam_params.pkl             fitted parameters, dict[95, 150]
    ├── beamsys/
    │   └── <mode>_<freq>G.txt             injected parameters (r, theta, value)
    └── fp/
        ├── r.npy, theta.npy, chi.npy      focal plane, per detector pair
        └── fp_r.npy, fp_theta.npy, fp_chi.npy   per detector
```

All maps are HEALPix `Nside = 1024`, RING ordered, in equatorial coordinates
and in `uK_CMB`. The `Ali_*` and `DEPROJ_Ali_*` maps are the 95 and 150 GHz
AliCPT bands; `LensWithBeamSys/<sim>/LensInfo/CMB_LEN.fits` is the beam-free
lensed CMB used to build the “no beam systematics” comparison map. See the
data paper for the exact content (CMB, foregrounds and noise) of each product.

To verify that all the files used by the examples are present, run:

```bash
python scripts/check_data.py
```

## Installation

The examples need `numpy`, `healpy`, `pymaster` and `matplotlib`:

```bash
conda env create -f environment.yml
conda activate alicpt-deproj
```

A pip-only alternative is

```bash
pip install numpy healpy pymaster matplotlib jupyterlab
```

The notebooks import `deproj_utils` from the repository root; either run them
in place (they add `..` to `sys.path`) or install the module:

```bash
pip install -e .
```

## Usage

```bash
export DEPROJ_DATA_ROOT=/sharefs/alicpt/users/doujiazheng/data/deproj
jupyter lab examples/
```

Each notebook writes its figures to `examples/output/` and keeps the numerical
results (spectra, parameters) in memory so that they can be inspected
interactively. No data or spectra files are inserted into the repository.

| Notebook | Reproduces | Output |
| --- | --- | --- |
| `plot_beam_params.ipynb` | Figure 2, plus input-vs-recovered for all detector pairs | `beam_params_150GHz.png`, `beam_params_scatter_150GHz.png` |
| `plot_maps_and_spectra.ipynb` | Figure 3 and Figure 4 | `Bmaps_150GHz.png`, `Dl_150GHz.png` |
| `generate_noise_maps.ipynb` | noise model of the dataset | `analysis_mask.png`, `noise_ALI_150.png` |

The band, detector pair and realization are set at the top of each notebook
(`FREQ`, `IDET`, `SIM`); both 95 and 150 GHz are supported.

## Implementation notes

* `deproj_utils.py` contains trimmed versions of three routines that were
  originally part of [cmbframe](https://github.com/1cosmologist/CMBframe):
  `plot_maps`, `iqu2teb` and `clean_bmap` (the `Emode_recycler` class). They
  were rewritten so that the examples do not depend on the full framework, and
  the hard-coded HEALPix search path of the original code was removed. No
  cmbframe source file is redistributed here.
* The spherical-harmonic transforms use `use_weights=False`, so no external
  HEALPix weight file is required.
* The power spectra are computed on the fly with NaMaster (apodized mask,
  decoupled bandpowers, beam deconvolved) and binned with a constant bin size
  of 30 up to `lmax = 700`. No spectrum files are part of the dataset.
* `deproj_beam_params.pkl` has shape `(106, 6, 1728)` per band: 106 time
  trunks, the six differential modes `(dg, dx, dy, dsigma, dp, dc)` and 1728
  detector pairs. `dx`/`dy` are in radians, `dsigma` in degrees.
* The `dp` parameter is biased by about `7 × 10^-3` because of the intrinsic
  CMB T–E correlation; this is expected and is corrected in the analysis
  (see Figure 2 of the data paper).

## Citation

If you use the dataset or this code, please cite the data paper

```bibtex
@article{Dou2026alicptdeproj,
  title  = {Simulated CMB dataset for beam systematics and deprojection studies in AliCPT-1},
  author = {Dou, Jiazheng and Li, Hong},
  year   = {2026},
  doi    = {<DATA_DOI>}
}
```

## License

GNU General Public License v3.0 or later — see [LICENSE](LICENSE).

## Contact

Jiazheng Dou — doujzh@ihep.ac.cn
