# OOT analysis (JULIA / Anton Paar DSC)

Oxidation onset temperature of polyolefins by dynamic DSC. The onset is found by
the extrapolated-tangent construction on the temperature axis (ASTM E2009,
dynamic variant of ISO 11357-6).

## Usage

Put all `.m` files of this folder on the MATLAB path (or `cd` into it) and run

```matlab
OOT_analysis
```

1. A file dialog opens. Select one or more exported `*.csv` files (Ctrl/Shift-click
   for several). The dialog remembers the last folder.
2. A run table opens. The **sample label** is read from the file name: the text between
   the dd/mm/yyyy date and `_OOT`.

   `20261024_0451_Measurement_JULIA_DSC_23102026_HDPE_M1053_virgin_OOT.csv` → `HDPE_M1053_virgin`

   Runs with the same label are treated as replicates (R1, R2, … in acquisition order).
   You can edit the label, replicate ID and polymer, untick runs, add more files and
   choose the output folder.
3. Press **Run analysis**.

Non-interactive (scripts / batch):

```matlab
d = dir('data/*_OOT*.csv');
R = OOT_analysis('Files', fullfile({d.folder}, {d.name}), 'Interactive', false);
```

## Output (default `<data folder>/OOT_results`)

| File | Content |
|---|---|
| `OOT_summary.csv` | one row per run: OOT, departure, inflection, peak, dHox, melting check, … |
| `OOT_by_sample.csv` | replicate mean and half-range per sample |
| `OOT_QC.txt` | QC notes (mass spread, open integrals, implausible crystallinity, …) |
| `fig/Fig1_OOT_ramps` | oxidising ramps, all runs overlaid |
| `fig/Fig2_onset_region` | onset region, ordinate expanded |
| `fig/Fig3_tangent_construction` | tangent construction, one panel per run |
| `fig/OOT_<label>_<rep>` | the same construction as a stand-alone figure per run, with the file name as footer |
| `fig/Fig4_OOT_by_sample` | OOT per sample, replicates and mean (only with 2 or more runs) |

Figures are saved as PNG (300 dpi) and vector PDF.

## Files

| File | Purpose |
|---|---|
| `OOT_analysis.m` | main entry point: UI, analysis loop, QC, figures |
| `readJuliaDSC.m` | fast reader for the JULIA export (header + one `textscan`) |
| `parseRunName.m` | sample label / acquisition time / polymer from the file name |
| `computeOOT.m` | baseline, tangent and onset on a uniform 0.1 K grid |
| `meltingCheck.m` | melting enthalpy and crystallinity of the inert pre-heat, used as a check on the entered mass |

Settings (palette, reference enthalpies per polymer, detection thresholds, QC
tolerances) are in `defaultConfig()` at the end of `OOT_analysis.m`.

Base MATLAB only, R2019b or later (`tiledlayout`). `exportgraphics` is used when
available (R2020a+), otherwise `print`.
