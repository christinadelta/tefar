# TEFAR — TMS-EEG / EEG FieldTrip Artefact Removal

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21844563.svg)](https://doi.org/10.5281/zenodo.21844563)

TEFAR is a lightweight, configurable MATLAB framework for semi-automatic, ICA-based artefact-component classification in EEG and TMS–EEG, built on [FieldTrip](https://www.fieldtriptoolbox.org/). A single scoring engine (`tefar_core`) is exposed through two profile wrappers, one for ordinary EEG, one for TMS–EEG, and detects line noise, blinks, lateral eye movements, cranial muscle, cardiac activity, and the TMS-specific decay and recharge transients from established independent-component signatures.

TEFAR is **semi-automatic by design**: it suggests components for removal and reports why, but the final decision rests with the researcher. Visual inspection of component topographies, time courses, and spectra remains an essential step of the workflow.

**Key properties**

- One engine, two profiles: `TEFAR_eeg` for resting/task EEG, `TEFAR_tms` for TMS–EEG.
- FieldTrip is the **only** dependency. Kurtosis, robust statistics, autocorrelation, and the periodogram are implemented internally. No MATLAB toolboxes required.
- Robust (median/MAD) thresholds that remain stable even when a single artefact dominates the decomposition.
- Fully configurable through a single `cfg` structure: every band, threshold, window, and channel set has a documented, overridable default.
- Trial-structure agnostic: temporal detectors locate their windows from each component's own time axis.
- Ground-truth validated on simulated data, with a self-contained, reproducible benchmark that runs in base MATLAB/Octave.

## Requirements

- MATLAB R2023b or later (tested on R2023b) 
- [FieldTrip](https://www.fieldtriptoolbox.org/download/) (tested with fieldtrip-20250414) on the MATLAB path
- No additional MATLAB toolboxes
- The simulation and logic-verification scripts (`simulate_tefar_data`, `verify_tefar_logic`) run in base MATLAB or GNU Octave, without FieldTrip

## Installation

```matlab
% 1. Get the code
%    git clone https://github.com/christinadelta/tefar.git
% 2. Add TEFAR and FieldTrip to your path
addpath('/path/to/tefar');
addpath('/path/to/fieldtrip');
ft_defaults
```

## Quick start

```matlab
% Ordinary EEG (FieldTrip data structure from ft_preprocessing):
[comp, artifacts] = TEFAR_eeg(data);

artifacts.reject    % components suggested for removal
artifacts.score     % number of detectors flagging each component
artifacts.metrics   % all per-component metric vectors

% Inspect before rejecting (always!), then:
cfg           = [];
cfg.component = artifacts.reject;
data_clean    = ft_rejectcomponent(cfg, comp);
```

For a complete worked example on simulated data (including topoplots of the flagged components) see [`demo_tefar_eeg.m`](demo_tefar_eeg.m) and [`demo_tefar_tms.m`](demo_tefar_tms.m).

## Usage

```matlab
% TMS-EEG (per-detector component lists as separate outputs):
[comp, line_c, musc_c, decay_c, addmusc_c, rech_c, blink_c] = ...
    TEFAR_tms(data_filtered, trl);

% ...or return the full result (scores, metrics, populated cfg):
[comp, ~,~,~,~,~,~, artifacts] = TEFAR_tms(data_filtered, trl);

% Override any default (e.g. US mains, stricter line rule, score-based rejection):
cfg = struct('line_freq',60, 'reject_rule','score', 'reject_threshold',2);
cfg.line.z = 6;
[comp, artifacts] = TEFAR_eeg(data, cfg);
```

`reject_rule` is `'union'` by default (a component flagged by any detector is suggested for removal). Set `'score'` to require agreement across at least `reject_threshold` detectors. The full configuration reference — every detector's bands, thresholds, and windows — is documented in the header of [`tefar_core.m`](tefar_core.m).

## Detectors

Each detector targets an established independent-component artefact signature:

| Detector | Feature | Rule | Profile |
|---|---|---|---|
| `line` | mains fundamental / sideband prominence | robust-z > 5 | both |
| `muscle` | high-frequency power fraction (mains excluded) | robust-z > 5 and fraction > 0.5 | both |
| `muscle_topo` | spatial focality (peak / Σ\|topo\|) | robust-z > 5 | both |
| `blink` | temporal kurtosis + frontal topography | kurt > 4 and frontal ratio > 0.5 | both |
| `eyemove` | fronto-lateral anti-symmetry | asym > 0.15 and frontal ratio > 0.5 | EEG |
| `cardiac` | autocorrelation periodicity + spikiness | robust-z > 5, autocorr > 0.15, kurt > 3.5 | both |
| `decay` | early post-pulse RMS vs baseline SD | RMS > 2 × baseline SD | TMS |
| `recharge` | later-window RMS vs baseline SD | RMS > 4 × baseline SD | TMS |

Two design choices matter in practice:

- **Robust z-scores (median/MAD)** replace the conventional `mean + k·SD` rule. When a single component carries almost all of an artefact (typical for line noise or a large muscle component), the outlier inflates its own mean and SD, and the maximum achievable classic z-score is only `(N−1)/√N`. Median/MAD thresholds are immune to this.
- **Temporal detectors are trial-structure agnostic**: the decay/recharge windows are located from `comp.time` per trial, so no epoching convention is assumed.

## Validation

TEFAR is validated on ground-truth simulations generated under the ICA model `x = A·s`, in which every latent source carries a known label, so every flag can be scored against the truth.

| Profile | Sensitivity | Specificity | Accuracy |
|---|---|---|---|
| TMS (30 random seeds) | 0.956 | 1.000 | 0.985 |
| EEG (30 random seeds) | 0.973 | 1.000 | 0.992 |

Near-perfect **specificity** is the property that matters most: brain components are essentially never flagged for removal. The occasional missed artefact is precisely what the semi-automatic visual-inspection step exists to catch.

Reproduce the validation with:

```matlab
verify_tefar_logic('tms');  verify_tefar_logic('eeg');  % toolbox-free detector check
validate_tefar('tms','components');   % detectors on known sources (needs FieldTrip)
validate_tefar('tms','endtoend');     % full pipeline: mix -> FastICA -> classify
make_tefar_figures                    % regenerates the figures in the paper
```

## Repository contents

| File | Purpose |
|---|---|
| `tefar_core.m` | The engine: runs ICA (FastICA by default) and scores every component against the configured detectors. |
| `TEFAR_eeg.m` | EEG profile wrapper (resting/task data, no TMS). |
| `TEFAR_tms.m` | TMS–EEG profile wrapper. |
| `demo_tefar_eeg.m` | Worked example: simulate EEG, classify, plot flagged components. |
| `demo_tefar_tms.m` | Worked example: simulate TMS–EEG, classify, plot flagged components. |
| `simulate_tefar_data.m` | Ground-truth simulation (`x = A·s`) with labelled sources; base MATLAB/Octave only. |
| `validate_tefar.m` | Scores TEFAR against the simulation via the real `tefar_core` (needs FieldTrip). |
| `verify_tefar_logic.m` | Independent, FieldTrip-free reimplementation of the detector maths for cross-checking. |
| `make_tefar_figures.m` | Regenerates the validation figures reported in the paper. |
| `CITATION.cff` | Citation metadata (used by GitHub's "Cite this repository" and Zenodo). |

## Citation

If you use TEFAR in your work, please cite the archived software version (v1.0.0 DOI: 10.5281/zenodo.22706419) or use the Cite this repository button above, which draws on CITATION.cff. The DOI 10.5281/zenodo.21844563 represents all versions and always resolves to the latest release.

A methods paper describing TEFAR is in preparation; this section will be
updated when the preprint is available.

## License

Released under the [MIT License](LICENSE).

## Contact & contributions

Bug reports and feature requests are welcome via [GitHub issues](https://github.com/christinadelta/tefar/issues). TEFAR is under active development — feedback from real-world EEG and TMS–EEG pipelines is especially appreciated.
