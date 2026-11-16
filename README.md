# Hippocampal Memory Development and Emotional Abuse
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24272624.svg)](https://doi.org/10.5281/zenodo.24272624)
[![Code: MIT](https://img.shields.io/badge/Code-MIT-yellow.svg)](LICENSE)
[![Data: CC BY 4.0](https://img.shields.io/badge/Data-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

This repository contains the analysis code and data resources for the paper:

Choi, K., Min, J. J., Rah, Y. J., Kim, J. Y., Choi, S. Y., Hong, S. H., Jang, Y., Han, S., Lee, K. H., Kim, J.-W., Chung, G. H., Chung, D., Cha, J., Lee, S. A. (accepted). **Functional development of spatiotemporal hippocampal memory binding in adolescence and its disruption in emotional abuse**. *Nature Communications*.

## Requirements

- **MATLAB R2024a** or later, with Statistics and Machine Learning Toolbox and Image Processing Toolbox
- **Psychtoolbox-3**: task presentation
- **SPM12**: fMRI preprocessing and GLM analysis
- **FreeSurfer**: hippocampal segmentation
- **CONN**: functional connectivity analysis

## Code

### 1. Episodic memory task

The task was presented using MATLAB R2020b and Psychtoolbox-3.

- `matlab/TASK_episodic_memory.m`: task presentation

### 2. Behavioral analysis

- `matlab/behavior_preprocessing.m`: behavioral preprocessing
- `matlab/behavior_demographics.m`: demographics and survey measures
- `matlab/behavior_memory.m`: episodic memory performance and error patterns

### 3. fMRI analysis: GLM and ROIs

- `matlab/preprocessing.m`: fMRI preprocessing
- `matlab/generate_regressors.m`: trial- and run-level task regressors
- `matlab/extract_ROIs.m`: ROI extraction from FreeSurfer output
- `matlab/GLM_first_level.m`: first-level GLM analysis
- `matlab/GLM_second_level.m`: group-level GLM analysis
- `matlab/ROI_activation_analysis.m`: hippocampal activation, longitudinal peak position, retrieval similarity, and volume analyses

### 4. fMRI analysis: connectivity

- `matlab/connectivity_script.m`: CONN batch setup and first-level analysis
- `matlab/connectivity_analysis.m`: connectivity statistics and encoding-window sensitivity

## Data

Source data, trial-level behavioral accuracy, and group-level statistical maps are available in `data/`. See `data/README.md` for details.

## License

Code is released under the MIT License (`LICENSE`).
Data files are released under CC BY 4.0.

## Contact

Kahyun Choi (kahyu3@snu.ac.kr).
