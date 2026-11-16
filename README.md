# Hippocampal Memory Development and Emotional Abuse

Code and data accompanying:

Choi, K., Min, J. J., Rah, Y. J., Kim, J. Y., Choi, S. Y., Hong, S. H., Jang, Y., Han, S., Lee, K. H., Kim, J.-W., Chung, G. H., Chung, D., Cha, J., Lee, S. A. (accepted). Functional development of spatiotemporal hippocampal memory binding in adolescence and its disruption in emotional abuse. *Nature Communications*.

## Requirements

- **MATLAB R2024a** or later, with Statistics and Machine Learning Toolbox and Image Processing Toolbox
- **Psychtoolbox-3**: task presentation
- **SPM12**: fMRI preprocessing and GLM analysis
- **FreeSurfer**: hippocampal segmentation
- **CONN**: functional connectivity analysis

## Code

### 1. Episodic memory task

The task was presented using MATLAB R2020b and Psychtoolbox-3.

- [TASK_episodic_memory.m](matlab/TASK_episodic_memory.m): task presentation

### 2. Behavioral analysis

- [behavior_preprocessing.m](matlab/behavior_preprocessing.m): behavioral preprocessing
- [behavior_demographics.m](matlab/behavior_demographics.m): demographics and survey measures
- [behavior_memory.m](matlab/behavior_memory.m): episodic memory performance and error patterns

### 3. fMRI analysis: GLM and ROIs

- [preprocessing.m](matlab/preprocessing.m): fMRI preprocessing
- [generate_regressors.m](matlab/generate_regressors.m): trial- and run-level task regressors
- [extract_ROIs.m](matlab/extract_ROIs.m): ROI extraction from FreeSurfer output
- [GLM_first_level.m](matlab/GLM_first_level.m): first-level GLM analysis
- [GLM_second_level.m](matlab/GLM_second_level.m): group-level GLM analysis
- [ROI_activation_analysis.m](matlab/ROI_activation_analysis.m): hippocampal activation, longitudinal peak position, retrieval similarity, and volume analyses

### 4. fMRI analysis: connectivity

- [connectivity_script.m](matlab/connectivity_script.m): CONN batch setup and first-level analysis
- [connectivity_analysis.m](matlab/connectivity_analysis.m): connectivity statistics and encoding-window sensitivity

## Data

Source data, trial-level behavioral accuracy, and group-level statistical maps are available in [data/](data/). See the [data description](data/README.md).

## Contact

Kahyun Choi (kahyu3@snu.ac.kr).
