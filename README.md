## Hippocampal Development and Effects of Emotional Abuse

This repository contains the code used to collect data and generate analyses for the following paper:

*Choi, K., Min, J.J., Rah, Y.J., Kim, J.Y., ..., Lee, S.A. (in prep). 
Functional Development of Spatiotemporal Hippocampal Memory Binding in Adolescence and Its Alteration by Emotional Abuse.*

<br>

## Prerequisites

Before running the scripts, ensure you have the following software and toolboxes installed:

* **MATLAB R2024a** (or later)
* **Psychophysics Toolbox Version 3 (PTB-3)**: For task presentation
* **SPM12**: For fMRI preprocessing and GLM analysis
* **FreeSurfer**: For structural segmentation (Hippocampal subfields)
* **CONN Toolbox**: For functional connectivity analysis

<br>

## Scripts details

### 1. Episodic memory task
*The task was presented using MATLAB R2020b (MathWorks, USA) with Psychophysics Toolbox Version 3.*

- ``matlab/TASK_episodic_memory.m``

<br>

### 2. Behavioral analysis

- ``matlab/behavior_preprocessing.m``: preprocess raw behavioral data

- ``matlab/behavior_demographics.m``: demographics (age, sex, survey)

- ``matlab/behavior_memory.m``: episodic memory performance

<br>

### 3. fMRI analysis: GLM
*Structural and functional MRI data were preprocessed using SPM12 and FreeSurfer.*

- ``matlab/preprocessing.m``: fMRI preprocessing pipeline

- ``matlab/generate_regressors.m``: regressors for design matrix

- ``matlab/extract_ROIs.m``: ROI extraction from FreeSurfer output
<br>

*General Linear Model (GLM) analyses were done using SPM12.*

- ``matlab/GLM_first_level.m``: GLM individual level analysis

- ``matlab/GLM_second_level.m``: GLM group level analysis

- ``matlab/ROI_activation_analysis.m``: ROI-based activation analysis

<br>

### 4. fMRI analysis: connectivity
*Seed-based functional connectivity analyses were done using the CONN toolbox.*

- ``matlab/connectivity_script.m``: seed-based functional connectivity batch script for CONN toolbox

- ``matlab/connectivity_analysis.m``: statistical analysis of connectivity

<br>

## Contact
All scripts were written by **Kahyun Choi** (kahyu3@snu.ac.kr).

For any questions or issues, please contact via email.
