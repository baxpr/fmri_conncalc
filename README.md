# fmri_conncalc
Preprocessing and functional connectivity computation for fMRI

## Quickstart

Here is an example for the "jsins" version of the processor, as described in
`conncalc_jsins_v1.0.0.yaml`.

```
singularity
  run
  --bind <INDIR>:/INPUTS
  --bind <OUTDIR>:/OUTPUTS
  baxpr-fmri_conncalc-master-v1.0.0.simg
  magick_path /usr/bin
  param_file params_JSins.csv
  wroi_file rois_JSins.nii.gz
  roi_file ''
  roiinfo_file rois_JSins.csv
  coregmat_file /INPUTS/coreg_mat.txt \
  deffwd_file /INPUTS/y_deffwd.nii.gz \
  ct1_file /INPUTS/ct1.nii.gz \
  wgm_file /INPUTS/wgm.nii.gz \
  wcseg_file /INPUTS/wcseg.nii.gz \
  func_file /INPUTS/fmri.nii.gz \
  project PROJECT_LABEL \
  subject SUBJECT_LABEL \
  session SESSION_LABEL \
  scan SCAN_LABEL \
  out_dir /OUTPUTS
```

The outputs are:

```
fmri_conncalc.pdf    Report
params.csv           Parameters used in the analysis
FD.txt               Framewise displacement time series
DVARS.txt            Framewise variance time series
badvols.txt          Scrubbed volumes indicator time series
rp_adfunc.txt        Realignment (motion) values
wmeanadfunc.nii.gz   Mean functional image in standard space
wadfunc.nii.gz       Slice time corrected and realigned functional images in standard space
rroi_labels.nii.gz   Region of interest label image
roi_snr.nii.gz       ROI SNR image
roi_info.csv         ROI info
roi_labels.csv       ROI names (if available)

Series of results repeated for each of the four processing streams
(keep or remove mean gray matter; scrub or no scrub):

  confounds_removegm_noscrub.txt               Confound (filter) matrix
  connectivity_matrix_R_removegm_noscrub.csv   Connectivity matrix
  filtered_removegm_noscrub.nii.gz             Filtered functional images
  roi_timeseries_removegm_noscrub.csv          Filtered ROI time series
  stats_removegm_noscrub.txt                   Various statistics
  Zmap_removegm_noscrub.nii.gz                 Unsmoothed ROI connectivity maps
  sZmap_removegm_noscrub.nii.gz                Smoothed ROI connectivity maps
```


## Dependencies

The built singularity container `baxpr-fmri_conncalc-master-v1.0.0.simg` (URL is shub://baxpr/fmri_conncalc:v1.0.0) is stand-alone with no external dependencies. The compiled matlab `bin/run_fmri_conncalc.sh` requires only the appropriate MATLAB Runtime to execute. To build these there are two stages:

1. Compile the MATLAB code into a stand-alone executable, using `compile_matlab.sh`. This requires a full MATLAB installation (R2017a, v92) and SPM12 (https://www.fil.ion.ucl.ac.uk/spm/).

2. Build the singularity container. In addition to a few specific OS packages, this requires the MATLAB Compiled Runtime. All are specified to be downloaded during the build in the singularity recipe `Singularity.v1.0.0`. The container help text gives build instructions. Alternatively the built container can be obtained from singularity-hub:
   `singularity pull shub://baxpr/fmri_conncalc:v1.0.0`


## Peculiarities of specific pipelines

Some critical analysis parameters are specified in the `param_file`, e.g. `params_JSins.csv`. This is a reference to a file that's in the built container, but these can also be viewed in the code repository e.g. `src/params/params_JSins.csv`. The parameters get as detailed as the repetition time of the fMRI scans. If the needed parameter file is not in the container already:
- Add the new parameter file ``src/params'
- Update the matlab compilation code to include it
- Recompile the matlab
- Commit to github
- Rebuild the container (increment the patch version e.g. 1.0.0 to 1.0.1)
- Create an updated YAML file appropriate for the parameter set

### jsins version

`conncalc_jsins_v1.0.0.yaml`

Standard space regions of interest are used, `params/JS_insula/rois_JSins.nii.gz`, identical for every subject.

Connectivity matrix is computed (Pearson bivariate correlation R). A connectivity map is computed for each ROI (Fisher Z transform applied to Pearson bivariate correlation). Spatial smoothing is applied to the connectivity maps only.

Parameter settings in `params_JSins.csv`:

- FMRI repetition time (TR) is assumed to be 2.000 sec
- Use all fMRI volumes (none dropped)
- No slice timing correction
- 6mm FWHM Gaussian spatial smoothing applied to connectivity maps
- Filter settings (confound regressor matrix):
  * 0.01 Hz - 0.10 Hz bandpass filter (Fourier basis)
  * 6 motion parameters (translation and rotation)
  * 6 first differences of motion parameters
  * First 6 principal components of voxel time series from the eroded white matter/CSF compartment
- For scrubbed results, volumes before and after an excursion of FD > 0.5 are removed. DVARS is not used for scrubbing.
- Connectivity maps are saved for each ROI.


### szhab version

No YAML available yet.

Subject-specific regions of interest are used, as described in the native space ROI image supplied as input. This image must be in the same space as the subject's native space structural.

Connectivity matrix is computed (Pearson bivariate correlation R of filtered time series). Spatial smoothing is not used.

Parameter settings in `params_SZhab.csv`:

- FMRI repetition time (TR) is assumed to be 2.000 sec
- 5 initial volumes are dropped, and the following 60 volumes are used for the analysis
- No slice timing correction
- Filter settings (confound regressor matrix):
  * 0.01 Hz - 0.15 Hz bandpass filter (Fourier basis)
  * 6 motion parameters (translation and rotation)
  * First 3 principal components of voxel time series from the eroded white matter/CSF compartment
- For scrubbed results, volumes before and after an excursion of FD > 0.5 are removed. DVARS is not used for scrubbing.


## General pipeline

Other than the above, processing proceeds as follows.

1. Drop functional volumes as specified.

2. Perform slice timing correction as specified. (SPM12 slice timing correction)

3. Perform motion realignment: two-stage alignment to mean image. (SPM12 realignment)

4. Coregister the mean functional image to the T1 weighted structural using a rigid body transform. The structural is first skull-stripped by zeroing all voxels that were not labeled by the multiatlas segmentation. The transformation is then applied to all functional volumes. (SPM12 coregistration)

5. Quality parameters are computed: framewise displacement FD and framewise signal variance DVARS. Volumes exceeding scrubbing criteria are marked ("badvols").

5. The functional and structural images are warped to standard space using the supplied nonlinear transform (forward deformation image). (SPM12 deformation tools)

6. The supplied standard space ROI image file is resampled to match the standard space fMRI geometry. (SPM12 reslice)

8. Connectivity computation. All filtering is done in a single step: a design matrix of confounds is created (see lists above), it is fit to each voxel time series, and the residuals are extracted. Then bivariate Pearson correlation is computed between ROI residual time series to produce the connectivity matrix. Fisher transformed correlation between ROIs/voxel residual time series is used to produce connectivity maps if that option is selected.

