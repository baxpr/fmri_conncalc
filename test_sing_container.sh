
# jsins
singularity \
run \
--cleanenv --contain \
--bind INPUTS:/INPUTS \
--bind OUTPUTS:/OUTPUTS \
baxpr-fmri_conncalc-master-v1.0.0.simg \
magick_path /usr/bin \
param_file params_JSins.csv \
wroi_file rois_JSins.nii.gz \
roi_file '' \
roiinfo_file rois_JSins.csv \
coregmat_file ../INPUTS/coreg_mat.txt \
deffwd_file ../INPUTS/y_deffwd.nii.gz \
ct1_file ../INPUTS/ct1.nii.gz \
wgm_file ../INPUTS/wgm.nii.gz \
wcseg_file ../INPUTS/wcseg.nii.gz \
func_file ../INPUTS/fmri.nii.gz \
project UNK_PROJ \
subject UNK_SUBJ \
session UNK_SESS \
scan UNK_SCAN \
out_dir /OUTPUTS
