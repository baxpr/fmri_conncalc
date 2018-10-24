#!/bin/sh
#
# Compile the matlab code so we can run it without a matlab license. Required:
#     Matlab 2017a, including compiler, with license
#     Installation of SPM12, https://www.fil.ion.ucl.ac.uk/spm/software/spm12/

# Where to find SPM12
spm_dir=/opt/spm12

# For the compiler to find all SPM dependencies, we need to do some stuff from 
# spm_make_standalone.m in our SPM installation. This only needs to be done once
# for a given installation, but it won't make changes if it finds it has already
# run.
matlab -nodesktop -nosplash -sd src -r "prep_spm_for_compile('${spm_dir}'); exit"

# Run the matlab compiler from the linux command line - it's easier to get the 
# correct punctuation.
#
# We need to -I include the root spm directory and the paths that SPM needs - 
# the paths it auto-adds to the matlab path when you run it interactively: 
# config, matlabbatch, matlabbatch/cfg_basicio
#
# We need to -a add SPM's Contents.txt file (produced by prep_spm_for_compile.m)
# so it can find its version number.
#
# We need to -a add non-code image and template files in some SPM directories:
#    canonical, EEGtemplates, toolbox, tpm
#
# Need to -I include some local subdirectories:
#    PCT, canny
#
# Need to -a add some local non-code files that the compiler won't find on its own:
#    pdf_connfilter_figure.fig
#    pdf_connmat_figure.fig
#    params/params_JSins.csv
#    params/JS_insula/rois_JSins.nii.gz
#    params/params_SZhab.csv
mcc -m -v matlab/fmri_conncalc.m \
-I ${spm_dir} \
-I ${spm_dir}/config \
-I ${spm_dir}/matlabbatch \
-I ${spm_dir}/matlabbatch/cfg_basicio \
-a ${spm_dir}/Contents.txt \
-a ${spm_dir}/canonical \
-a ${spm_dir}/EEGtemplates \
-a ${spm_dir}/toolbox \
-a ${spm_dir}/tpm \
-I src/PCT \
-I src/canny \
-a src/pdf_connfilter_figure.fig \
-a src/pdf_connmat_figure.fig \
-a src/params/params_SZhab.csv \
-a src/params/params_JSins.csv \
-a src/params/JS_insula/rois_JSins.nii.gz \
-a src/params/JS_insula/rois_JSins.csv \
-d bin

# Grant lenient execute permissions to the matlab executable and runscript
chmod go+rx bin/fmri_conncalc
chmod go+rx bin/run_fmri_conncalc.sh

