function [coregmat_file,deffwd_file,ct1_file,wgm_file,wcseg_file,func_file,roi_file] = ...
	prep_files( ...
	out_dir,coregmat_file,deffwd_file,ct1_file,wgm_file,wcseg_file,func_file,roi_file)

copyfile(coregmat_file,[out_dir '/coregmat.txt']);
coregmat_file = [out_dir '/coregmat.txt'];

copyfile(wgm_file,[out_dir '/wgm.nii.gz']);
system(['gunzip -f ' out_dir '/wgm.nii.gz']);
wgm_file = [out_dir '/wgm.nii'];

copyfile(ct1_file,[out_dir '/ct1.nii.gz']);
system(['gunzip -f ' out_dir '/ct1.nii.gz']);
ct1_file = [out_dir '/ct1.nii'];

copyfile(wcseg_file,[out_dir '/wcseg.nii.gz']);
system(['gunzip -f ' out_dir '/wcseg.nii.gz']);
wcseg_file = [out_dir '/wcseg.nii'];

copyfile(deffwd_file,[out_dir '/y_deffwd.nii.gz']);
system(['gunzip -f ' out_dir '/y_deffwd.nii.gz']);
deffwd_file = [out_dir '/y_deffwd.nii'];

copyfile(func_file,[out_dir '/func.nii.gz']);
system(['gunzip -f ' out_dir '/func.nii.gz']);
func_file = [out_dir '/func.nii'];

copyfile(roi_file,[out_dir '/roi_labels.nii.gz']);
system(['gunzip -f ' out_dir '/roi_labels.nii.gz']);
roi_file = [out_dir '/roi_labels.nii'];

