function wimage_file = warp_images( ...
	fwddef_file, ...
	image_file, ...
	ref_image, ...
	interp, ...
	out_dir ...
	)

clear job
job.comp{1}.def = {fwddef_file};
job.comp{2}.id.space = {ref_image};
job.out{1}.pull.fnames = {image_file};
job.out{1}.pull.savedir.saveusr = {out_dir};
job.out{1}.pull.interp = interp;
job.out{1}.pull.mask = 0;
job.out{1}.pull.fwhm = [0 0 0];

[~,n,e] = fileparts(image_file);
wimage_file = fullfile(out_dir,['w' n e]);

spm_deformations(job);

