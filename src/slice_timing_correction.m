function afunc_file = slice_timing_correction(func_file,params)
% Slice timing correction module for connectivity preproc spider. Slice
% timing correction is applied in the third dimension (Z, often axial).

% Number of slices
V = spm_vol(func_file);
nslices = V(1).dim(3);

% Slice order
switch params.slorder
	
	case 'ascending'
		slorder = 1:nslices;
		
	case 'descending'
		slorder = nslices:-1:1;

	case 'ascending_interleaved'
		slorder = [1:2:nslices 2:2:nslices];
		
	case 'descending_interleaved'
		slorder = [nslices:-2:1 (nslices-1):-2:1];

	case 'none'
		warning('Skipping slice timing correction')
		[func_p,func_n,func_e] = fileparts(func_file);
		afunc_file = fullfile(func_p,['a' func_n func_e]);
		copyfile(func_file,afunc_file);
		return
		
	otherwise
		error('Unknown slice order')
		
end


% Call the SPM routine
tr = params.tr;
ta = params.tr - params.tr/nslices;
spm_slice_timing(func_vols,slorder,1,[ta/nslices-1 tr-ta],'a');


% Filename for slice time corrected images. The file doesn't
% exist, so we have to predict its name
clear afunc_filename
[func_p,func_n,func_e] = fileparts(func_file);
afunc_file = fullfile(func_p,['a' func_n func_e]);

