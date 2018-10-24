function dfunc_file = drop_volumes(func_file,params)

% Drop initial volumes to account for saturation effects. Assumes 4D Nifti.

% Load images
V = spm_vol(func_file);
Y = spm_read_vols(V);

% If we don't have a drop count, make it zero
if isempty(params.drop_initialvols)
	warning('No dropvols - assuming zero')
	params.drop_initialvols = 0;
end

% If we don't have a vol limit, use them all
if isempty(params.vols_kept)
	warning('No vols_kept - keeping all')
	params.vols_kept = size(Y,4) - params.drop_initialvols;
end	

% Check that things make sense
if (params.drop_initialvols + params.vols_kept) > size(Y,4)
	error('More volumes requested than exist in %s',func_file)
end

% Drop volumes
keeps = (1:params.vols_kept) + params.drop_initialvols;
outV = V(keeps);
outY = Y(:,:,:,keeps);

% Output filename and updated indices, and write
[func_p,func_n,func_e] = fileparts(func_file);
dfunc_file = fullfile(func_p,['d' func_n func_e]);
for v = 1:length(outV)
	outV(v).fname = dfunc_file;
	outV(v).n(1) = v;
	spm_write_vol(outV(v),outY(:,:,:,v));
end
