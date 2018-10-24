function [roi_file,roisnr] = roi_snr(out_dir,swfunc_file,roi_file)

% We compute SNR with reference to the variance of the *global* signal
% rather than region-wise.

% Resample ROI image to fmri geometry and load
flags = struct( ...
	'mask',true, ...
	'mean',false, ...
	'interp',0, ...
	'which',1, ...
	'wrap',[0 0 0], ...
	'prefix','r' ...
	);
spm_reslice({[swfunc_file ',1']; roi_file},flags);
[~,n,e] = fileparts(roi_file);
roi_file = fullfile(out_dir,['r' n e]);

Vroi = spm_vol(roi_file);
Yroi = spm_read_vols(Vroi);

% We'll also make an ROI SNR image
Ysnr = zeros(size(Yroi));

% Load fmri
Vfmri = spm_vol(swfunc_file);
Yfmri = spm_read_vols(Vfmri);
spm_check_orientations([Vroi;Vfmri]);

% Compute mean fMRI, identify in-brain voxels, compute global scaling
% factor, rescale to global mean 1
Ym = mean(Yfmri,4);
t = spm_antimode(Ym(:));
inbrain = Ym(:)>t;
sf = mean(Ym(inbrain));
Yfmri = Yfmri / sf;

% Reshape to vector
Yfmri = reshape(Yfmri,[],size(Yfmri,4))';

% Compute noise std dev
noise = std( mean(Yfmri(:,inbrain),2) );

% Extract ROI values and compute SNR
roi_label = unique(Yroi);
roi_label = roi_label(roi_label~=0);
roisnr = table(roi_label);
for r = 1:height(roisnr)
	inds = Yroi(:)==roi_label(r);
	roidata = Yfmri(:,inds);
	roidata = mean(roidata,2);
	roisnr.SNR(r,1) = mean(roidata) / noise;
	Ysnr(inds) = roisnr.SNR(r,1);
	roisnr.localSNR(r,1) = mean(roidata) / std(roidata);
	roisnr.SSR(r,1) = mean(mean(Yfmri(:,inds)));
end

Vsnr = rmfield(Vroi,'pinfo');
Vsnr.fname = [out_dir '/roi_snr.nii'];
spm_write_vol(Vsnr,Ysnr);

% Save SNR info with labels
roisnr.Properties.VariableNames{'roi_label'} = 'Label';
writetable(roisnr,[out_dir '/roi_info.csv'])

