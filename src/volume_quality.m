function [FD,DVARS,badvols] = ...
	volume_quality( ...
	out_dir, ...
	meanfmri_file, ...
	fmri_file, ...
	rp_file, ...
	FDthresh, ...
	DVARSthresh ...
	)


%% Framewise displacement

% We'll calculate the frame-to-frame displacement of points on a sphere of
% radius 50. The matrix at a point gives its displacement relative to zero
% such that
%     Xt = Mt * X0
% We also know
%     Xtm1 = Mtm1 * X0
% We want M such that
%     Xt  = M * Xtm1
% So we get
%     Mt * X0 = M * Mtm1 * X0
%     (Mt*X0) * inv(Mtm1*X0) = M = Mt * X0 * inv(X0) * inv(Mtm1) = Mt/Mtm1

rp = load(rp_file);
nvol = size(rp,1);

% Initialize some variables used for the Power method
amm = nan(nvol,1);
bmm = nan(nvol,1);
gmm = nan(nvol,1);
FD = zeros(nvol,1);

% Examine each time point separately
for t = 2:nvol
	
	% Displacement from previous vol, Power et al 2012 method. This is a
	% sum of displacements over the six degrees of freedom.
	radius = 50;
	fdM = spm_matrix(rp(t,:)) / spm_matrix(rp(t-1,:));
	param = spm_imatrix(fdM);
	qa = spm_matrix([0 0 0 param(4) 0 0]) * [0 radius 0 0]';
	amm(t) = sqrt(sum( (qa-[0 radius 0 0]').^2 ));
	qb = spm_matrix([0 0 0 0 param(5) 0]) * [0 0 radius 0]';
	bmm(t) = sqrt(sum( (qb-[0 0 radius 0]').^2 ));
	qg = spm_matrix([0 0 0 0 0 param(6)]) * [radius 0 0 0]';
	gmm(t) = sqrt(sum( (qg-[radius 0 0 0]').^2 ));
	
	FD(t) = sum(abs( ...
		[rp(t,1)-rp(t-1,1) rp(t,2)-rp(t-1,2) rp(t,3)-rp(t-1,3) ...
		amm(t) bmm(t) gmm(t)] ));
	
end

% Save to file
save(fullfile(out_dir,'FD.txt'),'FD','-ascii');
FD(1) = nan;


%% DVARS

% Use mean fMRI to get brain voxel mask
Vmean = spm_vol(meanfmri_file);
Ymean = spm_read_vols(Vmean);
threshold = spm_antimode(Ymean(:));
mask = Ymean > threshold;

% Frame by frame RMS difference in voxel intensities, percent change units
% relative to the mean intensity of the brain. Compute for each volume
% separately to save loading the whole fmri into memory
V1 = spm_vol([fmri_file ',1']);
Y1 = spm_read_vols(V1);
DVARS = zeros(nvol,1);
for t = 2:nvol
	Vt = spm_vol([fmri_file ',' num2str(t)]);
	spm_check_orientations([V1;Vt]);
	Yt = spm_read_vols(Vt);
	DVARS(t) = 100  ./ mean(Ymean(mask)) .* sqrt(mean( (Yt(mask)-Y1(mask)).^2 ));
	V1 = Vt;
	Y1 = Yt;
end
save(fullfile(out_dir,'DVARS.txt'),'DVARS','-ascii')
DVARS(1) = nan;



%% Compute bad volumes for scrubbing
badvols = (FD>FDthresh) | (DVARS>DVARSthresh);
badvols = badvols | [badvols(2:end); false];
writetable(table(badvols),fullfile(out_dir,'badvols.txt'),'WriteVariableNames',false)


