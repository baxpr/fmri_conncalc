%% Generate insula ROIs.

% Erosion of 1 voxel from the edges of each sub-region mask to ensure they
% are not touching

% Requires SPM12

roi_files = {
	'InsulaCluster_K3_L-dAI.nii.gz'
	'InsulaCluster_K3_L-PI.nii.gz'
	'InsulaCluster_K3_L-vAI.nii.gz'
	'InsulaCluster_K3_R-dAI.nii.gz'
	'InsulaCluster_K3_R-PI.nii.gz'
	'InsulaCluster_K3_R-vAI.nii.gz'
	};

% Set up the ROI label table
rois = table(roi_files,'VariableNames',{'file'});
rois.label = (1:height(rois))';
for h = 1:height(rois)
	q = strrep(rois.file{h},'InsulaCluster_K3_','');
	q = strrep(q,'.nii.gz','');
	rois.name(h) = {q};
end

% Load ROI masks
for h = 1:height(rois)
	system(['gunzip -fk ' rois.file{h}]);
	V = spm_vol(rois.file{h}(1:end-3));
	if h==1
		Y = spm_read_vols(V);
	else
		Y(:,:,:,h) = spm_read_vols(V);
	end
	system(['rm -f ' rois.file{h}(1:end-3)]);
end

% Check for ROI overlap
if any(any(any(sum(Y,4)>1)))
	warning('ROI overlap')
end

% Single-voxel erosion neighborhood
nhood = zeros(3,3,3);
nhood(:,:,1) = [0 0 0; 0 1 0; 0 0 0];
nhood(:,:,2) = [0 1 0; 1 1 1; 0 1 0];
nhood(:,:,3) = [0 0 0; 0 1 0; 0 0 0];
eY = zeros(size(Y));
for h = 1:height(rois)
	eY(:,:,:,h) = imerode(Y(:,:,:,h),nhood);
end

% Collapse into single volume
s = size(Y);
Yall = zeros(s(1:3));
eYall = zeros(s(1:3));
for h = 1:height(rois)
	q = Y(:,:,:,h);
	Yall(q==1) = h;
	q = eY(:,:,:,h);
	eYall(q==1) = h;
end

% Write images to file
Vout = V;
Vout.pinfo(1:2) = [0 1];
Vout.fname = 'rois_JSins.nii';
spm_write_vol(Vout,Yall);
system(['gzip -f ' Vout.fname]);
Vout.fname = 'eroded_rois_JSins.nii';
spm_write_vol(Vout,eYall);
system(['gzip -f ' Vout.fname]);

% Save ROI info
writetable(rois,'rois_JSins.csv');
