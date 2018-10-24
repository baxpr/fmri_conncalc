function [filtered_fmri_file,roidata_file,connmat_file] = ...
	connectivity_filter( ...
	out_dir, ...
	badvols, ...
	rp_file, ...
	wcseg_file, ...
	wfmri_file, ...
	roi_file, ...
	params, ...
	remove_gm, ...
	FD, ...
	DVARS, ...
	project, ...
	subject, ...
	session, ...
	scan ...
	)

% NOTE: tr is given in sec.
%
% Bandpass filter is implemented as a set of sine and cosine basis
% functions in the confound matrix.

%% Filename tag based on filtering params
if remove_gm==0 && isempty(badvols)
	filetag = 'keepgm_noscrub';
elseif remove_gm==0 && ~isempty(badvols)
	filetag = 'keepgm_scrub';
elseif remove_gm==1 && isempty(badvols)
	filetag = 'removegm_noscrub';
elseif remove_gm==1 && ~isempty(badvols)
	filetag = 'removegm_scrub';
end	


%% Motion and first differences regressors
rp = load(char(rp_file));
if params.mot_PCs == 0
	mot_regr = [];
elseif params.mot_PCs == 6
	mot_regr = rp;
else
	[~,P] = pca(zscore(rp));
	mot_regr = P(:,1:params.mot_PCs);
end

rp_delta = [zeros(1,6); diff(rp)];
if params.motderiv_PCs == 0
	motderiv_regr = [];
elseif params.motderiv_PCs == 6
	motderiv_regr = rp_delta;
else
	[~,P] = pca(zscore(rp_delta));
	motderiv_regr = P(:,1:params.motderiv_PCs);
end


%% Resample the w seg image to the fmri voxel geometry
flags = struct( ...
	'mask',true, ...
	'mean',false, ...
	'interp',0, ...
	'which',1, ...
	'wrap',[0 0 0], ...
	'prefix','r' ...
	);
spm_reslice({[wfmri_file ',1']; wcseg_file},flags);
[p,n,e] = fileparts(wcseg_file);
rwcseg_file = fullfile(p,['r' n e]);


%% Gray, white, CSF regressors from unsmoothed fmri

% Read images and reshape
wcsegV = spm_vol(char(rwcseg_file));
wcsegY = spm_read_vols(wcsegV);

fmriV = spm_vol(char(wfmri_file));
fmriY = spm_read_vols(fmriV);
fmriY = reshape(fmriY,[],size(fmriY,4))';

% Define GM and WMCSF compartments
wmcsfY = ismember(wcsegY,[40 41 44 45   4 11 49 50 51 52]);
gmY = (wcsegY>0) & ~wmcsfY;

% Erode WMCSF
nhood = nan(3,3,3);
nhood(:,:,1) = [0 0 0; 0 1 0; 0 0 0];
nhood(:,:,2) = [0 1 0; 1 1 1; 0 1 0];
nhood(:,:,3) = [0 0 0; 0 1 0; 0 0 0];
erodeY = wmcsfY;
erodeY = imerode(erodeY,nhood);
erodeV = wcsegV;
erodeV.pinfo(1:2) = [1 0];
erodeV.fname = [out_dir '/wmcsf_eroded.nii'];
spm_write_vol(erodeV,erodeY);

% Gray matter regressor
if remove_gm
	gray_signals = fmriY(:,gmY(:)>0);
	gray_regr = mean(gray_signals,2);
else
	gray_regr = [];
end

% WMCSF regressors 
wmcsf_signals = fmriY(:,erodeY(:)>0);
if params.wmcsf_PCs == 0
	wmcsf_regr = [];
elseif params.wmcsf_PCs == 1
	wmcsf_regr = mean(wmcsf_signals,2);
else
	[~,wmcsf_PCs] = pca(zscore(wmcsf_signals));
	wmcsf_regr = wmcsf_PCs(:,1:params.wmcsf_PCs);
end

% Save memory
clear wcsegY fmriY


%% Bandpass filter regressors
bp_regr = fourier_filter_basis( ...
	size(rp,1),params.tr,params.bandpasslo_hz,params.bandpasshi_hz);


%% Scrub regressors
badvol_regr = zeros(length(badvols),sum(badvols));
badinds = find(badvols);
for b = 1:length(badinds)
    badvol_regr(badinds(b),b) = 1;
end


%% Create and apply confound (filter) matrix

% Read unsmoothed images
fmriV = spm_vol(char(wfmri_file));
fmriY = spm_read_vols(fmriV);
o = size(fmriY);
fmriY = reshape(fmriY,[],o(4))';

% Scale image data to percent of global mean
meanfmri = mean(fmriY,1);
thresh = spm_antimode(meanfmri);
globalmean = mean(meanfmri(meanfmri>thresh));
fmriY = 100 * fmriY / globalmean;

% Build confound matrix
confounds = [ ...
    zscore(bp_regr) ...
    zscore(mot_regr) ...
    zscore(motderiv_regr) ...
	badvol_regr ...
    zscore(wmcsf_regr) ...
    zscore(gray_regr) ...
    ];
confound_matrix_file = fullfile(char(out_dir),['confounds_' filetag '.txt']);
save(confound_matrix_file,'confounds','-ascii')

% Regress out the confounds from the images
desmtx = [confounds ones(size(confounds,1),1)];
beta = lscov(desmtx, fmriY);
fmriYc = fmriY - desmtx * beta;

% In-brain R values for random sample of grey matter voxels with adequate
% fMRI signal
keeps = find(gmY(:)>0 & meanfmri(:)>thresh);
randkeeps = keeps(randperm(length(keeps),2500));
Rpre = ( corr(fmriY(:,randkeeps)) );
Rpost = ( corr(fmriYc(:,randkeeps)) );

% Save memory
clear fmriY


% Extract ROI signals
roiV = spm_vol(roi_file);
roiY = spm_read_vols(roiV);

roi_label = unique(roiY);
roi_label = roi_label(roi_label~=0);

roi_data = nan(size(fmriYc,1),length(roi_label));

for r = 1:length(roi_label)
	inds = roiY(:)==roi_label(r);
	tmp = fmriYc(:,inds);
	roi_data(:,r) = mean(tmp,2);
end

% Compute connectivity Z maps
if params.output_maps
	
	Z = atanh(corr(roi_data,fmriYc)) * sqrt(size(roi_data,1)-3);
	
	% Write out connectivity maps
	Z = reshape(Z',[o(1:3) size(roi_data,2)]);
	filtered_Z_file = fullfile(out_dir,['Zmap_' filetag '.nii']);
	for v = 1:size(roi_data,2)
		thisV = rmfield(fmriV(1),'pinfo');
		thisV.dt(1) = spm_type('float32');
		thisV.n(1) = v;
		thisV.fname = filtered_Z_file;
		spm_write_vol(thisV,Z(:,:,:,v));
	end
	
	% Smooth Z maps
	smooth(filtered_Z_file,params);
	
end

% Write out the filtered unsmoothed images
fmriYc = reshape(fmriYc',o);
filtered_fmri_file = fullfile(out_dir,['filtered_' filetag '.nii']);
for v = 1:numel(fmriV)
    thisV = rmfield(fmriV(v),'pinfo');
    thisV.dt(1) = spm_type('float32');
    thisV.fname = filtered_fmri_file;
    spm_write_vol(thisV,fmriYc(:,:,:,v));
end

% Save ROI time series to file, marked by the image label
roi_label_str = strcat('r', strtrim(cellstr(num2str(roi_label))') );
roi_data = array2table(roi_data,'VariableNames',roi_label_str);
roidata_file = [out_dir '/roi_timeseries_' filetag '.csv'];
writetable(roi_data,roidata_file);

% Save memory
clear fmriYc Z


%% Compute connectivity R
% Mark rows and columns by ROI image label
connmat = corr(table2array(roi_data));
connmat = array2table(connmat,'VariableNames',roi_label_str, ...
	'RowNames',roi_label_str);
connmat_file = [out_dir '/connectivity_matrix_R_' filetag '.csv'];
writetable(connmat,connmat_file,'WriteRowNames',true)


%% Save a stats file
statsstr = sprintf( ...
    [ ...
    'FD_median=%0.2f\n' ...
    'DVARS_median=%0.2f\n' ...
    'bandpass_lo=%0.3f Hz\n' ...
    'bandpass_hi=%0.3f Hz\n' ...
    'TR=%0.3f\n' ...
    'mot_PCs=%d\n' ...
    'motderiv_PCs=%d\n' ...
    'wmcsf_PCs=%d\n' ...
    'remove_gray=%d\n' ...
    'slorder=%s\n' ...
    'init_dropvols=%d\n' ...
	'kept_vols=%d\n' ...
	'scrubbed_vols=%d\n' ...
	'Final_DOF=%d\n' ...
    ], ...
    nanmedian(FD), nanmedian(DVARS), ...
    params.bandpasslo_hz,params.bandpasshi_hz, params.tr, ...
    params.mot_PCs, params.motderiv_PCs, ...
	params.wmcsf_PCs, remove_gm, ...
	params.slorder, params.drop_initialvols, params.vols_kept, ...
	sum(badvols), size(desmtx,1)-size(desmtx,2) );
f = fopen([out_dir '/stats_' filetag '.txt'],'wt');
fprintf(f,statsstr);
fclose(f);


%% We also generate a report for PDF:

% Figure out screen size so the figure will fit
ss = get(0,'screensize');
ssw = ss(3);
ssh = ss(4);
ratio = 8.5/11;
if ssw/ssh >= ratio
	dh = ssh;
	dw = ssh * ratio;
else
	dw = ssw;
	dh = ssw / ratio;
end

% Create figure
pdf_figure = openfig('pdf_connfilter_figure.fig','new');
set(pdf_figure,'Tag','pdf_connfilter');
set(pdf_figure,'Units','pixels','Position',[0 0 dw dh]);
figH = guihandles(pdf_figure);

% Summary
set(figH.summary_text, 'String', sprintf( ...
    [ ...
    'FD_median %0.2f\n' ...
    'DVARS_median %0.2f\n' ...
    'bandpass_lo %0.3f Hz\n' ...
    'bandpass_hi %0.3f Hz\n' ...
    'TR %0.3f s\n' ...
    'mot_PCs %d\n' ...
    'motderiv_PCs %d\n' ...
    'wmcsf_PCs %d\n' ...
    'remove_gray %d\n' ...
    'slorder %s\n' ...
    'init_dropvols %d\n' ...
	'kept_vols %d\n' ...
	'scrubbed_vols %d\n' ...
	'Final DOF %d\n' ...
    ], ...
    nanmedian(FD), nanmedian(DVARS), ...
    params.bandpasslo_hz,params.bandpasshi_hz, params.tr, ...
    params.mot_PCs, params.motderiv_PCs, ...
	params.wmcsf_PCs, remove_gm, ...
	params.slorder, params.drop_initialvols, params.vols_kept, ...
	sum(badvols), size(desmtx,1)-size(desmtx,2) ))

% Scan info
set(figH.scan_info, 'String', sprintf( ...
	'%s: %s, %s, %s, %s', ...
	filetag, project, subject, session, scan));
set(figH.date,'String',['Report date: ' date]);
set(figH.version,'String',['Matlab version: ' version]);


% Compute COM of image in each axis
fmriV = spm_vol([char(wfmri_file) ',1']);
fmriY = spm_read_vols(fmriV);
q3 = squeeze(nansum(nansum(fmriY,1),2));
fmricom3 = round(sum((1:length(q3))' .* q3) / sum(q3));

% Image slice
axes(figH.slice)
img = fmriY(:,:,fmricom3);
imagesc(imrotate(abs(img),90))
colormap(gray)
axis image
axis off
title('Raw image')

% GM mask
axes(figH.gm)
img = gmY(:,:,fmricom3);
imagesc(imrotate(abs(img),90))
colormap(gray)
axis image
axis off
title('Grey matter')

% WM/CSF mask
axes(figH.wmcsf)
img = erodeY(:,:,fmricom3);
imagesc(imrotate(abs(img),90))
colormap(gray)
axis image
axis off
title('WM + CSF eroded')

% FD
axes(figH.FD)
plot(FD,'b')
hold on
plot(find(badvols),FD(badvols),'ro')
title('Framewise Displacement (FD)')
set(gca,'XTick',[],'XLim',[0 length(FD)+1])

% DVARS
axes(figH.DVARS)
plot(DVARS,'b')
hold on
plot(find(badvols),DVARS(badvols),'ro')
title('Signal Change (DVARS)')
xlabel('Volume')
set(gca,'XLim',[0 length(DVARS)+1])

% Histogram plot
bins = -1:0.01:1;
Hpre = hist(Rpre(:),bins);
Hpost = hist(Rpost(:),bins);
axes(figH.histograms)
plot(bins,Hpre/length(Rpre(:)),'r')
hold on
plot(bins,Hpost/length(Rpost(:)),'b')
plot([0 0],get(gca,'Ylim'),':k')
legend({'Unfiltered' 'Filtered'},'Interpreter','None')
xlabel('Voxel-voxel correlation (R)')
ylabel('Frequency in random gray matter sample')

% Print to PNG
print(gcf,'-dpng',fullfile(out_dir,['connectivity_' filetag '.png']))
close(gcf);


%% Second page of PDF with the connectivity matrix

% Create figure
pdf_figure = openfig('pdf_connmat_figure.fig','new');
set(pdf_figure,'Tag','pdf_connmat');
set(pdf_figure,'Units','pixels','Position',[0 0 dw dh]);
figH = guihandles(pdf_figure);

% Summary
set(figH.summary_text, 'String', sprintf( ...
    'remove_gray %d; scrubbed_vols %d\n', ...
	remove_gm, sum(badvols) ))

% Scan info
set(figH.scan_info, 'String', sprintf( ...
	'%s: %s, %s, %s, %s', ...
	filetag, project, subject, session, scan));
set(figH.date,'String',['Report date: ' date]);
set(figH.version,'String',['Matlab version: ' version]);

% Conn matrix
axes(figH.connmatrix)
imagesc(table2array(connmat),[-1 1])
colormap(jet)
axis image
set(gca,'XTickLabelRotation',90)

% Print to PNG
print(gcf,'-dpng',fullfile(out_dir,['connmatrix_' filetag '.png']))
close(gcf);


