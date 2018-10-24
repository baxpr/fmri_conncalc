function PCT(SPM,c,ta)
% Create percent change threshold image
% function PCT
% 
% Sub function
% function PCT(SPM,c,ta)
% SPM    - Pathname of SPM.mat; image put in this directory by default
% c      - Contrast to define percent change image of interest
% ta     - t activation threshold
%
%
% It is assumed that c has been scaled such that c*beta represents a unit
% change effect.
%
%____________________________________________________________________________
% @(#)PCT.m	1.14 T. Nichols 02/12/19
%

if nargin==0

  pos = 1;

  %-Select SPM.mat
  %-----------------------------------------------------------------------
  SPM   = spm_get(1,'SPM.mat','Select SPM.mat');
  swd   = spm_str_manip(SPM,'H');
  xSDM  = load(SPM);
  if exist(fullfile(swd,'xCon.mat'),'file')
    load(fullfile(swd,'xCon.mat'));
  else
    xCon = [];
  end

  %-Get contrast
  %-------------------------------------------------------------------
  c = [];
  if ~isempty(xCon)
    Ts = find([xCon(:).STAT] == 'T');
    if length(Ts)
      C    = [xCon(Ts).c]';
      nVa  = diag(C*xSDM.xX.Bcov*C');       % Normalized variance
      I    = Ts(min(find(max(nVa)==nVa)));  % Take largest (first if multiple)
      c    = xCon(I).c';

      str = [sprintf('%3g ',c) sprintf('\n')];
      spm('alert',str,sprintf('Default Contrast (%d)',I),1,0);
    end
  end

  if length(c)==0 | ...
	spm_input(sprintf('Use default (%d) contrast?',I),...
		  pos,'y/n',[],1) == 'n'

    [I,xCon] = spm_conman(xSDM.xX,xCon,'T',1,'Enter contrast for PCT','',1);
    
    %-Save contrast structure  - ensures new contrasts are saved
    save(fullfile(swd,'xCon.mat'),'xCon');
    
    MyxCon = xCon(I);
    c      = MyxCon.c';
  else
    MyxCon = xCon(I);
  end

  % Enforce contrast scaling
  cPos =  sum(c(c>0));  
  cNeg = -sum(c(c<0));
  c = c/max([cPos cNeg]);
  if (cPos>0) & (cNeg>0) & (abs(cPos-cNeg)>1e-3)
    warning('Differencing contrast that doesn''t sum to zero!\nHope that''s OK')
  end

  %-Check that design matrix is appropriately scaled?
  % Useful for fMRI designs created before my modifications.
  % Won't work for event related designs with closely spaced events.
  %-------------------------------------------------------------------
  if  spm_input('Enforce scaling of design matrix?','+1','y/n',[],1)=='y'

    % Find D:
    % D is the scaling required s.t. X*diag(D)'s columns give the beta's
    % the interpretation of unit change in the data.
    D       = ones(size(c));
    D(c~=0) = 1./(max(xSDM.xX.X(:,c~=0))-min(xSDM.xX.X(:,c~=0)));

    % Instead of scaling X by D, can just divide c by D
    c       = c./D;

    if exist('xCon')
      xCon(I).c = c';
      try, 
	save(fullfile(swd,'xCon.mat'),'xCon');
	spm('alert','Rescaled contrast saved','Contrast rescaling',1)
      catch
	spm('alert','Rescaled contrast *not* saved','Contrast rescaling',1)
      end
    end

  end


  %-Get height threshold
  %-------------------------------------------------------------------
  n     = 1;
  STAT  = 'T';
  edf   = [1 xSDM.xX.erdf];
 
  if isfield(MyxCon,'Vspm') & ~isempty(MyxCon.Vspm)
    % have a statistic image
    str   = 'FWE|FDR|uncorrected';
  else
    str   = 'FWE|uncorrected'; 
  end
  d = spm_input('corrected height threshold','+1','b',str,[],3);
  if strcmp(d,'FWE')
  	u  = spm_input('FWE-corrected p value','+0','r',0.05,1,[0,1]);
  	u  = spm_uc(u,edf,STAT,xSDM.R,n,xSDM.S);
  elseif strcmp(d,'FDR')
	VspmSv = spm_vol(fullfile(swd,MyxCon.Vspm));
  	u  = spm_input('FDR-corrected p value','+0','r',0.05,1,[0,1]);
  	u  = spm_uc_FDR(u,edf,STAT,n,VspmSv,0);
  else
  	%-NB: Uncorrected p for conjunctions is p of the conjunction SPM
  	u  = spm_input(['threshold {',STAT,' or p value}'],'+0','r',0.05,1);
  	if u <= 1; u = spm_u(u^(1/n),edf,STAT); end
  end

  str = ['c = ' sprintf('%3g ',c) sprintf('\n')];
  str = [str sprintf('\tu = %g',u)];
  spm('alert',str,'PCT parameters',1,0);

  PCT(SPM,c,u);

else
  % We assume c is appropriately scaled


  PCTnm = 'ResRMS_PCT';

  swd   = spm_str_manip(SPM,'H');
  c     = c(:)';

  load(SPM,'xX');

  % No special treatment of fMRI is needed, since xX.Bcov accounts for
  % filtering.
  nsd   = sqrt(c*xX.Bcov*c');
  

  PMn  = fullfile(swd,'Mean.img');
  Vs   = spm_vol(fullfile(swd,'ResMS.img'));
  if ~exist(PMn,'file')
    VMn = spm_grand_mean(SPM,PMn);
  else
    VMn = spm_vol(PMn);
  end
  % If more than 10% of voxels null/NaN, don't bother with
  % antimode... background has already been cropped.
  d = spm_read_vols(VMn); dn=sum(isnan(d(:)))+sum(d(:)==0);
  if (dn > 0.1*length(d(:)))
    Mn = spm_mode(VMn);
  else
    Mn = spm_mode(VMn,[],[]);
  end
  clear d;
  spm('alert',num2str(Mn),'Modal/Global of Mean',1)
  
  fprintf('Mode = %g\n',Mn);

  PCTcalc(Vs,Mn,nsd,ta,fullfile(swd,'ResRMS_PCT'))
  PCTcalc(Vs,VMn,nsd,ta,fullfile(swd,'PCT'))

end

