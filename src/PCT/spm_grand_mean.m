function V = spm_grand_mean(SPM,MnNm)
% Creates grand mean using beta images
% V = GrandMean(SPMdir,MnNm)
%
% SPM    - Pathname of SPM.mat
% MnNm   - Mean filename
%
% V      - Mapped mean image
%
%
% Ideally grand mean would be consistently parameterized as a particular
% beta (e.g. beta_p).  Since we can't be sure of this, we just work out
%
%   1/n sum(Y) = 1/n sum(Yh) = 1/n 1'Yh = 1/n 1'X betah = 1/n (1'X) betah
%
%________________________________________________________________________
% @(#)spm_grand_mean.m	1.4 T. Nichols 02/11/17

load(SPM,'xX','Vbeta');

X = xX.X;
n = size(X,1);
p = size(X,2);

Vi = spm_vol(str2mat(Vbeta{:}));
Vo = struct(    'fname',        MnNm,...
                'dim',          [Vi(1).dim(1:3),spm_type('float')],...
                'mat',          Vi(1).mat,...
                'pinfo',        [1.0,0,0]',...
                'descrip',      'spm - grand mean image');
%-Create basic header
Vo = spm_create_image(Vo);

oXb = ones(n,1)'*X;

% Set scalefactors to perform weighted sum (see spm_add).
for i = 1:p
  Vi(i).pinfo(1:2,:)  = Vi(i).pinfo(1:2,:)*oXb(i)/n;
end

Vo.pinfo(1,1) = spm_add(Vi,Vo);

%-Write header (update with scaling information)
Vo = spm_create_image(Vo);

if (nargout)
  V = Vo;
end
