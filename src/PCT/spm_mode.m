function Mo = spm_mode(V,M,qT)
% Calculate modal value in image
% M = spm_mode(V,Vm,qT)
%
% V   - Mapped image
% M   - Mask; either
%          2-vector of upper & lower threshold values,
%              (values outside of (M(1),M(2)) are ignored), or
%          mapped mask image
% qT  - Quantile thresholds to use with spm_antimode to determine lower
%       threshold; can be [], indicating use of spm_antimode with defaults.
%       Useful for segmenting brain from nonbrain.
% 
% 
%
% Mode estimate is determined with a histogram using bin widths of
% 1.595*IQR/n^(1/5), as per Scott (1992), pg 100.
%
%_________________________________________________________________________
% @(#)spm_mode.m	1.7 T. Nichols 01/12/07

if nargin<2, 
  M = [];
end
if (nargin>2)
  [m,Mo] = spm_antimode(V,M,qT);
  return
end

%
% Histogram/frequency polygon bin width
% 
% Scott (1992), pg 100; 1.595 = 2.15/1.348
%
[IQR,Qs,n] = spm_iqr(V,M);
binw       = 1.595*IQR/n^(1/5);

[n,x] = spm_histvol(V,M,-binw);

% Mean taken in case of nonunique local modes
Mo = mean(x(n==max(n)));
