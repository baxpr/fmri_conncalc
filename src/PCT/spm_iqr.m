function [IQR,Qs,n] = spm_iqr(V,M)
% Calculate IQR of image
% [IQR Qs n] = spm_iqr(V[,M])
%
% V     - Mapped image or vector
% M     - Mask; either
%          2-vector of upper & lower exclusive threshold values, or
%              (values outside of (M(1),M(2)) are ignored), or
%          mapped mask image
%
% IQR   - Interquartile range
% Qs    - Quartile
%
%
% Data outside of [M(1) M(2)] are ignored.
%_________________________________________________________________________
% @(#)spm_iqr.m	1.2 T. Nichols 01/07/29

if (nargin<2)
  M  = [];
end

[Qs,n] = spm_prctile_vol(V,M,[25 75]');
IQR = Qs(2) - Qs(1);

