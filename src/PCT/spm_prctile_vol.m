function [x,n] = spm_prctile_vol(V,M,p)
% spm_prctile_vol: Calculate percentiles of image 
% FORMAT [x,n] = spm_prctile_vol(V,M,p)
% Input: 
% V     - Mapped image or vector
% M     - Mask; either2-vector of upper & lower exclusive threshold
%         values, (values outside of (M(1),M(2)) are ignored), or mapped
%         mask image (not currently supported) 
% p     - Vector of percentiles, (in [0,100]).
% Output:
% x     - Values of percentiles
% n     - Number of image voxels used. 
%_________________________________________________________________________
% @(#)spm_prctile_vol.m	1.2 02/11/29
% Based on spmd_prctile_vol.m 1.7 T. Nichols 02/10/03

%-----------------------------Function called ---------------------------
% spm_read_vols.m
%------------------------------------------------------------------------

if isstruct(V)
  Img = spm_read_vols(V);
  Img = Img(:);
else
  Img = V(:);
end

if isstruct(M)
  if ~isstruct(V) & length(V(:))~=prod(M.dim(1:3)),
    error('Length of vector V doesn''t match mask')
  end
  Img(~spm_read_vols(M)) = [];
elseif ~isempty(M)
  Img(Img<=M(1) | M(2)<=Img) = [];
end

Img(~isfinite(Img)) = [];

x  = spm_prctile(Img,p(:));
n  = length(Img);
