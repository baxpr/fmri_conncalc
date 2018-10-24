function [n, x]=spm_histvol(V, M, nbins)
% Create Histogram of an image volume
% FORMAT [n, x]=spm_histvol(V, M, nbins)
% V     - mapped image volume (see spm_vol)
% M     - Mask; values outside of the mask are ignored and not counted.
%         Can be either either
%          2-vector of upper & lower threshold values 
%              (values outside of (M(1),M(2)) are ignored), or
%          mapped mask image
% nbins - number of bins to use.
% n     - number of counts in each bin
% x     - position of bin centres
%
%
% If nbins is negative, it defines the bin width.
%
% If M is a 2-vector and M(1) is NaN or -Inf, then M(1) is set equal to
% minimum of data; same as for M(2) with the maximum of the data.
%_______________________________________________________________________
% @(#)spm_histvol.m	1.4 John Asburner 01/11/08

if (nargin<2), M = []; end
if (nargin<3), nbins = 256; end

% determine range...
if ~isempty(M) & ~isstruct(M) & all(isfinite(M))
  mn = M(1);
  mx = M(2);
else
  mx = -Inf;
  mn =  Inf;
  for p=1:V.dim(3),
    img = spm_slice_vol(V,spm_matrix([0 0 p]),V.dim(1:2),1);
    if isstruct(M)
      msk = find(isfinite(img) & ...
		 spm_slice_vol(M,spm_matrix([0 0 p]),M.dim(1:2),1)>0);
    else
      msk = find(isfinite(img));
    end
    mx  = max([max(img(msk)) mx]);
    mn  = min([min(img(msk)) mn]);
  end;
end
if ~isempty(M) & ~isstruct(M) & ~all(isfinite(M))
  if isfinite(M(1)), mn = M(1); end
  if isfinite(M(2)), mx = M(2); end
end

% compute histograms...
if nbins>0
  x = [mn:(mx-mn+1)/nbins:mx];
else
  binw = -nbins;
  x = [mn:binw:mx];
  if (mx-x(end))>=binw/2, x = [x x(end)+binw]; end
end
n = zeros(size(x));
for p=1:V.dim(3),
	img = spm_slice_vol(V,spm_matrix([0 0 p]),V.dim(1:2),1);
	if isstruct(M)
	  msk = find(isfinite(img) & ...
		     spm_slice_vol(M,spm_matrix([0 0 p]),M.dim(1:2),1)>0);
	elseif length(M)==2
	  msk = find(isfinite(img) & ...
		     (M(1)<img & img<M(2)));
	else
	  msk = find(isfinite(img));
	end
	n   = n+hist(img(msk),x);
end;
return;
