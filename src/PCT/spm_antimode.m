function [m,mm]= spm_antimode(V,M,qT)
% Calculate least frequent value in image
% [m,mm] = spm_antimode(V,M)
%
% V   - Mapped image or vector
% M   - Mask, used to eliminate bad/background regions
%          2-vector of upper & lower excluseive thresholds,
%              (values outside of (M(1),M(2)) are ignored), or
%          mapped mask image
% qT  - Quantile thresholds; applied after masking, used to eliminate
%       low-frequency tails of distribution of image values.
% 
% 
% m   - Estimate of location of minimum density.
% mm  - Estimate of location of maximum density (above m)
%
% qT defaults to [10 90], that is, by default values below 10%ile and
% above 90%ile are ignored.
%
%
% Based on the following result: If f is a density continuous on interval
% [a b], the location of the largest interval between observations is a
% consistent estimator of minimum of f on [a b] (it converges in
% probability).
%
% This method will fail, however, if image consists soley of ingegral
% values.  If this appears to happen, we then approximate the location of
% the largest inter-observation interval with the minimal histogram bin. 
%
% Since it is easily available, the mode of image intensities greater
% than the antimode is also returned.
%
%
%
% Reference: JA Hartigan, 1977, "Distribution Problems in Clustering,"
% in "Classification and Clustering", ed JV Ryzin, Academic Press, NY, pp
% 63-65.
%
%_________________________________________________________________________
% @(#)spm_antimode.m	1.10 T. Nichols 02/11/29

if nargin<2, M  = []; end
if nargin<3 | isempty(qT), qT = [10 90]; end
FindMode = (nargout>1);


%
% Load, mask, and massage image
%
if isstruct(V)
  Img  = spm_read_vols(V);
  Img  = Img(:);
else
  Img  = V(:);
end

% Load mask
if isstruct(M), 
  ImgMsk = spm_read_vols(M); 
  ImgMsk = ImgMsk(:);
end

%
% Apply mask and threshold
%
T = [];
if isempty(qT)
  if isstruct(M)
    % Apply mask 
    Img(~ImgMsk) = [];
  end
else
  if all(qT==[0 100])
    T  = [-Inf Inf];
  else
    T  = spm_prctile_vol(V,M,qT);
  end
  if FindMode
    % Save upper tail for later
    if isstruct(M)
      uTail = Img(T(2)<=Img & ImgMsk);
    else
      uTail = Img(T(2)<=Img);
    end      
  end
  % Threshold (and maybe mask)
  if isstruct(M)
    Img(Img<=T(1) | T(2)<=Img | ~ImgMsk) = [];
  else
    Img(Img<=T(1) | T(2)<=Img) = [];
  end
end
Img(~isfinite(Img)) = [];


%
% Hartigan approach
%
sImg  = sort(Img(:));
dsImg = diff(sImg);

% Carefully find location of order statistic pairs
dsIi  = find(dsImg==max(dsImg));
dsIip = min(dsIi+1,length(sImg));          %-Just in case hit end
ms    = mean([sImg(dsIi)';sImg(dsIip)']);  %-Middle of pair
m     = mean(ms);


if length(ms)>5

  %
  % If there were more than, say, 5 modes, we've probably got a discrete
  % image.  
  %
  % Use histogram approach 
  %

  % Find hstogram bin width
  [IQR,Qs,n] = spm_iqr(V,T);
  % Rule from Scott (1992), pg 100; 1.595 = 2.15/1.348
  binw       = 1.595*IQR/n^(1/5);

  nbin = ceil((max(Img)-min(Img))/binw);

  [n,x] = hist(Img,nbin);
  
  % Drop first and last histogram bins, just in case there are problems
  % with truncated images
  n([1 end]) = [];
  x([1 end]) = [];

  m     = mean(x(find(n==min(n))));

end


if (nargout>1)
  %
  % Find mode greater than antimode
  %
  
  % Toss out values less than antimode
  Img(Img<m) = [];

  % Add in upper tail
  Img = [Img; uTail(isfinite(uTail))];

  % Find hstogram bin width again
  [IQR,Qs,n] = spm_iqr(Img);
  % Rule from Scott
  binw       = 1.595*IQR/n^(1/5);

  nbin  = ceil((max(Img)-min(Img))/binw);

  [n,x] = hist(Img,nbin);
  
  mm    = mean(x(find(n==max(n))));

end
