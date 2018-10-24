function PCTcalc(Vs,Vmu,nsd,ta,fNm)
% function PCTcalc(Vs,Vmu,nsd,ta,fNm)
% Vs     - Residual variance image, mapped volume
% Vmu    - Basline/grand mean image, mapped volume -or- scalar
% nsd    - Normalized standard deviation, var(effect)/sqrt(sigs)
% ta     - t activation threshold
% fNm    - Output filename
%
%
%_____________________________________________________________________________
% @(#)PCTcalc.m	1.4 Thomas Nichols 02/12/03


Vo = struct(    'fname',        fNm,...
                'dim',          [Vs.dim(1:3),spm_type('float')],...
                'mat',          Vs.mat,...
                'pinfo',        [1.0,0,0]',...
                'descrip',      'Percent Change Threshold image');

spm_create_image(Vo);

if ~isstruct(Vmu), mu = Vmu; end

if isstruct(Vmu)
  fprintf('CoefVar Scafact = %g\n',100*ta*nsd);
else
  fprintf('Var Scafact = %g\n',100*ta*nsd/mu);
end


for p=1:Vo.dim(3)

    % Sample volume(s)
    M = spm_matrix([0 0 p]);
    Sigs = spm_slice_vol(Vs,M,Vs.dim(1:2),0);
    if isstruct(Vmu)
       mu = spm_slice_vol(Vmu,M,Vs.dim(1:2),0);
    end

    % Avoid division by zero warnings
    mu(mu==0) = NaN;

    % Compute
    PCT = 100*ta*nsd./mu.*sqrt(Sigs);

    % Write
    Vo = spm_write_plane(Vo,PCT,p);

end        
