function xY = create_voi(SPM,option,P)

% code from spm_regions.m

xY.name  = option.name;

%-Get raw data, whiten and filter 
%--------------------------------------------------------------------------
V        = hdr_read(SPM.xY.VY,P.fmripath);
y        = spm_get_data(V,option.vox');
y        = spm_filter(SPM.xX.K,SPM.xX.W*y);
 
 
%-Computation
%==========================================================================
 
%-Remove null space of contrast
%-------------------------------------------------------------------------- 
xY.Ic = option.adjust;
cd(P.glmpath)
if xY.Ic ~= 0
 
    %-Parameter estimates: beta = xX.pKX*xX.K*y
    %----------------------------------------------------------------------
    V        = hdr_read(SPM.Vbeta,P.glmpath);
    beta     = spm_get_data(V,option.vox');
 
    %-subtract Y0 = XO*beta,  Y = Yc + Y0 + e
    %----------------------------------------------------------------------
    if ~isnan(xY.Ic)
        y = y - spm_FcUtil('Y0',SPM.xCon(xY.Ic),SPM.xX.xKXs,beta);
    else
        y = y - SPM.xX.xKXs.X * beta;
    end
 
end
 
%-Confounds
%--------------------------------------------------------------------------
xY.X0     = SPM.xX.xKXs.X(:,[SPM.xX.iB SPM.xX.iG]);
 

% and add session-specific filter confounds and motion
%--------------------------------------------------------------------------
Xf = [];
for s=1:numel(SPM.nscan)
    Xf      = blkdiag(Xf,SPM.xX.K(s).X0);
    
end


Mc = [];
for s=1:numel(SPM.nscan)
    for mo = 1:6
        idmo = find(strcmp(SPM.xX.name,sprintf('Sn(%d) R%d',s,mo)));
        Mc = [Mc,SPM.xX.xKXs.X(:,idmo)];
    end
end


xY.X0 = [xY.X0 Xf Mc];

 

%-Remove null space of X0
%--------------------------------------------------------------------------
xY.X0     = xY.X0(:,any(xY.X0));
 
 
%-Compute regional response in terms of first eigenvariate
%--------------------------------------------------------------------------

% deal with nan
y(:,find(isnan(sum(y,1)))) = [];


[m,n]   = size(y);
if m > n
    [v,s,v] = svd(y'*y);
    s       = diag(s);
    v       = v(:,1);
    u       = y*v/sqrt(s(1));
else
    [u,s,u] = svd(y*y');
    s       = diag(s);
    u       = u(:,1);
    v       = y'*u/sqrt(s(1));
end
d       = sign(sum(v));
u       = u*d;
v       = v*d;
Y       = u*sqrt(s(1)/n);
 
%-Set in structure
%--------------------------------------------------------------------------
xY.y    = y;
xY.u    = Y;
xY.v    = v;
xY.s    = s;

xY.peakcenter = option.peakcenter;
    