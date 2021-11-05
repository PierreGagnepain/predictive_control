function xY = native_dcmvoi(P,roicoordinate,cname,nvoxel,analysiname)

% Function creating VOI for DCM :
%   - map MNI "roicoordinate" to subject's native space
%   - get voxel peak using T-value from the contrast defined in "cname"
%   - grow a new ROI within the subject's native space ROI mask, with
%   contiguous voxels with highest T-value
%   - extract 1st eigenvariate (and create VOI mask) -> output is stored in xY using "analysiname" 


spmpath = fullfile(P.glmpath,'SPM.mat');load(spmpath)

% Effects of interest contrast (to adjust VOI time-course)
%===========================================================================

% Create effect of interest contrast if does not exist
conname   = 'Effects of interest';
cf      = find(cellfun(@isempty,strfind({SPM.xCon.name},conname)) == 0);
ncond   = length(SPM.xCon);
if isempty(cf)
    nH      = length(SPM.xCon(1).c);
    noncon  = find(sum(horzcat(SPM.xCon(1:ncond).c),2) == 0);
    cn      = length(SPM.xCon);
    cn      = cn+1;
    c       = eye(nH);
    c(:,noncon) = 0;
    c(noncon,:) = [];
    if size(c,1)>1,  c=detrend(c,0); end
    
    SPM.xCon(cn)   = spm_FcUtil('Set',conname,'F','c',c',SPM.xX.xKXs);
    spm_contrasts(SPM,cn);
    cf = cn;
end
option.adjust = cf;
    
    
% loop through masks
%===========================================================================

clear xY

for m = 1:length(roicoordinate.coord)
    
    % contrast ?
    ctnt      = find(strcmp({SPM.xCon.name},cname{m}));
    
    % vol info
    vnatif      = spm_vol(fullfile(P.glmpath,sprintf('spmT_000%d.nii',ctnt)));
    
    % put mni mask in native space
    maskcoord   = roicoordinate.coord{m};
    native_c    = mni2native_v2(P.invfile,maskcoord');
    vox         = round(mm2vx_mat (native_c, vnatif.mat))';
    id          = unique(sub2ind(vnatif.dim,vox(:,1),vox(:,2),vox(:,3)));
    
    % get max peak coordinate
    [T, XYZmm]  = spm_read_vols(vnatif);
    Tmask       = T(id);
    [Ts idT]    = sort(Tmask,'descend');
    idx         = id(idT(1));% take max
    peakcenter  = XYZmm(:,idx);
    voxcenter   = round(mm2vx_mat (peakcenter, vnatif.mat))';
    
    % create adaptative voi mask
    mask                    = zeros(size(T));
    mask(id)                = 1;
    newRoi                  = GrowRoiFromPeak(voxcenter,T,nvoxel,mask);
    roi_index               = sub2ind(size(T),newRoi(:,1),newRoi(:,2),newRoi(:,3)); %single indices to MAP specifying voxels in the roi
    
    if m == 1
        newRoi_mask             = zeros(size(T));
        newRoi_mask(roi_index)  = 1;
    else
        newRoi_mask(roi_index)  = m;
    end
    
    % extract voi
    option.name             = sprintf('VOI_%s',roicoordinate.name{m});
    option.vox              = newRoi;
    option.peakcenter       = peakcenter;
    
    xY (m)                  = create_voi(SPM,option,P);
end

% save VOI.mat
fn = fullfile(P.voi_dir,sprintf('%s_VOI_%s.mat',P.subcode,analysiname));
save(fn,'xY')


% create voi mask
fn          = fullfile(P.voi_dir,sprintf('%s_VOI_%s.nii',P.subcode,analysiname));
vmask       = vnatif;
vmask.fname = fn;

spm_write_vol(vmask,newRoi_mask);
