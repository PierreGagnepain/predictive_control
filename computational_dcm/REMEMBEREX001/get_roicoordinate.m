
% atlas info
atlasfile           = '/Volumes/Transcend/Brainnetome_atlas/BN_Atlas_246_1mm.nii';
[atlasvol,XYZmm]    = spm_read_vols(spm_vol(atlasfile));
roiidx      = {};
roiidx{1}   = 216;
roiidx{2}   = 218;
roiidx{3}   = 152;
roiidx{4}   = [20 22];
roiidx{5}   = [16 24];

roicoordinate = {};
roicoordinate.name = {'rHip','cHip','PC','aMFG','pMFG'};
% get roi coordinates
for r = 1:length(roiidx)
    
    coord = XYZmm(:,find(ismember(atlasvol,roiidx{r})));
    
    if r == 4
        coord(:,find(coord(2,:)<35)) = [];
    elseif r == 5
        coord(:,find(coord(2,:)>25)) = [];
    end
    roicoordinate.coord{r} = coord;
end

save('roicoordinate','roicoordinate')
