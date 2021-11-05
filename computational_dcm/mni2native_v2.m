function native_c = mni2native_v2(norm_file,mni)

V   = spm_vol(norm_file);
M1  = V.mat; % Mapping from voxels in Template to MNI space
vox = round(M1\[mni'; ones(1,size(mni,1))]);
N   = nifti(V.fname);

x = N.dat(:,:,:,1,1);
y = N.dat(:,:,:,1,2);
z = N.dat(:,:,:,1,3);

native_c = [];
for c = 1:size(mni,1)
    native_c(c,:) = [x(vox(1,c),vox(2,c),vox(3,c)),...
        y(vox(1,c),vox(2,c),vox(3,c)),...
        z(vox(1,c),vox(2,1),vox(3,c))];
end

