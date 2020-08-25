
b0 = 'b0.nii';  % here I need to extract the first volume from dti_data.nii
anat_t1 = 'anat_t1.nii';
tmp = spm_coreg(anat_t1, b0);  % it would require an additional input (options structure), since here it is omitted, default options are considered
tmat = spm_matrix(tmp(:)');

% dti = spm_vol('eddy_unwarped_images.nii');
dti = spm_vol(b0);
anat = spm_vol(anat_t1);

dti(:).mat = spm_get_space(b0, tmat\dti.mat);
T = anat.mat\tmat\dti.mat;  % transformation from b0 to dti (according to reasoning on the interpretation of the output given in spm_coreg)

spm_write_vol(dti, spm_read_vols(dti));

tic
T = affine3d(T');
toc

tic
n_images = min([size(image_nr, 4), 20]);

bvals = importdata('bvals');
bvals = bvals(:,1:n_images);
save('bvals_short', 'bvals')

bvecs = importdata('bvecs');
bvecs = bvecs(:,1:n_images);
save('bvecs_short', 'bvecs')

r_dti = [];
for i=1:n_images
    r_dti(:,:,:,i) = imwarp(image_nr(:,:,:,i), T);
end

niftiwrite(r_dti, 'r_dti.nii')
toc
