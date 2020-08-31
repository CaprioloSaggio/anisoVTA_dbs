% Here's how to check the conductivity values of the STN, expected to be
% around 0.33, as it is grey matter

% find STN data
path = 'C:\Users\Notebook\Desktop\Downloads\DBS\04_Source\00_Master\templates\templates\space\MNI_ICBM_2009b_NLIN_ASYM\atlases\DISTAL Minimal (Ewert 2017)';
load([path, filesep, 'atlas_index.mat']);
stnl = atlases.fv{1, 1}.vertices;
stnr = atlases.fv{1, 2}.vertices;
save('stn.mat', 'stnl', 'stnr');

% load stn points
load stn

% find conductivity in the STN
test = unique(knnsearch(vol.ctr, stnl));
test_c = cond(test, :);
test_cm = mean(test_c);

% compute the isotropic conductivity in the STN, in order to compare it
% with the reference values for  GM in literature (0.2 or 0.3)
tensor = zeros([numel(r_cond)/6, 3, 3]);
for i=1:6
    tensor(:,1,1) = reshape(r_cond(:,:,:,1), [numel(r_cond)/6,1]);
    tensor(:,1,2) = reshape(r_cond(:,:,:,2), [numel(r_cond)/6,1]);
    tensor(:,2,1) = reshape(r_cond(:,:,:,2), [numel(r_cond)/6,1]);
    tensor(:,1,3) = reshape(r_cond(:,:,:,3), [numel(r_cond)/6,1]);
    tensor(:,3,1) = reshape(r_cond(:,:,:,3), [numel(r_cond)/6,1]);
    tensor(:,2,2) = reshape(r_cond(:,:,:,4), [numel(r_cond)/6,1]);
    tensor(:,2,3) = reshape(r_cond(:,:,:,5), [numel(r_cond)/6,1]);
    tensor(:,3,2) = reshape(r_cond(:,:,:,5), [numel(r_cond)/6,1]);
    tensor(:,3,3) = reshape(r_cond(:,:,:,6), [numel(r_cond)/6,1]);
end
tensor_vol = tensor(cond_i,:,:);
test = unique(knnsearch(vol.ctr, stnl));
tensor_stn = tensor(test,:,:);

eigval = []; eigvec = []; eigenvalue = []; iso_cond = []; 
for i=1:size(tensor_stn,1)
    [eigvec(i,:,:), eigval(i,:,:)] = eig(squeeze(tensor_stn(i,:,:)));
    eigenvalue(i,:) = diag(squeeze(eigval(i,:,:)));
    iso_cond(i) = (prod(eigenvalue(i,:)))^(1/3);
end

min(iso_cond)
max(iso_cond)
mean(iso_cond)
median(iso_cond)    