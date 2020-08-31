function nodes = build_grid(mri)
% INPUT:
% mri: 3D image data obtained by using MRIread
% 
% OUTPUT:
% nodes: nodes of a equispaced grid along the whole volume of the image

x_dim = mri.dim(1);
y_dim = mri.dim(2);
z_dim = mri.dim(3);

nodes = zeros(((x_dim)*(y_dim)*(z_dim)), 3);
% offset vector for node coordinates
b = 0:((x_dim)*(y_dim)*(z_dim)-1);

% assign coordinates within the box
nodes(:, 1) = mod(b, (x_dim));
nodes(:, 2) = mod(fix(b/(x_dim)), (y_dim));
nodes(:, 3) = fix(b/((x_dim)*(y_dim)));
nodes = nodes + .5;

% mri.mat is the field containing the transformation from voxel to patient space
T = mri.mat; 

% converting position of meshpoints to the head coordinate system
nodes = ft_warp_apply(T, nodes, 'homogeneous');  

end