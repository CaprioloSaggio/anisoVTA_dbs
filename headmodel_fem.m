%% initialize workspace and command line
clear
ft_defaults
clc
image = 'C:\Users\Notebook\Desktop\Downloads\DBS\03_Data\DICOM_raw_Alba\fsl_analyses\anat_t1_tra.nii';
disp('########## Starting FEM headmodel generation ##########')

%% read image
disp('########## Reading image ##########')
mri = ft_read_mri(image);

%% reslice the volume
disp('########## Reslicing the volume ##########')
cfg     = [];
cfg.dim = mri.dim;
mri     = ft_volumereslice(cfg, mri);

%% show intermediate result
% cfg=[];
% ft_sourceplot(cfg,mri);
% wf = input('Are you ok with the result? Answer: [yes], no, or stop')    % TODO: implement a difference in behaviour between stop (return) and no (repeat last step)
% if not(wf(1) == 'y' || wf(1) == 'Y') disp('Process interrupted by user'); return
% end

%% segment volume

% probably this part is not useful, since I don't use the tissue
% information to fetch the conductivity of the voxel
disp('########## Segmenting the volume ##########')
cfg           = [];
cfg.output    = {'gray','white','csf','skull','scalp'};
% cfg.output = {'head'};
segmentedmri  = ft_volumesegment(cfg, mri);

% save segmentedmri segmentedmri

%% put together the segmentation masks in a single indexed volume
% disp('########## Building an indexed volume comprehensive of all the segmented regions ##########')
% seg_i = ft_datatype_segmentation(segmentedmri,'segmentationstyle','indexed');

%% build the mesh
disp('########## Preparing the mesh ##########')
cfg        = [];
% cfg.shift  = 0.3;
cfg.shift = 0;

% build the 
cfg.method = 'rg_grid';
mesh_grid = rg_ft_prepare_mesh(cfg,segmentedmri);

% build hexahedral mesh
% cfg.method = 'hexahedral';
% mesh_hex = ft_prepare_mesh(cfg,segmentedmri);

% build tetrahedral mesh
cfg.method = 'tetrahedral';
mesh_tet = ft_prepare_mesh(cfg,segmentedmri);

%% find correspondence with conductivity values

% build the KDtree
disp('########## Converting grid to KD-Tree data structure ##########')
mesh_grid.ctr = KDTreeSearcher(mesh_grid.pos);

% compute center of each tetrahedral element in the mesh defining a new
% array of 3D points
disp('########## Computing grid based on tetrahedral mesh (looking for central node for each element) ##########')
mesh_tet.ctr = reshape(mesh_tet.pos(mesh_tet.tet(1:size(mesh_tet.tet,1),1:4),:), ...
                       [size(mesh_tet.tet, 1), 4, 3]);                  
mesh_tet.ctr = squeeze(mean(mesh_tet.ctr, 2));

% run Nearest Neighbours
disp('########## Running Nearest Neighbour to find correspondences ##########')
cond_i = knnsearch(mesh_grid.ctr, mesh_tet.ctr); % contains the index of the 
                                                 % conductivities that have
                                                 % a match in the mesh
mesh_tet = rmfield(mesh_tet, ctr);  % clean structure
                                                 

%% build headmodel

% load conductivity tensor
load('C:\Users\Notebook\Desktop\Downloads\DBS\03_Data\DICOM_raw_Alba\fsl_analyses\r_cond.mat')

% rearrange conductivity tensor
disp('########## Rearranging conductivity tensor ##########')
npnt = size(r_cond);
cond_t = zeros([prod(npnt(1:3)), 6]);
for i=1:6
    cond3d = r_cond(:,:,:,i);   % this passage is added for more readibility
    cond_t(:, i) = reshape(cond3d, [numel(cond3d), 1]);
end

% select only the conductivities that match elements in the tetrahedral mesh
cond = cond_t(cond_i, :);

% prepare headmodel
disp('########## Preparing the headmodel ##########')
cfg        = [];
cfg.method ='simbio';
% cfg.conductivity = [0.33 0.14 1.79 0.01 0.43];   % order follows mesh.tissuelabel
cfg.conductivity = cond;
mesh_tet.tissue = [];
mesh_tet.tissuelabel = [];

mesh_tet.tet(:, [3, 4]) = mesh_tet.tet(:, [4, 3]);  % necessary not to get 
                                                    % an error from sb_calc_stiff 
                                                    % relative to orientation
tic
vol = ft_prepare_headmodel(cfg, mesh_tet);
toc

%% show result
% disp('########## Showing results ##########')
% ft_plot_mesh(mesh, 'surfaceonly', 'yes');
% Achtung!!! This section is extremely heavy and MATLAB crashes (with 9.338M elements in the mesh)