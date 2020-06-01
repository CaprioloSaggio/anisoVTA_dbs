%% initialize workspace and command line
clear
ft_defaults
clc
image = 'C:\Users\Notebook\Desktop\Downloads\DBS\03_Data\DICOM_raw_Alba\fsl_analyses\anat_t1_tra.nii'
disp('########## Starting FEM headmodel generation ##########')

%% read image
disp('########## Reading image ##########')
mri = ft_read_mri(image);

%% reslice the volume
disp('########## Reslicing the volume ##########')
cfg     = [];
cfg.dim = mri.dim;
mri     = ft_volumereslice(cfg,mri);

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
cfg.shift  = 0.3;

% build hexahedral mesh
% cfg.method = 'rg_tetrahedral';
cfg.method = 'hexahedral';
% mesh = rg_ft_prepare_mesh(cfg,mri);
mesh_hex = ft_prepare_mesh(cfg,segmentedmri);

%build tetrahedral mesh
cfg.method = 'tetrahedral';
mesh_tet = ft_prepare_mesh(cfg,segmentedmri);

%% build headmodel
%load conductivity tensor
load('C:\Users\Notebook\Desktop\Downloads\DBS\03_Data\DICOM_raw_Alba\fsl_analyses\r_cond.mat')

disp('########## Rearranging conductivity tensor ##########')
cond_t = [];
for i=1:6
    cond3d = r_cond(:,:,:,i);
    cond_t = [ cond_t, reshape(cond3d, [numel(cond3d), 1]) ];
end

disp('########## Preparing the headmodel ##########')
cfg        = [];
cfg.method ='simbio';
% cfg.conductivity = [0.33 0.14 1.79 0.01 0.43];   % order follows mesh.tissuelabel
cfg.conductivity = cond_t;
vol        = ft_prepare_headmodel(cfg, mesh);

%% show result
% disp('########## Showing results ##########')
% ft_plot_mesh(mesh, 'surfaceonly', 'yes');
% Achtung!!! This section is extremely heavy and MATLAB crashes (with 9.338M elements in the mesh)