clear
clc

%% find correspondence between hexa- and tetrahedrals
% before this run headmodel_fem.m until section 'build the mesh' included
% with the settings to obtain both the hexahedral mesh and the tetrahedral
% one

% compute center of each hexahedral element in the mesh defining a new
% array of 3D points
disp('########## Computing grid based on hexahedral mesh ##########')
tic
mesh_hex.ctr = [];
for el = 1:size(mesh_hex.hex,1)  % very demanding computation
    m = [];                      % TODO: find a faster way!!
    for i = 1:8
        m = [m; mesh_hex.pos(mesh_hex.hex(el,i),:)];
    end
    m = mean(m);
    mesh_hex.ctr = [mesh_hex.pnt; m(1), m(2), m(3)];
end
toc

% build the KDtree
disp('########## Convert grid to KD-Tree data structure ##########')
tic
mesh_hex.ctr = KDTreeSearcher(mesh_hex.pnt)
toc
% mesh_hex.ctr = KDTreeSearcher(mesh_hex.ctr)

% compute center of each tetrahedral element in the mesh defining a new
% array of 3D points
disp('########## Computing grid based on tetrahedral mesh ##########')
tic
mesh_tet.ctr = [];
for el = 1:size(mesh_tet.tet,1)  % very demanding computation: takes some minutes
    m = [];
    for i = 1:4  % probably this cycle can be speeded up (e.g. by avoiding it)
        m = [m; mesh_tet.pos(mesh_tet.tet(el,i),:)];
    end
    m = mean(m);
    mesh_tet.ctr = [mesh_tet.pnt; m(1), m(2), m(3)];
end
toc

% run Nearest Neighbours
disp('########## Running Nearest Neighbour to find correspondences ##########')
tic
nn = knnsearch(mesh_hex.ctr, mesh_tet.ctr);
toc

%% resample electrode grid  % to be continued

% read data
[f, v] = stlread('C:\Users\Notebook\Desktop\Downloads\DBS\04_Source\00_Master\templates\electrode_models\bsc_vercise_directed.stl');

% define ranges
step = 0.1;
rx = min(v(:,1)):step:max(v(:,1));
ry = min(v(:,2)):step:max(v(:,2));
rz = min(v(:,3)):step:max(v(:,3));

% scan z dimension
levels = unique(v(:,3));
for lev = levels
    continue
end


%% represent electrode (lower part, so that the curve is visible)
fv = stlread('C:\Users\Notebook\Desktop\Downloads\DBS\04_Source\00_Master\templates\electrode_models\bsc_vercise_directed.stl');
points = fv.Points(find(fv.Points(:,3)<10),:);
plot3(points(:,1), points(:,2), points(:,3), 'o')


%% create hex mesh
load("C:\Users\Notebook\Desktop\Downloads\DBS\04_Source\00_Master\templates\electrode_models\boston_vercise_directed.mat")
[bb_x, bb_y, bb_z] = ind2sub(size(electrode.contacts(1).vertices), find(electrode.contacts(1).vertices));
shift_coord = [min(bb_x) - 2, min(bb_y) - 2, min(bb_z) - 2];
bb_x = [min(bb_x), max(bb_x)];
bb_y = [min(bb_y), max(bb_y)];
bb_z = [min(bb_z), max(bb_z)];
x_dim = size(bb_x(1)-1:bb_x(2)+1, 2);
y_dim = size(bb_y(1)-1:bb_y(2)+1, 2);
z_dim = size(bb_z(1)-1:bb_z(2)+1, 2);
% this is the version used in prepare_mesh_hexahedral --> build_mesh_hexahedral

% x_dim = max(electrode.contacts(1).vertices(:,1)) - min(electrode.contacts(1).vertices(:,1));
% y_dim = max(electrode.contacts(1).vertices(:,2)) - min(electrode.contacts(1).vertices(:,2));
% z_dim = max(electrode.contacts(1).vertices(:,3)) - min(electrode.contacts(1).vertices(:,3));
mesh.hex = create_elements(x_dim, y_dim, z_dim);
mesh.pos = create_nodes(x_dim, y_dim, z_dim);
ft_plot_mesh(mesh, 'surfaceonly', 'yes');

%% coregister before fsl pre-processing (anat to b0)  % to be continued

% load images
a_anat = 'anat_t1_tra.nii';
a_b0 = 'b0_AP.nii';  % TODO: try to start from the DICOM maybe

% load anatomical image
anat = spm_vol(a_anat)

% load diffusion image
b0 = spm_vol(a_b0);

% extract dimensions of the anatomical image


%% coregistration using SPM

% choose images
ref = 'C:\Users\Notebook\Desktop\Downloads\DBS\03_Data\DICOM_raw_Alba\fsl_analyses\anat_t1_tra.nii';
anat = niftiread(ref);
source = 'C:\Users\Notebook\Desktop\Downloads\DBS\03_Data\DICOM_raw_Alba\fsl_analyses\b0.nii';
dti = 'C:\Users\Notebook\Desktop\Downloads\DBS\03_Data\DICOM_raw_Alba\fsl_analyses\dti_tensor.nii';

% define configuration strucure
flags = [];
flags.cost_fun = 'nmi';
flags.sep = [4 2];
flags.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
flags.fwhm = [7 7];
flags.interp = 4;
flags.wrap = [0 0 0];
flags.mask = 1;
flags.prefix = 'reg_';

% get the diffusion data and convert them to conduction data
diff = niftiread(dti);
cond = 0.736 * diff;

% ____________________________________________________________________________________
%{
% run coregistration
x = spm_coreg(ref, source, flags);  % most of the computational burden lays here
T = spm_matrix(x(:)');

% get transformations for the image space of the 2 images
ref_T = spm_vol(ref).mat;        % a 4x4 affine transformation matrix mapping 
source_T = spm_vol(source).mat;  % from voxel coordinates to real world coordinates.

% finally obtain the transformation to apply
T = ref_T \ T * source_T;
T = affine3d(T');

% applt transformation to the 4D volume
for i=1:6
    cond3d = cond(:, :, :, i);
    if (i==1) niftiwrite(cond3d, 'r_cond1.nii'); end
    r_cond(:, :, :, i) = imwarp(cond3d, T);
end
clear cond3d i

% save as NIfTI
filename = 'C:\Users\Notebook\Desktop\Downloads\DBS\03_Data\DICOM_raw_Alba\fsl_analyses\r_cond.nii';
niftiwrite(r_cond, filename)
%}
% ____________________________________________________________________________________


% inspired to spm_run_coreg
x = spm_coreg(ref, source, flags);    
M  = spm_matrix(x);
% MM = zeros(4,4);
MM = spm_get_space(source);
spm_get_space(source, M\MM);
for i=1:6
    cond3d = cond(:, :, :, i);
    filename = sprintf('cond_%d.nii', i);
    niftiwrite(cond3d, filename);
    spm_get_space(filename, M\MM);
    P = {ref; filename};
    spm_reslice(P, flags);  % TODO: it saves a lot of files (2x6=12), I 
                            % should extract the useful code from spm_reslice 
                            % and use just that part 
end

for i=1:6
    filename = sprintf('reg_cond_%d.nii', i);
    cond3d = niftiread(filename);
    r_cond(:, :, :, i) = cond3d;
end
clear cond3d i

%% coregistration using Image Analysis Toolbox

% load images
b0 = niftiread('b0.nii');
b0 = fliplr(b0);
anat = niftiread('anat_t1_tra.nii');

% set parameters
imageSizeFixed = size(anat);
imageSizeMoving = size(b0);
    % --- this information must be extracted from DICOM headers 
pixelExtentInWorldZfixed = 1;
pixelExtentInWorldZmoving = 4.6;
pixelExtentInWorldX = 1;
pixelExtentInWorldXmoving = 2;
pixelExtentInWorldY = 1;
pixelExtentInWorldYmoving = 2;
    % ---
Rfixed = imref3d(imageSizeFixed,pixelExtentInWorldX,pixelExtentInWorldY,pixelExtentInWorldZfixed);
Rmoving = imref3d(imageSizeMoving,pixelExtentInWorldXmoving,pixelExtentInWorldYmoving,pixelExtentInWorldZmoving);
transformType = 'affine';
[optimizer, metric] = imregconfig('multimodal');
%optimizer.MaximumStepLength = 0.05;

% perform registration
fprintf('/n ########## Performing Coregistration ########## /n')
[reg, T] = imregister(b0,Rmoving,anat,Rfixed,transformType,optimizer,metric);

% datavis
i = 170;
figure
imshowpair(anat(:,:,i), reg(:,:,i),'Scaling','joint')


%% extract DICOM information
clear
clc

path = "C:\Users\Notebook\Desktop\Downloads\DBS\03_Data\DICOM_raw_Alba\DICOM";
cd(path)
folder = "0000FFD6";
directory = fullfile(path, folder);
[b_val, b_vect] = extract_b(directory, 1);



%% test fanDTasia  % to be continued
% image = niftiread("C:\Users\Notebook\Desktop\Downloads\DBS\03_Data\DW-MRI Alba\DTI_b1250\b1250_AP.nii")
S=openFDT('C:\Users\Notebook\Desktop\Downloads\DBS\03_Data\fandtasia_demo\fandtasia_demo.fdt');
params=[0.517646 -0.000000 0.855595 29.235000
0.000000 -0.000000 1.000000 1240.777997
0.348230 -0.000000 0.937409 1234.911011
0.108361 0.326375 0.939008 1234.871010
-0.281838 0.201529 0.938058 1234.735004
-0.281855 -0.201251 0.938113 1234.591007
0.108372 -0.326106 0.939101 1234.629006
0.476321 0.341043 0.810437 1221.494964
-0.183095 0.553680 0.812351 1221.224972
-0.586871 0.000000 0.809680 1220.737976
-0.183125 -0.553475 0.812484 1220.830990
0.476368 -0.340785 0.810518 1221.251968
0.651703 -0.000000 0.758474 1215.924988
0.203676 0.615476 0.761384 1216.274990
-0.528686 0.378981 0.759516 1216.029968
-0.528745 -0.378730 0.759600 1215.760971
0.203716 -0.615286 0.761527 1215.823971
0.766785 0.331304 0.549799 1194.685013
0.558232 0.620463 0.550821 1195.053009
-0.082953 0.829150 0.552836 1194.882018
-0.426261 0.717082 0.551449 1195.067994
-0.816442 0.177303 0.549533 1194.216034
-0.816486 -0.177008 0.549563 1194.086033
-0.426355 -0.716931 0.551571 1194.548035
-0.082976 -0.829053 0.552979 1194.282042
0.558339 -0.620273 0.550926 1194.601990
0.766862 -0.331032 0.549855 1194.443009
0.875633 -0.000000 0.482977 1188.291016
0.274327 0.829877 0.485848 1188.505997
-0.710602 0.510656 0.484021 1188.111999
-0.710713 -0.510427 0.484098 1187.740997
0.274396 -0.829783 0.485970 1187.925980
0.789811 0.567942 0.231606 1163.350007
-0.305163 0.923393 0.232853 1163.144005
-0.972927 -0.000000 0.231112 1162.549042
-0.305250 -0.923347 0.232920 1162.506981
0.789946 -0.567737 0.231647 1162.956025
0.957350 0.208323 0.200208 1160.565010
0.500841 0.841825 0.201219 1160.860012
0.097618 0.974486 0.202108 1160.720028
-0.654304 0.729044 0.200950 1160.135983
-0.899185 0.388953 0.200457 1159.795025
-0.899294 -0.388688 0.200482 1159.516033
-0.654454 -0.728897 0.200997 1159.616024
0.097649 -0.974470 0.202173 1160.027035
0.500974 -0.841733 0.201273 1160.260036
0.957409 -0.208037 0.200220 1160.421013];
GradientOrientations=params(:,1:3);
b_value=params(:,4);
