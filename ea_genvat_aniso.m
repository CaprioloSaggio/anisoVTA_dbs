function varargout = ea_genvat_aniso(varargin)
% this function is called at line 1057 of function ea_stimparams.m in the
% following fashion:
% [stimparams(1,side).VAT(el).VAT,volume]=feval(ea_genvat,elstruct(el).coords_mm,getappdata(handles.stimfig,'S'),side,options,stimname,options.prefs.machine.vatsettings.horn_ethresh,handles.stimfig);
% 
% TODO: implement a GUI to select the DTI to give in input to the function
% TODO: implement saving of the headmodel and the possibility to just load
%       it in case it is present
% TODO: check if it's better to run the function in MNI space or in patient
%       space and understand if the transformations I0m using concern
%       patient or MNI space

% TODO: input variables copy-pasted from ea_genvat_horn: check if they are 
% all useful after all the steps have been implemented
if nargin==5
    acoords=varargin{1};
    S=varargin{2};
    side=varargin{3};
    options=varargin{4};
    stimname=varargin{5};
    thresh=options.prefs.machine.vatsettings.horn_ethresh; %0.2;
    
elseif nargin==7
    acoords=varargin{1};
    S=varargin{2};
    side=varargin{3};
    options=varargin{4};
    stimname=varargin{5};
    thresh=varargin{6};
    lgfigure=varargin{7};

elseif nargin==1  % in this case I just return the name of the method
    if ischar(varargin{1}) 
        varargout{1} = 'Anisotropic';
        return
    end
end

dbg_vis = 0;
dbg = 1;
if dbg; tic; end


%% ________________________________________________________________________
%% BOUNDING CYLINDER
if max(S.amplitude{side})>4
    stretchfactor=0.75*(max(S.amplitude{side})/2.5);
else
    stretchfactor=0.5;
end

% define a bounding cylinder to which restrict the mesh
upper_bound = 15+(20*stretchfactor);
lower_bound = -20*stretchfactor;
el_o_orig = [0, 0, upper_bound];
el_o_etop = [0, 0, lower_bound];
cyltrisize = 0.01;              % the maximum triangle size of the bounding cyl
cyltetvol = 1;                  % the maximum tetrahedral volume in the mesh --> CRITICAL PARAMETER FOR THE RESOLUTION OF THE FIELD
cylradius = 40*stretchfactor;   % define the radius of the bounding cylinder
ndiv=50;                        % division of circle for the bounding cylinder

% define and transform the cylinder, directly obtaining the mesh
[cylinder.pos,cylinder.face,cylinder.tet]= meshacylinder(el_o_etop,el_o_orig,cylradius,cyltrisize,cyltetvol,ndiv);
cylinder.tet = cylinder.tet(:,1:4);  % eliminate last column, that is the one used for homogeneous representation
cylinder.ctr = tetctr(cylinder.pos, cylinder.tet);


%% ________________________________________________________________________
%% ELECTRODES MODEL
disp('########## Building electrode model ##########')

%% get the electrode model and the relative patient-specific transformation
% TODO: copy-pasted from ea_genvat_horn, check it!
resultfig=getappdata(lgfigure,'resultfig');

% Important to load in reco from a new since we need to decide whether to
% use native or template coordinates. Even when running in template space,
% the native coordinates are sometimes used (VTA is then calculated in native space and ported to template).
options.loadrecoforviz=1;
[coords_mm,trajectory,markers]=ea_load_reconstruction(options);
elstruct(1).coords_mm=coords_mm;
elstruct(1).coords_mm=ea_resolvecoords(markers,options);
elstruct(1).trajectory=trajectory;
elstruct(1).name=options.patientname;
elstruct(1).markers=markers;

elspec=getappdata(resultfig,'elspec');
coords=acoords{side};
setappdata(resultfig,'elstruct',elstruct);

% Add stretchfactor to elstruct simply for purpose of checking if headmodel
% changed. Larger stim amplitudes need larger bounding boxes so
% stretchfactor must be incorporated here.
if max(S.amplitude{side})>4
    elstruct.stretchfactor=0.75; %(max(S.amplitude{side})/10);
else
    elstruct.stretchfactor=0.5;
end

% compute transformation from general to patient specific electrode model
% (surface, containing info on insulation or contact)
[elfv,~,T,electrode]=ea_buildelfv(elspec,elstruct,side);

% load and trnasform volumetric mesh of the electrode (w/o information about 
% insulation or contact)
elmodel_path=[ea_getearoot,'templates',filesep,'electrode_models',filesep,elspec.matfname,'_vol.mat'];
elmodel = load(elmodel_path);

if dbg_vis
    % old_elmodel = elmodel;
    figure
    plot3(elmodel.node(:,1), elmodel.node(:,2), elmodel.node(:,3));
end


%% cut the electrode at the size of the bounding cylinder

outside_bounding_cylinder = [find(elmodel.node(:,3)>upper_bound); find(elmodel.node(:,3)<lower_bound)];  
face_outside_bounding_cylinder = sum(ismember(elmodel.face, outside_bounding_cylinder), 2);  % find faces with at least one node out of the boundaries
elmodel.face(face_outside_bounding_cylinder>0,:) = [];  % eliminate faces with nodes out of the boundaries
% it is unnecessarily complicated to remove also the unused nodes, so I leave
% them there, this way the correspondence with the elmodel.face indices is
% maintained

% TODO: give a try to inpolyhedron.m function, that tests if points are
% inside a mesh or outside (probably it is slower than my solution, but my
% solution doesn't look to return the correct result
% inside_bounding_cylinder = inpolyhedron(cylinder.tet, cylinder.pos, elmodel.node(:,1), elmodel.node(:,2), elmodel.node(:,3));  % it tells me I'm testing too many points and (39371) and it can make MATLAB crash

% ########## deprecated ##########
% inside_bounding_cylinder = find(elmodel.node((elmodel.node(:,3)<15+20*stretchfactor) & (-20*stretchfactor<elmodel.node(:,3))));  
% f_inside_bounding_cylinder = sum(~ismember(elmodel.face, inside_bounding_cylinder), 2);  % find faces with nodes out of the boundaries
% elmodel.face(f_inside_bounding_cylinder>0,:) = [];  % eliminate faces with nodes out of the boundaries
% % elmodel.node = elmodel.node(inside_bounding_cylinder,:);  % keep only nodes inside the boundary

if dbg_vis
    elmodel.ctr = tetctr(elmodel.node, elmodel.face+1);  % find the centroids of the electrode elements
    figure
    plot3(elmodel.ctr(:,1), elmodel.ctr(:,2), elmodel.ctr(:,3), 'bx');
    hold on
    plot3(cylinder.pos(:,1), cylinder.pos(:,2), cylinder.pos(:,3), 'ro')
end

%% transform the electrode model into patient (MNI?) space
elmodel.node = T*[elmodel.node, ones(size(elmodel.node,1),1)]';
elmodel.node = elmodel.node(1:3,:)';
elmodel.face = elmodel.face + 1;  % here the index starts from 0, while matlab starts from 1
elmodel.ctr = tetctr(elmodel.node, elmodel.face);  % find the centroids of the electrode elements

% ########## deprecated ##########
% % maybe the part concerning the insulation is useless, since I can put the
% % all electrode at insulating conductivity and then apply the mask only for
% % the contacts
% % put together the mesh containing all the insulation in the electrode
% insulation.faces = electrode.insulation(:).faces;
% insulation.vertices = electrode.insulation(:).vertices;

% put together the mesh containing all the contacts in the electrode
% contacts.faces = electrode.contacts(:).faces;
contacts_vertices = [];
for i=1:length(electrode.contacts)
    contacts_vertices = [contacts_vertices; electrode.contacts(i).vertices]; 
end

% check what mesh elements are inside the contact regions
% insulation.vol = ea_intriangulation(insulation.vertices, insulation.faces, elmodel.ctr);  % ########## deprecated ##########
% contacts.vol = ea_intriangulation(contacts.vertices, contacts.faces, elmodel.ctr);

% insulation.elements = elmodel.face(insulation.vol);  % ########## deprecated ##########
% contacts.elements = elmodel.face(contacts.vol);

% ########## deprecated ##########
% insulation_ctr = tetctr(insulation.vertices, insulation.faces);
% contacts_ctr = tetctr(contacts.vertices, contacts.faces);

% ########## deprecated ##########
% if dbg_vis
%     figure
%     plot3(insulation.vertices(:,1),insulation.vertices(:,2),insulation.vertices(:,3),'r*');
%     hold on
%     plot3(contacts.vertices(:,1), contacts.vertices(:,2), contacts.vertices(:,3), 'bo')
%     plot3(elmodel.ctr(:,1),elmodel.ctr(:,2),elmodel.ctr(:,3),'g.');
%     drawnow
% end

%% ________________________________________________________________________
%% ENCAPSULATION LAYER
% it is modelled assuming it as a 0.5mm thick layer around the whole length of
% the electrode lead
cylinder.r = vecnorm(cylinder.ctr(:,1:2)')';  % find radial distance of each point in the bounding cylinder from the center of the electrode
elrad = elspec.lead_diameter / 2;  % find radius of the electrode lead
encapsulation_thickness = 0.5;  % 0.5 mm according to Gunalan et al 2017
encaps_index = find(elrad<cylinder.r<(elrad+encapsulation_thickness));  % find all elements in the cylinder mesh that correspond to encapsulation tissue

if dbg_vis
    encaps_ctr = cylinder.ctr(cylinder.r<(elrad+encapsulation_thickness) & cylinder.r>elrad, :);
    figure
    plot3(cylinder.ctr(:,1), cylinder.ctr(:,2), cylinder.ctr(:,3), 'r.')
    hold on
    plot3(encaps_ctr(:,1), encaps_ctr(:,2), encaps_ctr(:,3), 'b.')
    legend('bounding cylinder nodes', 'fibrotic tissue nodes', 'location', 'northeast')
    title('encapsulation layer')
end

%% ________________________________________________________________________
%% BOUNDING CYLINDER (transform)

% transform using the same trnasformation of the electrodes (from voxel to patient space)
cylinder.pos = T*[cylinder.pos, ones(length(cylinder.pos), 1)]';
cylinder.pos = cylinder.pos(1:3,:)';
cylinder.ctr = T*[cylinder.ctr, ones(length(cylinder.ctr), 1)]';
cylinder.ctr = cylinder.ctr(1:3,:)';

if dbg_vis
    figure
    title('## bounding cylinder ##');
    plot3(cylinder.pos(:,1),cylinder.pos(:,2),cylinder.pos(:,3),'r*');
    hold on
    plot3(contacts_vertices(:,1), contacts_vertices(:,2), contacts_vertices(:,3), 'bo')
    plot3(elmodel.ctr(:,1),elmodel.ctr(:,2),elmodel.ctr(:,3),'g.');
end


%% ________________________________________________________________________
%% HEADMODEL
disp('########## Starting FEM headmodel generation ##########')


%% initialize (with cartoon data in debug mode)
% insert here the path of the coregistered anatomical image
% TODO: here I'm assuming that the normalized and coregistered image is the
% anat_*.nii, check if it's true
anat = dir([options.root, options.patientname, filesep, options.prefs.prenii_searchstring]);
anat = [options.root,options.patientname,filesep,anat(1).name];
if dbg
    anat = 'C:\Users\Notebook\Desktop\Downloads\DBS\03_Data\DICOM_raw_Alba\fsl_analyses\r_anat_t1_tra.nii';
end
% insert here the path of the diffusion tensor
dti = 'C:\Users\Notebook\Desktop\Downloads\DBS\03_Data\DICOM_raw_Alba\fsl_analyses\dti_tensor.nii';


%% read image
disp('########## Reading image ##########')
mri = MRIread(anat);  % using this FreeSurfer function I also directly extract the transformation from voxel space to patient space
dti = niftiread(dti);


%% build the grid of points corresponding to the voxels of the volume
% build grid from volume 
mesh_grid.pos = build_grid(mri);
% ########## Ask for confirmation: I'm using anatomical data to obtain
% ########## the grid correspondent to the conducitivity data, because
% ########## before I've performed coregistration. Is it correct? The size
% ########## of the 2 images is the same


%% find correspondence with conductivity values
% build the KDtree
disp('########## Converting grid to KD-Tree data structure ##########')
mesh_grid.ctr = KDTreeSearcher(mesh_grid.pos);

% run Nearest Neighbours
disp('########## Running Nearest Neighbour to find correspondences ##########')
cond_i = knnsearch(mesh_grid.ctr, cylinder.ctr); % contains the index of the 
                                                 % conductivities that have
                                                 % a match in the mesh                                                

%% build headmodel
% load conductivity tensor
% load('C:\Users\Notebook\Desktop\Downloads\DBS\03_Data\DICOM_raw_Alba\fsl_analyses\r_cond.mat')

% get the diffusion data and convert them to conduction data
r_cond = 0.736 * dti;  % TODO: add a reference for the value that I used for scaling

% rearrange conductivity tensor
disp('########## Defining the conductivity tensor ##########')
npnt = size(r_cond);
cond_t = zeros([prod(npnt(1:3)), 6]);
for i=1:6
    cond3d = r_cond(:,:,:,i);   % this passage is added for more readibility
    cond_t(:, i) = reshape(cond3d, [numel(cond3d), 1]);
end

% select only the conductivities that match elements in the tetrahedral mesh
cond = cond_t(cond_i, :);  % TODO: plot the mesh_tet.ctr to see if there is also background


%% find correspondence with electrodes and include them into the model
% I set the conductivity of all voxels being close to the centroid of an
% electrode element to the insulating value of 10^-16 S/m and then I apply the
% contact conductivity of 10^8 S/m to all the voxels near to the centroid
% of a contact element.
cond_contacts = 1e8;  % Salvador et al 2012 states this should be 2 S/m, but here it is transcranial stimulation
cond_insulation = 1e-16;

el_cond = knnsearch(cylinder.ctr, elmodel.ctr);
cond(unique(el_cond),:) = cond_insulation;  % insulation conductivity (isotropic)
% el_cond = knnsearch(cylinder.ctr, elmodel.ctr(contacts.vol,:));  % ########## deprecated ##########
el_cond = knnsearch(cylinder.ctr, contacts_vertices);
cond(unique(el_cond),:) = cond_contacts;  % contact conductivity (isotropic)

if dbg; toc; end


%% add encapsulation layer (fibrotic tissue forming around the electrodes)
cond_encapsulation = 0.07;  % found in "isotropic conductivities", but according to Gunalan et al 2017 it can be 0.05±0.2 S/m
cond(encaps_index,:) = cond_encapsulation; %#ok<FNDSB>


%% ________________________________________________________________________
%% COMPUTE CONDUCTION MODEL
disp('########## Computing the conduction model ##########')
cfg        = [];
cfg.method ='simbio';
cfg.conductivity = cond;


if dbg; tic; end
try
    vol.stiff = simbio_stiff_matrix(cond, cylinder);
catch
    cylinder.tet(:, [3, 4]) = cylinder.tet(:, [4, 3]);  % necessary not to get 
                                                        % an error from sb_calc_stiff 
                                                        % relative to orientation
    vol.stiff = simbio_stiff_matrix(cond, cylinder);
%     cylinder.tissue = [];
%     cylinder.tissuelabel = [];
%     vol = ft_prepare_headmodel(cfg, cylinder);
end
vol.type = 'simbio';

if dbg; toc; end

% TODO: save the conductivity model and then use it if the function manages
% to retrieve it


%% ________________________________________________________________________
%% COMPUTE POTENTIAL
disp('########## Computing the potential based on stimulation ##########')

%% find boundary points in the volume of interest
cylinder.boundary = cylinder.pos(boundary(cylinder.pos(:,1), cylinder.pos(:,2), cylinder.pos(:,3)),:);
if dbg_vis
    figure 
    plot3(cylinder.pos(:,1), cylinder.pos(:,2), cylinder.pos(:,3), 'b.')
    hold on
    plot3(cylinder.boundary(:,1), cylinder.boundary(:,2), cylinder.boundary(:,3), 'r.')
    title('boundary nodes')
    legend('nodes inside the volume of interest','boundary nodes', 'location', 'northeast')
end

%% organize information to feed into the ea_apply_dbs function for potential computation
switch side
    case 1
        sidec='R';
        cnts={'k0','k1','k2','k3','k4','k5','k6','k7'};
    case 2
        sidec='L';
        cnts={'k8','k9','k10','k11','k12','k13','k14','k15'};
end

if ~isfield(S, 'sources')
    S.sources=1:4;
end

% organize information about active contacts
for con = find(S.activecontacts{side})+8*(side-1)-1
    for source=S.sources
        if S.([sidec,'s',num2str(source)]).amp && S.([sidec,'s',num2str(source)]).(['k',num2str(con+8*(side-1))]).perc  % find active contact corresponding to active source
            % find the nodes of the active contact in scope
            active_contacts = electrode.contacts(con).vertices;
            
            % find elements in mesh corresponding to nodes of the active
            % contact in scope
            cylinder.active = unique(knnsearch(cylinder.ctr, active_contacts));
            
            % define the activeidx structure, that organizes the
            % information for stimulation in a way that fits ea_apply_dbs
            activeidx(source).con(con).ix = cylinder.ctr(cylinder.active,:);
            activeidx(source).con(con).pol=S.([sidec,'s',num2str(source)]).(['k',num2str(con+8*(side-1))]).pol;
            activeidx(source).con(con).perc=S.([sidec,'s',num2str(source)]).(['k',num2str(con+8*(side-1))]).perc;
            
            stimsource=S.([sidec,'s',num2str(source)]);
            constvol=stimsource.va==1; % constvol is 1 for constant voltage and 0 for constant current.
        end
    end
end
% I assume the active_contacts nodes should be dense enough not to leave 
% gaps in the mesh

% ########## I have to start from here today ##########

    for cnt=1:length(cnts)
        if constvol
            U(cnt)=(logical(stimsource.(cnts{cnt}).perc))*stimsource.amp; % do not split amplitude in constant voltage setting.
        else
            U(cnt)=(stimsource.(cnts{cnt}).perc/100)*stimsource.amp;
        end
        if stimsource.(cnts{cnt}).pol==1
            U(cnt)=U(cnt)*-1;
        end
    end

    Acnt=find(U); % active contact
    if ~isempty(Acnt)

        dpvx=coords(Acnt,:);

        volts=U(U~=0);

        % calculate voltage distribution based on dipole
        ea_dispt('Calculating voltage distribution...');
        SIfx=1000;

        if any(volts>0)
            unipolar=0;
            U=U/2;
            %volts=volts/2;
        else
            unipolar=1;
        end
        ix=[];
        voltix=[];

        cnt=1;
        for ac=Acnt
            ix=[ix;activeidx(source).con(ac).ix];  % in ea_genvat_horn activeidx comes from ea_mesh_electrode
            voltix=[voltix;repmat(U(ac),length(activeidx(source).con(ac).ix),1),...
                repmat(cnt,length(activeidx(source).con(ac).ix),1)];
            cnt=cnt+1;

        end

        if isempty(ix)
           rmdir([options.root,options.patientname,filesep,'current_headmodel'],'s'); % the least I can do at this point is to clean up the faulty headmodel.
           ea_error('Something went wrong. Active vertex index not found.');
        end

        if ~constvol
            voltix(:,1)=voltix(:,1)/1000; % from mA to A
            %voltix=voltix;
        end

        potential = ea_apply_dbs(vol,ix,voltix,unipolar,constvol,wmboundary); % output in V. 4 indexes insulating material.
        % save('results','mesh','vol','ix','voltix','unipolar','constvol','wmboundary','potential3v','potential3ma','gradient3v','gradient3ma');
    end
end  % function



function ctr = tetctr(pos, tet)
% INPUT
% mesh_tet: tetrahedral mesh containing at least the fields 'pos' (nodes 
%           coordinates, of size Nx3) and 'tet' (elements, of size Mx4)
% OUTPUT
% ctr: centroids of the tetrahedral elements in a mesh (size Mx3)

% ctr = reshape(pos(tet(1:lentet,1:4),:), [lentet, 4, 3]);
ctr = reshape(pos(tet,:), [size(tet, 1), 4, 3]);
ctr = squeeze(mean(ctr, 2));
end  % subfunction



% the following 2 functions are used in ea_apply_dbs in order to obtain the
% potential distribution
function [stiff,rhs] = ea_dbs(stiff,rhs,dirinodes,dirival)

diagonal = diag(stiff);
stiff = stiff + stiff';
rhs = rhs - stiff*dirival;
stiff(dirinodes,:) = 0.0;
stiff(:,dirinodes) = 0.0;
diagonal = -diagonal;
%diagonal(:) = 0;
diagonal(dirinodes) = 1.0;
stiff = stiff + spdiags(diagonal(:),0,length(diagonal),length(diagonal));
rhs(dirinodes) = dirival(dirinodes);
end  % subfunction


function x = ea_sb_solve(sysmat,vecb)

% SB_SOLVE
%
% $Id: sb_solve.m 8776 2013-11-14 09:04:48Z roboos $
try
    L = ichol(sysmat);
catch
    alpha = max(sum(abs(sysmat),2)./diag(sysmat))-2;
    L = ichol(sysmat, struct('type','ict','droptol',1e-3,'diagcomp',alpha));
end

%scalen
[~,x]=evalc('pcg(sysmat,vecb,10e-10,5000,L,L'',vecb)');
end  % subfunction

