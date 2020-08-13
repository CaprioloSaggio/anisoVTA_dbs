function varargout = ea_genvat_aniso(varargin)
% this function is called at line 1059 of function ea_stimparams.m in the
% following fashion:
% [stimparams(1,side).VAT(el).VAT,volume]=feval(ea_genvat,elstruct(el).coords_mm,getappdata(handles.stimfig,'S'),side,options,stimname,options.prefs.machine.vatsettings.horn_ethresh,handles.stimfig);
% 
% TODO: implement a GUI to select the DTI to give in input to the function
% TODO: implement saving of the headmodel and the possibility to just load
%       it in case it is present
% TODO: check if it's better to run the function in MNI space or in patient
%       space and understand if the transformations I'm using concern
%       patient or MNI space
%

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


%% debug settings
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
[vol.pos,vol.face,vol.tet]= meshacylinder(el_o_etop,el_o_orig,cylradius,cyltrisize,cyltetvol,ndiv);
vol.tet = vol.tet(:,1:4);  % eliminate last column, that is the one used for homogeneous representation
vol.ctr = tetctr(vol.pos, vol.tet);


%% ________________________________________________________________________
%% ELECTRODES MODEL
ea_dispt('Building electrode model...')


%% get the electrode model and the relative patient-specific transformation
% Important to load in reco from a new since we need to decide whether to
% use native or template coordinates. Even when running in template space,
% the native coordinates are sometimes used (VTA is then calculated in 
% native space and ported to template).

resultfig=getappdata(lgfigure,'resultfig');
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
[~,~,T,electrode]=ea_buildelfv(elspec,elstruct,side);  % T should be the transformation between model space and patient space

% load and transform volumetric mesh of the electrode (w/o information about 
% insulation or contact)
elmodel_path=[ea_getearoot,'templates',filesep,'electrode_models',filesep,elspec.matfname,'_vol.mat'];
elmodel = load(elmodel_path);

if dbg_vis
    % old_elmodel = elmodel;
    figure
    plot3(elmodel.node(:,1), elmodel.node(:,2), elmodel.node(:,3));
end


%% cut the electrode at the size of the bounding cylinder
outside_bounding_vol = [find(elmodel.node(:,3)>upper_bound); find(elmodel.node(:,3)<lower_bound)];  
face_outside_bounding_vol = sum(ismember(elmodel.face, outside_bounding_vol), 2);  % find faces with at least one node out of the boundaries
elmodel.face(face_outside_bounding_vol>0,:) = [];  % eliminate faces with nodes out of the boundaries
% it is unnecessarily complicated to remove also the unused nodes, so I leave
% them there, this way the correspondence with the elmodel.face indices is
% maintained

% TODO: give a try to inpolyhedron.m function, that tests if points are
% inside a mesh or outside (probably it is slower than my solution, but my
% solution doesn't look to return the correct result
% inside_bounding_vol = inpolyhedron(vol.tet, vol.pos, elmodel.node(:,1), elmodel.node(:,2), elmodel.node(:,3));  % it tells me I'm testing too many points and (39371) and it can make MATLAB crash

% ########## deprecated ##########
% inside_bounding_vol = find(elmodel.node((elmodel.node(:,3)<15+20*stretchfactor) & (-20*stretchfactor<elmodel.node(:,3))));  
% f_inside_bounding_vol = sum(~ismember(elmodel.face, inside_bounding_vol), 2);  % find faces with nodes out of the boundaries
% elmodel.face(f_inside_bounding_vol>0,:) = [];  % eliminate faces with nodes out of the boundaries
% % elmodel.node = elmodel.node(inside_bounding_vol,:);  % keep only nodes inside the boundary

if dbg_vis
    elmodel.ctr = tetctr(elmodel.node, elmodel.face+1);  % find the centroids of the electrode elements
    figure
    plot3(elmodel.ctr(:,1), elmodel.ctr(:,2), elmodel.ctr(:,3), 'bx');
    hold on
    plot3(vol.pos(:,1), vol.pos(:,2), vol.pos(:,3), 'ro')
end


%% transform the electrode model into patient space
elmodel.node = T*[elmodel.node, ones(size(elmodel.node,1),1)]';
elmodel.node = elmodel.node(1:3,:)';
elmodel.face = elmodel.face + 1;  % here the index starts from 0, while matlab starts from 1
elmodel.ctr = tetctr(elmodel.node, elmodel.face);  % find the centroids of the electrode elements

% I don't need to take into account the insulation, since I can put the
% all electrode at insulating conductivity and then apply the mask only for
% the contacts

% ########## deprecated ##########
% % put together the mesh containing all the insulation in the electrode
% insulation.faces = electrode.insulation(:).faces;
% insulation.vertices = electrode.insulation(:).vertices;
% insulation_ctr = tetctr(insulation.vertices, insulation.faces);
% contacts_ctr = tetctr(contacts.vertices, contacts.faces);


% put together the mesh containing all the contacts in the electrode
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
vol.r = vecnorm(vol.ctr(:,1:2)')';  % find radial distance of each point in the bounding cylinder from the center of the electrode
elrad = elspec.lead_diameter / 2;  % find radius of the electrode lead
encapsulation_thickness = 0.5;  % 0.5 mm according to Gunalan et al 2017
encaps_index = find(vol.r>elrad & vol.r<(elrad+encapsulation_thickness));  % find all elements in the cylinder mesh that correspond to encapsulation tissue

if dbg_vis
    encaps_ctr = vol.ctr(vol.r<(elrad+encapsulation_thickness) & vol.r>elrad, :);
    figure
    plot3(vol.ctr(:,1), vol.ctr(:,2), vol.ctr(:,3), 'r.')
    hold on
    plot3(encaps_ctr(:,1), encaps_ctr(:,2), encaps_ctr(:,3), 'b.')
    legend('bounding cylinder nodes', 'fibrotic tissue nodes', 'location', 'northeast')
    title('encapsulation layer')
end


%% ________________________________________________________________________
%% BOUNDING CYLINDER (transform)

% transform using the same trnnsformation of the electrodes (from model to patient space)
vol.pos = T*[vol.pos, ones(length(vol.pos), 1)]';
vol.pos = vol.pos(1:3,:)';
vol.ctr = T*[vol.ctr, ones(length(vol.ctr), 1)]';
vol.ctr = vol.ctr(1:3,:)';

if dbg_vis
    figure
    title('## bounding cylinder ##');
    plot3(vol.pos(:,1),vol.pos(:,2),vol.pos(:,3),'r*');
    hold on
    plot3(contacts_vertices(:,1), contacts_vertices(:,2), contacts_vertices(:,3), 'bo')
    plot3(elmodel.ctr(:,1),elmodel.ctr(:,2),elmodel.ctr(:,3),'g.');
end


%% ________________________________________________________________________
%% HEADMODEL
ea_dispt('Starting FEM headmodel generation...')


%% initialize (with cartoon data in debug mode)
% insert here the path of the coregistered anatomical image
% TODO: here I'm assuming that the normalized and coregistered image is the
% anat_*.nii, check if it's true
anat = dir([options.root, options.patientname, filesep, options.prefs.prenii_searchstring]);
anat = [options.root,options.patientname,filesep,anat(1).name];  % anat is the path to the file called anat_t1.nii in the patient folder
dti = [options.root,options.patientname,filesep,'dti.nii'];  % dti is the path to the file called dti.nii in the patient folder
if dbg
    anat = 'C:\Users\Notebook\Desktop\Downloads\DBS\03_Data\DICOM_raw_Alba\fsl_analyses\r_anat_t1_tra.nii';
end
% insert here the path of the diffusion tensor
dti = 'C:\Users\Notebook\Desktop\Downloads\DBS\03_Data\DICOM_raw_Alba\fsl_analyses\dti_tensor.nii';


%% read image
ea_dispt('Reading image...')
mri = MRIread(anat);  % using this FreeSurfer function I also directly extract the transformation from voxel space to patient space
dti = niftiread(dti);


%% build the grid of points corresponding to the voxels of the volume
% build grid from volume 
mesh_grid.pos = build_grid(mri);
% ########## Ask for confirmation: I'm using anatomical data to obtain
% ########## the grid correspondent to the conductivity data, because
% ########## before I've performed coregistration. Is it correct? The size
% ########## of the 2 images is the same


%% find correspondence with conductivity values
ea_dispt('Defining the conductivity tensor...')

% build the KDtree
mesh_grid.ctr = KDTreeSearcher(mesh_grid.pos);

% run Nearest Neighbours
% disp('########## Running Nearest Neighbour to find correspondences ##########')
cond_i = knnsearch(mesh_grid.ctr, vol.ctr); % contains the index of the 
                                            % conductivities that have
                                            % a match in the mesh                                                

%% build headmodel
% load conductivity tensor
% load('C:\Users\Notebook\Desktop\Downloads\DBS\03_Data\DICOM_raw_Alba\fsl_analyses\r_cond.mat')

% get the diffusion data and convert them to conduction data
r_cond = abs(dti) * 0.736;  % McIntyre et al 2004
r_cond = r_cond * 1e3;  % from S/mm to S/m

% rearrange conductivity tensor
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
cond_contacts = 1e8;  % Salvador et al 2012 states this should be 2 S/m, but in that case it is transcranial stimulation
cond_insulation = 1e-16;

el_cond = unique(knnsearch(vol.ctr, elmodel.ctr));  % index of all the elements in cylinder corresponding to 
cond(el_cond,:) = cond_insulation;  % insulation conductivity (isotropic)

con_cond = unique(knnsearch(vol.ctr, contacts_vertices));
cond(con_cond,:) = cond_contacts;  % contact conductivity (isotropic)


%% add encapsulation layer (fibrotic tissue forming around the electrodes)
cond_encapsulation = 0.07;  % found in "isotropic conductivities", but according to Gunalan et al 2017 it can be 0.05±0.2 S/m
cond(encaps_index,:) = cond_encapsulation; %#ok<FNDSB>


%% ________________________________________________________________________
%% COMPUTE CONDUCTION MODEL
ea_dispt('Computing the conduction model...')
cfg        = [];
cfg.method ='simbio';
cfg.conductivity = cond;

try
    vol.stiff = simbio_stiff_matrix(cond, vol);
catch
    vol.tet(:, [3, 4]) = vol.tet(:, [4, 3]);  % necessary not to get 
                                                        % an error from sb_calc_stiff 
                                                        % relative to orientation
    vol.stiff = simbio_stiff_matrix(cond, vol);
end
vol.type = 'simbio';
vol.stiff = abs(vol.stiff);  % TODO: check it, I had to add this line in order to make ea_apply_dbs execute ichol


%% ________________________________________________________________________
%% COMPUTE POTENTIAL AND ITS GRADIENT
ea_dispt('Computing the potential based on stimulation...')


%% find boundary points in the volume of interest
% TODO: include this section in the function to build the cylinder mesh (it should speed up a little)
vol.boundary = unique(boundary(vol.pos(:,1), vol.pos(:,2), vol.pos(:,3)));
if dbg_vis
    figure 
    plot3(vol.pos(:,1), vol.pos(:,2), vol.pos(:,3), 'b.')
    hold on
    plot3(vol.pos(vol.boundary, 1), vol.pos(vol.boundary, 2), vol.pos(vol.boundary, 3), 'r.')
    title('boundary nodes')
    legend('nodes inside the volume of interest','boundary nodes', 'location', 'northeast')
end


%% define sources and contacts
if ~isfield(S, 'sources')
    S.sources=1:4;
end

switch side
    case 1
        sidec='R';
        cnts={'k0','k1','k2','k3','k4','k5','k6','k7'};
    case 2
        sidec='L';
        cnts={'k8','k9','k10','k11','k12','k13','k14','k15'};
end

% preallocate memory to store the gradient
gradient = cell(length(S.sources), 1);


%%
dpvx = [];
for source = S.sources  % ##### PAY ATTENTION: this block allows for multipolar stimulation, although it was not required by the project #####
    active_contacts = [];
    ix = [];
    voltix = [];
    stimsource = S.([sidec,'s',num2str(source)]);
    constvol = stimsource.va==1; % constvol is 1 for constant voltage and 0 for constant current.
    ac = find(S.activecontacts{side});
    U = zeros(length(ac),1);
    for con = ac  
        if S.([sidec,'s',num2str(source)]).amp && S.([sidec,'s',num2str(source)]).(['k',num2str(con+8*(side-1)-1)]).perc  % find active contact corresponding to active source
            %% organize information to feed into the ea_apply_dbs function for potential computation
            
            % find the nodes of the active contact in scope
            active_contacts = electrode.contacts(con).vertices;

            % find elements in mesh corresponding to nodes of the active
            % contact in scope
            vol.active = unique(knnsearch(vol.pos, active_contacts));  % be careful to if they leave a gap in between. In case maybe it is detrimental
            
            % define the activeidx structure, that organizes the
            % information for stimulation in a way that fits ea_apply_dbs
            activeidx(source).con(con).ix = vol.active;
            activeidx(source).con(con).pol = S.([sidec,'s',num2str(source)]).(['k',num2str(con+8*(side-1)-1)]).pol;
            activeidx(source).con(con).perc = S.([sidec,'s',num2str(source)]).(['k',num2str(con+8*(side-1)-1)]).perc;
            
            if constvol
                U(con)=(logical(stimsource.(cnts{con}).perc))*stimsource.amp; % do not split amplitude in constant voltage setting.
            else
                U(con)=(stimsource.(cnts{con}).perc/100)*stimsource.amp;
            end
            
            dpvx = [dpvx; coords(boolean(U), :)];  % find coordinates where the contacts are active

            
            if stimsource.(cnts{con}).pol==1
                U(con)=U(con)*-1;
            end
            
            ix = [ix;activeidx(source).con(con).ix]; %#ok<*AGROW>
            voltix_new = [repmat(U(con), length(activeidx(source).con(con).ix), 1), ...
                      repmat(con, length(activeidx(source).con(con).ix), 1)];   
            if ~constvol
                voltix_new(:,1) = voltix_new(:,1) / 1000;  % from mA to A
            end
            voltix = [voltix; voltix_new];
            
            if dbg_vis
                figure
                plot3(elmodel.ctr(:,1), elmodel.ctr(:,2), elmodel.ctr(:,3), 'b.')
                hold on
                plot3(vol.pos(ix,1), vol.pos(ix,2), vol.pos(ix,3), 'r.', 'MarkerSize', 20)
                plot3(vol.pos(:,1), vol.pos(:,2), vol.pos(:,3), 'g.')
                legend('electrode', 'active', 'volume', 'location', 'northeast')
                title('active elements')
            end
        end
    end

    
    %% compute potential distribution and gradient
    if ~(isempty(active_contacts))
        if any(U>0) 
            unipolar=0;
            U=U/2;
        else
            unipolar=1;
        end
        
        potential = ea_apply_dbs(vol,ix,voltix,unipolar,constvol,vol.boundary); % output in V. 4 indexes insulating material.
        gradient{source} = ea_calc_gradient(vol,potential);


        %% get high EF values for active electrodes
        % this can be adjusted by assigning all tetrahedar belonging to the
        % active electrode a new value: gradient{source}(elec_tet_ix,:) = new_value;
        elec_tet_ix = sub2ind(size(vol.pos), vertcat(ix,ix,ix), vertcat(ones(length(ix),1), ones(length(ix),1).*2, ones(length(ix),1).*3));
        elec_tet_ix = find(sum(ismember(vol.tet,elec_tet_ix),2)==4);

        % gradient{source}(elec_tet_ix,:) = repmat(max(gradient{source}),[length(elec_tet_ix),1]); %assign maximum efield value
        tmp = sort(abs(gradient{source}),'descend');
        gradient{source}(elec_tet_ix,:) = repmat(mean(tmp(1:ceil(length(tmp(:,1))*0.001),:)),[length(elec_tet_ix),1]); % choose mean of highest 0.1% as new efield value
        clear tmp


    else
        % if the source is not active, the solution is trivial
        gradient{source} = zeros(size(vol.tet,1),3);
    end
end

% combine gradients from all sources
gradient=gradient{1}+gradient{2}+gradient{3}+gradient{4}; 

% convert back to mm
vol.pos=vol.pos*1000; 

vatgrad=getappdata(resultfig,'vatgrad');
if isempty(vatgrad)
    clear('vatgrad');
end


%% ########################################################################
%% TEMPORARY
% thresh = 1e-8;


%% ________________________________________________________________________
%% FLOWFIELD VISUALIZATION
ea_dispt('Calculating quiver field of gradient for display purposes...');
midpts = vol.ctr;


%% select subset of points to use for quiver representation
reduc=10;
indices = (1:reduc:length(midpts) + round(randn(1)*(reduc/3)))';
indices=unique(indices(2:end-1));
indices(indices==0)=[];
indices(indices>length(midpts))=[];


%% transform to template space if necessary
if options.native==1 && options.orignative==0 % case if we are visualizing in MNI but want to calculate VTA in native space -> now transform back to MNI
    c=midpts';
    [~,anatpresent]=ea_assignpretra(options);
    V=ea_open_vol([options.root,options.patientname,filesep,anatpresent{1}]);  % this function obtains image volume information (header) of the files just found in the previous line
    c=V.mat\[c;ones(1,size(c,2))];  % V.mat is the pose of the image (wrt patient)
    midpts=ea_map_coords(c(1:3,:), ...  % Coordinates mapping
        [options.root,options.patientname,filesep,anatpresent{1}], ...
        [options.root,options.patientname,filesep,'y_ea_inv_normparams.nii'], ...  % this file are generated by SPM
        '')';
    midpts=midpts(:,1:3);
    options.native=options.orignative; % go back to template space
end

if dbg_vis
    figure
    quiver3(midpts(:,1),midpts(:,2),midpts(:,3),gradient(:,1),gradient(:,2),gradient(:,3))
    title('quiver field')
end


%% compute normalized gradient
% define midpoints of quiver field
vatgrad(side).x=midpts(indices,1); vatgrad(side).y=midpts(indices,2); vatgrad(side).z=midpts(indices,3);

gradvis=gradient(indices,:);  % subsample gradient
mag_gradvis=sqrt(sum(gradvis'.^2,1))';  % obtain magnitude of the gradient
nmag_gradvis=mag_gradvis;  % copy to get normalized version
nmag_gradvis(nmag_gradvis>thresh)=thresh;

nmag_gradvis(nmag_gradvis>thresh)=1;  % TEMPORARY: to check visualization (then just DELETE THIS LINE)

nmag_gradvis=(nmag_gradvis-min(nmag_gradvis(:)))/max(nmag_gradvis-min(nmag_gradvis(:)));  % norm from 1 - 0
nmag_gradvis=nmag_gradvis/5; % largest grad vector will be 1/50 mm long

% apply scaling to gradvis
gradvis=gradvis.*repmat(nmag_gradvis,1,3);
gradvis=gradvis./repmat(mag_gradvis,1,3);
vatgrad(side).qx=gradvis(:,1); vatgrad(side).qy=gradvis(:,2); vatgrad(side).qz=gradvis(:,3);

% share normalized gradient with the workflow
setappdata(resultfig,'vatgrad',vatgrad);


%% prepare vat structure
% calculate electric field ET by calculating midpoints of each
% mesh-connection and setting difference of voltage to these points.
vat.pos=midpts;

if dbg_vis && exist('c', 'var')
    figure
    plot3(midpts(:,1),midpts(:,2),midpts(:,3),'b.');
    title('Volume warped to native space')
end

vat.ET = sqrt(sum(gradient'.^2,1)); % vol.cond.*sqrt(sum(gradient'.^2,1)); would be stromstaerke (=> magnitude).

% remove electrode nodes from VAT if specified in the options
if options.prefs.machine.vatsettings.horn_removeElectrode
    where_electrode = unique(vol.tet(el_cond,:));
    vat.pos(where_electrode,:) = [];
    vat.ET(where_electrode) = [];
end

ea_dispt('Preparing VAT...');
vat.tET=vat.ET>thresh;
vat.tpos=vat.pos(vat.tET,:);
% nvat.tpos=nvat.pos(vat.tET,:);
outliers=ea_removeoutliers(vat.tpos,mean(dpvx,1),voltix,constvol);
vat.tpos(outliers,:)=[];
% nvat.tpos(outliers,:)=[];
if dbg_vis
    figure, plot3(vat.tpos(:,1),vat.tpos(:,2),vat.tpos(:,3),'r.');
    [~,trajectory,~,~,~,~]=ea_load_reconstruction(options);
    hold on
    plot3(trajectory{side}(:,1),trajectory{side}(:,2),trajectory{side}(:,3),'k*');
%     plot3(nvat.tpos(:,1),nvat.tpos(:,2),nvat.tpos(:,3),'m.');
    toptions=options; toptions.native=1;
    [~,trajectory,~,~,~,~]=ea_load_reconstruction(toptions);
    plot3(trajectory{side}(:,1),trajectory{side}(:,2),trajectory{side}(:,3),'g*');

end

% the following will be used for volume-to-isosurf creation as well as
% volumetrics of the vat in mm^11
ea_dispt('Calculating interpolant on scattered FEM mesh data...');
F=scatteredInterpolant(vat.pos(:,1),vat.pos(:,2),vat.pos(:,3),vat.ET','linear','none');


%% convert to equispaced image
ea_dispt('Converting to equispaced image data...');
res=100;
gv=cell(3,1); spacing=zeros(3,1);
try
    for dim=1:3  
        % in each dimension create an array of res=100 equispaced points 
        % from the minimum position of the activated area to the maximum one
        gv{dim}=linspace(min(round(vat.tpos(:,dim)))-5,max(round(vat.tpos(:,dim)))+5,res);
        spacing(dim)=abs(gv{dim}(1)-gv{dim}(2));
    end
catch  % if linspace fails it means the vat is empty or that it has problems -> no vat is returned
    varargout{1}=nan;
    varargout{2}=nan;
    varargout{3}=nan;
    return
end


%% ________________________________________________________________________
%%  NIfTI AND DATA EXPORT
ea_dispt('Creating nifti header for export...');

chun1=randperm(res); chun2=randperm(res); chun3=randperm(res);
Vvat.mat=mldivide([(chun1);(chun2);(chun3);ones(1,res)]',[gv{1}(chun1);gv{2}(chun2);gv{3}(chun3);ones(1,res)]')';
Vvat.dim=[res,res,res];
Vvat.dt=[4,0];
Vvat.n=[1 1];
Vvat.descrip='lead dbs - vat';
if ~exist([options.root,options.patientname,filesep,'stimulations',filesep,ea_nt(options)],'file')
    mkdir([options.root,options.patientname,filesep,'stimulations',filesep,ea_nt(options)]);
end

ea_dispt('Filling data with values from interpolant...');
eeg = F(gv);
eeg(isnan(eeg))=0;
eeg(eeg>options.prefs.vat.efieldmax)=options.prefs.vat.efieldmax; % upperlimit files to 10000.

if dbg_vis  % TODO: check this block, since at the moment F is exactly superimposed to vat.pos (interpolant did nothing)
    figure
    plot3(F.Points(:,1),F.Points(:,2),F.Points(:,3),'r.')
    hold on
    plot3(vat.pos(:,1),vat.pos(:,2),vat.pos(:,3),'b.')
    title('E-Field')
end

ea_dispt('Calculating output file data...');
eg=eeg;
eg=eg>thresh;
% binary e-field - "vat"

neeg=eeg;
neeg(~eg)=nan;

neeg(neeg>0)=ea_normal(neeg(neeg>0),1,0,' ',0,1,'TRUE');%
% normalized e-field (zscored).
neeg(~isnan(neeg))=neeg(~isnan(neeg))-min(neeg(~isnan(neeg)));
neeg(~isnan(neeg))=neeg(~isnan(neeg))/sum(neeg(~isnan(neeg))); % 0-1 distributed.

[xg,yg,zg] = meshgrid(gv{1},gv{2},gv{3});

XYZmax=[max(yg(eg>0)),max(xg(eg>0)),max(zg(eg>0))]; % x and y need to be permuted here (should be correct but wouldnt matter anyways since only serves to calc radius)
try
    radius=ea_pdist([XYZmax;dpvx]);
catch
     keyboard
end


%% compute VAT volume
%eg=smooth3(eg,'gaussian',[25 25 25]);
ea_dispt('Calculating volume...');

vatvolume=sum(eg(:))*spacing(1)*spacing(2)*spacing(3); % returns volume of vat in mm^3
S.volume(side)=vatvolume;


%% export files and variables
ea_dispt('Writing files...');

% determine stimulation name:
if ~exist([options.root,options.patientname,filesep,'stimulations',filesep,ea_nt(options),stimname],'file')
    mkdir([options.root,options.patientname,filesep,'stimulations',filesep,ea_nt(options),stimname]);
end

switch side
    case 1
        Vvat.fname=[options.root,options.patientname,filesep,'stimulations',filesep,ea_nt(options),stimname,filesep,'vat_right.nii'];
        Vvate=Vvat; Vvatne=Vvat;
        Vvate.fname=[options.root,options.patientname,filesep,'stimulations',filesep,ea_nt(options),stimname,filesep,'vat_efield_right.nii'];
        Vvatne.fname=[options.root,options.patientname,filesep,'stimulations',filesep,ea_nt(options),stimname,filesep,'vat_efield_gauss_right.nii'];
    case 2
        Vvat.fname=[options.root,options.patientname,filesep,'stimulations',filesep,ea_nt(options),stimname,filesep,'vat_left.nii'];
        Vvate=Vvat; Vvatne=Vvat;
        Vvate.fname=[options.root,options.patientname,filesep,'stimulations',filesep,ea_nt(options),stimname,filesep,'vat_efield_left.nii'];
        Vvatne.fname=[options.root,options.patientname,filesep,'stimulations',filesep,ea_nt(options),stimname,filesep,'vat_efield_gauss_left.nii'];
end
%save(stimfile,'S');
ea_savestimulation(S,options);
% setappdata(lgfigure,'curS',S);

%spm_write_vol(Vvat,flipdim(eg,3));

Vvate.img=eeg; %permute(eeg,[2,1,3]);
Vvate.dt=[16,0];
ea_write_nii(Vvate);

Vvatne.img=neeg; %permute(neeg,[2,1,3]);
ea_write_nii(Vvatne);

Vvat.img=eg; %permute(eg,[1,2,3]);
ea_write_nii(Vvat);

ea_dispt('Calculating isosurface to display...');
vatfv=isosurface(xg,yg,zg,permute(Vvat.img,[2,1,3]),0.75);

caps=isocaps(xg,yg,zg,permute(Vvat.img,[2,1,3]),0.5);

vatfv.faces=[vatfv.faces;caps.faces+size(vatfv.vertices,1)];
vatfv.vertices=[vatfv.vertices;caps.vertices];


try
    vatfv=ea_smoothpatch(vatfv,1,35);
catch
    try
        cd([ea_getearoot,'ext_libs',filesep,'smoothpatch']);
        mex ea_smoothpatch_curvature_double.c -v
        mex ea_smoothpatch_inversedistance_double.c -v
        mex ea_vertex_neighbours_double.c -v
        vatfv=ea_smoothpatch(vatfv);
    catch
        warndlg('Patch could not be smoothed. Please supply a compatible Matlab compiler to smooth VTAs.');
    end
end
% new save by Till to save VAT and quiver in seperate .mat-file for quick
% visualization
switch side
    case 1
        vatfvname=[options.root,options.patientname,filesep,'stimulations',filesep,ea_nt(options),stimname,filesep,'vat_right.mat'];
    case 2
        vatfvname=[options.root,options.patientname,filesep,'stimulations',filesep,ea_nt(options),stimname,filesep,'vat_left.mat'];
end
vatgrad = vatgrad(side);
save(vatfvname,'vatfv','vatgrad','vatvolume');

% vta.nii save, filled and eroded/dilated by 3 voxels.
Vvat.img=imfill(Vvat.img,'holes');
SE = strel('sphere',3);
Vvat.img = imerode(Vvat.img,SE);
Vvat.img = imdilate(Vvat.img,SE);
ea_write_nii(Vvat);


%% ________________________________________________________________________
%% OUTPUT
varargout{1}=vatfv;
varargout{2}=vatvolume;
varargout{3}=radius;
ea_dispt(''); % stop chain of timed processes.


%% 
if dbg; toc; end


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



function potential = ea_apply_dbs(vol,elec,val,unipolar,constvol,boundarynodes)
if constvol
    if unipolar
        dirinodes = [boundarynodes; elec];
    else
        dirinodes = elec;
    end

    rhs = zeros(length(vol.pos),1);
    dirival = zeros(size(vol.pos,1),1);
    dirival(elec) = val(:,1);
else

    if unipolar
        dirinodes = boundarynodes;
    else
        dirinodes = 1;
    end
    dirival = zeros(size(vol.pos,1),1);

    rhs = zeros(size(vol.pos,1),1);
    uvals=unique(val(:,2));
    if unipolar && length(uvals)==1
        elec_center_id = ea_find_elec_center(elec,vol.pos);
        rhs(elec_center_id) = val(1,1);
    else

        for v=1:length(uvals)
            elec_center_id = ea_find_elec_center(elec(val(:,2)==uvals(v)),vol.pos);
            thesevals=val(val(:,2)==uvals(v),1);
            rhs(elec_center_id) = thesevals(1);
        end

        %warning('Bipolar constant current stimulation currently not implemented!');
    end
end

[stiff, rhs] = ea_dbs(vol.stiff,rhs,dirinodes,dirival);

potential = ea_sb_solve(stiff,rhs);
end  % subfunction



% the following 3 functions are used in ea_apply_dbs in order to obtain the
% potential distribution
function center_id = ea_find_elec_center(elec, pos)

center = mean(pos(elec,:));

dist_center = sqrt(sum((pos(elec,:)-repmat(center,length(elec),1)).^2,2));
[~, elec_id] = min(dist_center);
center_id = elec(elec_id);
end  % subfunction



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



function x = ea_sb_solve(sysmat,vecb) %#ok<INUSD>

% SB_SOLVE
%
% $Id: sb_solve.m 8776 2013-11-14 09:04:48Z roboos $
try
    L = ichol(sysmat); %#ok<NASGU>
catch
    alpha = max(sum(abs(sysmat),2)./diag(sysmat))-2;
    L = ichol(sysmat, struct('type','ict','droptol',1e-3,'diagcomp',alpha)); %#ok<NASGU>
end

%scalen
[~,x]=evalc('pcg(sysmat,vecb,10e-10,5000,L,L'',vecb)');
end  % subfunction



% once I've computed the potential distribution with the 3 functions above, 
% I need to compute its gradient
function gradient = ea_calc_gradient(vol,potential)
normal = cross(vol.pos(vol.tet(:,4),:)-vol.pos(vol.tet(:,3),:),vol.pos(vol.tet(:,3),:)-vol.pos(vol.tet(:,2),:));
gradient = repmat(potential(vol.tet(:,1))./sum(normal.*(vol.pos(vol.tet(:,1),:)-(vol.pos(vol.tet(:,2),:)+vol.pos(vol.tet(:,3),:)+vol.pos(vol.tet(:,4),:))/3),2),1,3).*normal;
normal = cross(vol.pos(vol.tet(:,1),:)-vol.pos(vol.tet(:,4),:),vol.pos(vol.tet(:,4),:)-vol.pos(vol.tet(:,3),:));
gradient = gradient + repmat(potential(vol.tet(:,2))./sum(normal.*(vol.pos(vol.tet(:,2),:)-(vol.pos(vol.tet(:,3),:)+vol.pos(vol.tet(:,4),:)+vol.pos(vol.tet(:,1),:))/3),2),1,3).*normal;
normal = cross(vol.pos(vol.tet(:,2),:)-vol.pos(vol.tet(:,1),:),vol.pos(vol.tet(:,1),:)-vol.pos(vol.tet(:,4),:));
gradient = gradient + repmat(potential(vol.tet(:,3))./sum(normal.*(vol.pos(vol.tet(:,3),:)-(vol.pos(vol.tet(:,4),:)+vol.pos(vol.tet(:,1),:)+vol.pos(vol.tet(:,2),:))/3),2),1,3).*normal;
normal = cross(vol.pos(vol.tet(:,3),:)-vol.pos(vol.tet(:,2),:),vol.pos(vol.tet(:,2),:)-vol.pos(vol.tet(:,1),:));
gradient = gradient + repmat(potential(vol.tet(:,4))./sum(normal.*(vol.pos(vol.tet(:,4),:)-(vol.pos(vol.tet(:,1),:)+vol.pos(vol.tet(:,2),:)+vol.pos(vol.tet(:,3),:))/3),2),1,3).*normal;
end  % subfunction



% the following 2 equations are used to get rid of outlier results in the
% data visualization step
function outliers=ea_removeoutliers(pointcloud,mp,voltix,constvol)
% using the heuristic Maedler/Coenen model to detect outliers
vmax=max(abs(voltix));

if ~constvol
   vmax=vmax*1000;
end
r=maedler12_eq3(vmax,1000);
r=r*4.5;

D=pointcloud-repmat(mp,size(pointcloud,1),1);
S=[r,r,r];  % std(D,[],1);  % fixed 50 cm
outliers=D>repmat(S,size(D,1),1);
outliers=any(outliers,2);  % for each row see if there is at least an outlier 
                           % coordinate and in case mark that row (point in 
                           % the point cloud) as outlier
end  % subfunction



function r=maedler12_eq3(U,Im)  
% INPUT
% U: stimulation voltage
% Im: impedance
% 
% OUTPUT
% r: estimated VAT radius according to Maedler and Coenen 2012
% 
% This function radius of Volume of Activated Tissue for stimulation settings U and Ohm. See Maedler 2012 for details.
% Clinical measurements of DBS electrode impedance typically range from
% 500?1500 Ohm (Butson 2006).
r=0; %
if U %(U>0)

    k1=-1.0473;
    k3=0.2786;
    k4=0.0009856;

    r=-(k4*Im-sqrt(k4^2*Im^2  +   2*k1*k4*Im    +   k1^2 +   4*k3*U)   +   k1)/(2*k3);
end
end  % subfunction