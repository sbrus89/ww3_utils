% Example_7_Gvlobal: Make a global mesh
clearvars; clc;
addpath(genpath('utilities/'));
addpath(genpath('datasets/'));
addpath(genpath('m_map/'));

%% STEP 1: set mesh extents and set parameters for mesh. 
%% The greater US East Coast and Gulf of Mexico region
min_el    = 2e3;  	             % minimum resolution in meters.
bbox      = [-180 180; -88  85]; % lon min lon max; lat min lat max
max_el    = 50e3; 		         % maximum resolution in meters. 
wl        = 30;                  % 30 elements resolve M2 wavelength.
dt        = 0;                   % Only reduces res away from coast
grade     = 0.0;                % mesh grade in decimal percent. 
R         = 3; 			         % Number of elements to resolve feature.
slp       = 10;                  % slope of 10

outname = 'Global_DE_2km'

%% STEP 2: specify geographical datasets and process the geographical data
%% to be used later with other OceanMesh classes...
dem       = 'bathy_SSG_1_120.nc'; 
coastline = {'GSHHS_h_L1','GSHHS_h_L6'};
%gdat1 = geodata('shp',coastline,'dem',dem,...
%                'bbox',bbox,'h0',min_el);
gdat1 = geodata('shp',coastline,'dem',dem,...
                'bbox',bbox,'h0',min_el);
            
%% STEP 3: create an edge function class
%fh1 = edgefx('geodata',gdat1,...
%             'fs',R,'wl',wl,'max_el',max_el,...
%             'slp',slp,'dt',dt,'g',grade); 

fh1 = edgefx('geodata',gdat1,'g',grade,...
             'usr','cellWidthVsLatLon_2km.nc'); 
                
%% STEP 4: Pass your edgefx class object along with some meshing options 
%% and build the mesh...
mshopts = meshgen('ef',fh1,'bou',gdat1,'nscreen',10,'plot_on',1,...
                  'proj','stereo','cleanup',0);
mshopts = mshopts.build; 

%% STEP 5: Match points and edges across boundary
m = mshopts.grd; % get out the msh object 
m = clean(m,'ds',0,'db',0.5);

% Plotting the triangulation on Robinson
plot(m,'tri',1,'Robinson');
save([outname '.mat'],'m'); 
write(m,outname);
exit;
