close all; clear;

centers =[0, 0;
    100 100;
    100 -100;
    -100 -100;
    -100 100];
rs = [5,90;
    5,90;
    5,50;
    5,90];

%% setup the simulation
physical_constants;
freq = 2e9;
unit = 1e-3; % all length in mm
f_start = 0.9 * freq;
f_stop = freq;
max_res = floor(c0 / (f_stop) / unit / 20);% cell size: lambda/30
padding = max_res *20;
hs = 1.6; % substrate thickness
h = 30; % cavity height

Np = 1; % simulation index

% size of the simulation box
SimBox = [padding+max(rs)*2+max(centers) padding+max(rs)*2+max(centers) padding+h+hs];

%% setup FDTD parameter & excitation function
FDTD = InitFDTD('EndCriteria', 1e-4,'NrTs',50000);
FDTD = SetGaussExcite(FDTD,0.5*(f_start+f_stop),0.5*(f_stop-f_start));
BC = {'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PEC' 'PML_8'}; % boundary conditions
FDTD = SetBoundaryCond( FDTD, BC );

%% setup CSXCAD geometry & mesh
% currently, openEMS cannot automatically generate a mesh

CSX = InitCSX();

%create fixed lines for the simulation box and port
mesh.x = [-SimBox(1)/2 SimBox(1)/2];
mesh.y = [-SimBox(2)/2 SimBox(2)/2];
mesh.z = [0 h-hs h+hs SimBox(3)];

%% Generate spirals
for idx=1:size(rs)
    this_r = rs(idx,:);
    this_center = centers(idx,:);
    if idx==1
        excite=true;
    else
        excite=false;
    end
    [CSX,mesh] = make_spiral(CSX,mesh,this_center, ...
        this_r,0.32,h,hs,idx,excite);
end

mesh.x = SmoothMeshLines(mesh.x, max_res, 1.4); % create a smooth mesh between specified fixed mesh lines
mesh.y = SmoothMeshLines(mesh.y, max_res, 1.4 );
mesh.z = SmoothMeshLines(mesh.z,max_res,1.4);

%% Detect edges and smooth mesh
mesh = AddPML(mesh,8);
CSX = DefineRectGrid(CSX, unit, mesh);

%% nf2ff dump box
start = [mesh.x(9) mesh.y(9) mesh.z(9)];
stop  = [mesh.x(end-9) mesh.y(end-9) mesh.z(end-9)];
[CSX, nf2ff] = CreateNF2FFBox(CSX, 'nf2ff', start, stop, 'Directions', [1 1 1 1 0 1]);

%% prepare simulation folder
Sim_Path = sprintf('tmp_spiral%d', Np);
Sim_CSX = 'spiral_ant.xml';

[status, message, messageid] = rmdir( Sim_Path, 's' ); % clear previous directory
[status, message, messageid] = mkdir( Sim_Path ); % create empty simulation folder

%% write openEMS compatible xml-file
WriteOpenEMS([Sim_Path '/' Sim_CSX], FDTD, CSX);

%% show the structure
CSXGeomPlot([Sim_Path '/' Sim_CSX]);

rmdir( Sim_Path, 's' ); % clear previous directory