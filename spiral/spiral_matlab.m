clear;
%disp(size(args));
freq = 0.5e9;
r0 = 5; % minimum inner spiral radius
alpha = 0.32; % exponential spiral coefficient
h = 20; % height above gorund plane
rmax = 70; % max radius
nP = 0; % simeoultaneous sim index
r0
rmax
%% setup the simulation
physical_constants;
unit = 1e-3; % all length in mm
f_start = 0.9 * freq;
f_stop = freq;
N = log(rmax/r0)/(2*pi*alpha)
max_res = floor(c0 / (f_stop) / unit / 20)% cell size: lambda/30
padding = max_res *20
nP
hs = 1.6; % substrate thickness

%% Calculate top spiral
a0 = pi/4;
a1 = linspace(a0,N*2*pi,100*N);
a2 = a1 + pi/2;
a3 = a2 + pi/2;
a4 = a3 + pi/2;
r = r0 .* exp(alpha.*a1);
tp = zeros([2 2*numel(r)+52]);
bp = tp;
length(r)
tp(1,:) = [0,r,r(end).*ones([1 50]),flip(r),0];
tp(2,:) = [0,a1,linspace(a1(end),a2(end),50),flip(a2),0];
bp(1,:) = [0,r,r(end).*ones([1 50]),flip(r),0];
bp(2,:) = [0,a3,linspace(a3(end),a4(end),50),flip(a4),0];
% size of the simulation box
SimBox = [padding+max(r)*2 padding+max(r)*2 padding+h+hs];

%% setup FDTD parameter & excitation function
FDTD = InitFDTD('EndCriteria', 1e-4,'NrTs',50000);
FDTD = SetGaussExcite(FDTD,0.5*(f_start+f_stop),0.5*(f_stop-f_start));
BC = {'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PEC' 'PML_8'}; % boundary conditions
FDTD = SetBoundaryCond( FDTD, BC );

%% setup CSXCAD geometry & mesh
% currently, openEMS cannot automatically generate a mesh

CSX = InitCSX();

%create fixed lines for the simulation box and port
mesh.x = [-SimBox(1)/2 -1*r0 -1 1 r0 SimBox(1)/2];
mesh.x = SmoothMeshLines(mesh.x, max_res, 1.4); % create a smooth mesh between specified fixed mesh lines

mesh.y = [-SimBox(2)/2 -1*r0 -1 1 r0 SimBox(2)/2];
mesh.y = SmoothMeshLines( mesh.y, max_res, 1.4 );

%create fixed lines for the simulation box and given number of lines inside the substrate
mesh.z = [0 h-hs h+hs SimBox(3)];
mesh.z = SmoothMeshLines( mesh.z, max_res, 1.4 );


tp_cart(1,:) = tp(1,:) .* cos(tp(2,:));
tp_cart(2,:) = tp(1,:) .* sin(tp(2,:));
CSX = AddMetal(CSX, 'top_spiral');
CSX = AddLinPoly( CSX, 'top_spiral',10,'z',h+hs,tp_cart,0,'CoordSystem',0);
bp_cart(1,:) = bp(1,:) .* cos(bp(2,:));
bp_cart(2,:) = bp(1,:) .* sin(bp(2,:));
CSX = AddMetal(CSX, 'bot_spiral');
CSX = AddLinPoly( CSX, 'bot_spiral',10,'z',h-hs,bp_cart,0,'CoordSystem',0);

CSX = AddMaterial(CSX,'substrate');
CSX = SetMaterialProperty(CSX,'substrate','Epsilon',4.2);
start = [0 0 h-hs];
stop = [max(r) 2*pi h+hs];
CSX = AddBox(CSX, 'substrate', 0, start, stop,'CoordSystem',1);


%% apply the excitation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

start=[-1 -1 h-hs];
stop =[1 1 h+hs];
[CSX,port] = AddLumpedPort(CSX,5,1,50,start,stop,[0 0 1],true);

%% Detect edges and smooth mesh
mesh = AddPML(mesh,8);
CSX = DefineRectGrid(CSX, unit, mesh);

%% nf2ff dump box
start = [mesh.x(9) mesh.y(9) mesh.z(9)];
stop  = [mesh.x(end-9) mesh.y(end-9) mesh.z(end-9)];
[CSX, nf2ff] = CreateNF2FFBox(CSX, 'nf2ff', start, stop, 'Directions', [1 1 1 1 0 1]);

%% prepare simulation folder
Sim_Path = sprintf('tmp_spiral%d',nP);
Sim_CSX = 'spiral_ant.xml';

[status, message, messageid] = rmdir( Sim_Path, 's' ); % clear previous directory
[status, message, messageid] = mkdir( Sim_Path ); % create empty simulation folder

%% write openEMS compatible xml-file
WriteOpenEMS([Sim_Path '/' Sim_CSX], FDTD, CSX);

%% show the structure
CSXGeomPlot([Sim_Path '/' Sim_CSX]);

%% run openEMS
Settings.Silent = 0;
RunOpenEMS(Sim_Path, Sim_CSX);

% calculate the far field at phi=0 degrees and at phi=90 degrees
disp( 'calculating far field' );
thetas = linspace(-pi/2, pi/2, 181);
phis = linspace(0,pi,91);
nf2ff = CalcNF2FF(nf2ff, Sim_Path, freq, thetas, phis,'Verbose',1);
G = cell2mat(nf2ff.E_norm);
csvwrite(sprintf('farfieldspiral_rad_%d_freq_%d.csv', rmax,freq/1e6), G);
[el,az] = meshgrid(thetas,phis);
[X,Y,Z] = sph2cart(az,(el-pi/2),G');
surf(X,Y,Z);
[status, message, messageid] = rmdir( Sim_Path, 's' ); % clear previous directory
