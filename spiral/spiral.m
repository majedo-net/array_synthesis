function spiral
addpath('/opt/openEMS/share/openEMS/matlab');
addpath('/opt/openEMS/share/CSXCAD/matlab');
args = argv();
%disp(size(args));
f_start = str2num(args{1,1});
f_stop = str2num(args{2,1});
r0 = str2num(args{3,1}); % minimum inner spiral radius
alpha = str2num(args{4,1}); % exponential spiral coefficient
h = str2num(args{5,1}); % circular slot radius
rmax = str2num(args{6,1}); % back length
nP = str2num(args{7,1}); % particle index
nF = args{8,1}; % frequency index
%% setup the simulation
physical_constants;
unit = 1e-3; % all length in mm
N = log(rmax/r0)/(2*pi*alpha);
nP
nF
% frequency range of interest
padding = (3e8/f_start)/(2*unit);
max_res = floor(c0 / (f_stop) / unit / 15);% cell size: lambda/15
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
%mesh = DetectEdges(CSX, mesh,'SetPropertyType','Excitation');
%mesh = SmoothMesh(mesh,[max_res max_ang max_res],1.4);
CSX = DefineRectGrid(CSX, unit, mesh);

%% nf2ff calc
start = [mesh.x(2) mesh.y(2) mesh.z(2)];
stop  = [mesh.x(end-8) mesh.y(end-8) mesh.z(end-8)];
[CSX, nf2ff] = CreateNF2FFBox(CSX, 'nf2ff', start, stop, 'Directions', [1 1 1 1 0 1]);

%% prepare simulation folder
Sim_Path = sprintf('tmp_spiral%d%c',nP,nF);
Sim_CSX = 'spiral_ant.xml';

[status, message, messageid] = rmdir( Sim_Path, 's' ); % clear previous directory
[status, message, messageid] = mkdir( Sim_Path ); % create empty simulation folder

%% write openEMS compatible xml-file
WriteOpenEMS([Sim_Path '/' Sim_CSX], FDTD, CSX);

%% show the structure
%CSXGeomPlot([Sim_Path '/' Sim_CSX]);

%% run openEMS
Settings.Silent = 0;
RunOpenEMS(Sim_Path, Sim_CSX,'--numThreads=5',Settings);

%% postprocessing & do the plots
freq = linspace(f_start,f_stop,20);

port = calcPort(port, Sim_Path, freq);
Zin = port.uf.tot ./ port.if.tot;
s11 = port.uf.ref ./ port.uf.inc;


% calculate the far field at phi=0 degrees and at phi=90 degrees
disp( 'calculating far field' );
thetas = 0;
phis = 0;
f0 = linspace(f_start,f_stop,20);
nf2ff = CalcNF2FF(nf2ff, Sim_Path, f0, thetas*pi/180, phis*pi/180,'Verbose',0,'Radius',(3e8/f_start)*10);
R_g = 10.*log10(transpose(nf2ff.Dmax.'.*(1-abs(s11))));
G = 10*log10(nf2ff.Dmax);
Zmag = abs(Zin);
Zmag = Zmag';
%save('Zmag', 'Zmag')
Zave = mean(Zmag);
Zstd = std(Zmag);
Gave = mean(G);
Gstd = std(G);
S11 = abs(s11');

%sim_output = [Zstd, Zave, Gave, Gstd];

Zreal = real(Zin);
Zimag = imag(Zin);

Zreal_std = std(Zreal);
Zimag_std = std(Zimag);
Zreal_ave = mean(Zreal);
Zimag_ave = mean(Zimag);
Gmin = min(G);

sim_output = [Zreal_std, Zimag_std, Zreal_ave, Zimag_ave, Gmin,nP];
GZresults = [freq', Zreal', Zimag', G, S11];
csvwrite(sprintf('results/GZresults%d%c.csv', nP,nF), GZresults);
[status, message, messageid] = rmdir( Sim_Path, 's' ); % clear previous directory
