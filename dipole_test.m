
close all
clear
clc

%% switches & options...
postprocessing_only = 0;

%% prepare simulation folder
Sim_Path = 'tmp_Dipole_SAR';
Sim_CSX = 'Dipole_SAR.xml';

%% setup the simulation
physical_constants;
unit = 1e-3; % all lengths in mm

feed.R = 70; % feed resistance

%% setup FDTD parameter & excitation function
f0 = 6e9; % center frequency
lambda0 = c0/f0;

f_stop = 6.5e9; % 20 dB corner frequency
lambda_min = c0/f_stop;

mesh_res_air = lambda_min/20/unit;
mesh_res_phantom = 0.5;

dipole_length = 0.435*lambda0/unit;
disp(['Lambda-half dipole length: ' num2str(dipole_length) 'mm'])

%%
FDTD = InitFDTD();
FDTD = SetGaussExcite( FDTD, 0, f_stop );
% apply PML-8 boundary conditions in all directions
BC = {'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8'};
FDTD = SetBoundaryCond( FDTD, BC );

%% setup CSXCAD geometry & mesh
CSX = InitCSX();

%% Dipole
CSX = AddMetal( CSX, 'Dipole' ); % create a perfect electric conductor (PEC)
CSX = AddBox(CSX, 'Dipole', 1, [0 0 -dipole_length/2], [0 0 dipole_length/2]);

% mesh lines for the dipole
mesh.x = 0;
mesh.y = 0;
mesh.z = [-dipole_length/2-[-1/3 2/3]*mesh_res_phantom dipole_length/2+[-1/3 2/3]*mesh_res_phantom];

%% apply the excitation & resist as a current source
[CSX, port] = AddLumpedPort(CSX, 100, 1, feed.R, [-0.1 -0.1 -mesh_res_phantom/2], [0.1 0.1 +mesh_res_phantom/2], [0 0 1], true);

% mesh lines for the port
mesh.z = [mesh.z -mesh_res_phantom/2 +mesh_res_phantom/2];

%% smooth the mesh over the dipole and phantom
mesh = SmoothMesh(mesh, mesh_res_phantom);

%% add lines for the air-box
mesh.x = [mesh.x -50 50];
mesh.y = [mesh.y -50 50];
mesh.z = [mesh.z -50 50];

% smooth the final mesh (incl. air box)
mesh = SmoothMesh(mesh, mesh_res_air);


%% nf2ff calc
start = [mesh.x(1)   mesh.y(1)   mesh.z(1)];
stop  = [mesh.x(end) mesh.y(end) mesh.z(end)];
[CSX nf2ff] = CreateNF2FFBox(CSX, 'nf2ff', start, stop, 'OptResolution', lambda_min/15/unit);

%%
% add 10 equidistant cells (air)
% around the structure to keep the pml away from the nf2ff box
mesh = AddPML( mesh, 10 );

% Define the mesh
CSX = DefineRectGrid(CSX, unit, mesh);

%%
if (postprocessing_only==0)
    [status, message, messageid] = rmdir( Sim_Path, 's' ); % clear previous directory
    [status, message, messageid] = mkdir( Sim_Path ); % create empty simulation folder

    % write openEMS compatible xml-file
    WriteOpenEMS( [Sim_Path '/' Sim_CSX], FDTD, CSX );

    % show the structure
    CSXGeomPlot( [Sim_Path '/' Sim_CSX] );

    % run openEMS
    RunOpenEMS( Sim_Path, Sim_CSX );
end


%% postprocessing & make the plots
freq = linspace(5500e6, 6500e6, 501 );
port = calcPort(port, Sim_Path, freq);

s11 = port.uf.ref./port.uf.inc;
Zin = port.uf.tot./port.if.tot;

Pin_f0 = interp1(freq, port.P_acc, f0);

%%
% plot feed point impedance
figure
plot( freq/1e6, real(Zin), 'k-', 'Linewidth', 2 );
hold on
grid on
plot( freq/1e6, imag(Zin), 'r--', 'Linewidth', 2 );
title( 'feed point impedance' );
xlabel( 'frequency f / MHz' );
ylabel( 'impedance Z_{in} / Ohm' );
legend( 'real', 'imag' );

% plot reflection coefficient S11
figure
plot( freq/1e9, 20*log10(abs(s11)), 'k-', 'Linewidth', 2 );
grid on
title( 'reflection coefficient' );
xlabel( 'frequency f / MHz' );
ylabel( 'S_{11} (dB)' );