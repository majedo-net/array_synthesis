clear; close all;

%% Z Dipole
% load single
eth = h5read('results/aws-tess-array/ct/single_z.hdf5','/eth');
eth = eth.r + 1j.*eth.i;
eph = h5read('results/aws-tess-array/ct/single_z.hdf5','/eph');
eph = eph.r + 1j.*eph.i;
enorm_single = sqrt(abs(eth).^2 + abs(eph).^2);
%enorm_single = enorm_single./max(enorm_single,[],'all');

dl = [];
dfz = [];

for ii = 1:(numel(dir('results/aws-tess-array/ct/z'))-2)/2
    fc = h5readatt(sprintf('results/aws-tess-array/ct/z/fc%d.hdf5',ii-1),'/','fc');
    l = 3e8/fc/2;
    l0 = 3e8/6e9;
    dl(end+1) = l/l0;
    eth = h5read(sprintf('results/aws-tess-array/ct/z/fc%d.hdf5',ii-1),'/eth');
    eth = eth.r + 1j.*eth.i;
    eph = h5read(sprintf('results/aws-tess-array/ct/z/fc%d.hdf5',ii-1),'/eph');
    eph = eph.r + 1j.*eph.i;
    enorm = sqrt(abs(eth).^2 + abs(eph).^2);
    %enorm = enorm./max(enorm,[],'all');
    dfz(end+1) = norm(enorm-enorm_single,"fro");
end


%% Y Dipole
% load single
eth = h5read('results/aws-tess-array/ct/single_y.hdf5','/eth');
eth = eth.r + 1j.*eth.i;
eph = h5read('results/aws-tess-array/ct/single_y.hdf5','/eph');
eph = eph.r + 1j.*eph.i;
enorm_single = sqrt(abs(eth).^2 + abs(eph).^2);
%enorm_single = enorm_single./max(enorm_single,[],'all');

dl = [];
dfy = [];
for ii = 1:(numel(dir('results/aws-tess-array/ct/y'))-2)/2
    fc = h5readatt(sprintf('results/aws-tess-array/ct/y/fc%d.hdf5',ii-1),'/','fc');
    l = 3e8/fc/2;
    l0 = 3e8/6e9;
    dl(end+1) = l/l0;
    eth = h5read(sprintf('results/aws-tess-array/ct/y/fc%d.hdf5',ii-1),'/eth');
    eth = eth.r + 1j.*eth.i;
    eph = h5read(sprintf('results/aws-tess-array/ct/y/fc%d.hdf5',ii-1),'/eph');
    eph = eph.r + 1j.*eph.i;
    enorm = sqrt(abs(eth).^2 + abs(eph).^2);
    %enorm = enorm./max(enorm,[],'all');
    dfy(end+1) = norm(enorm-enorm_single,'fro');
end

%% Patch
% load single
eth = h5read('results/aws-tess-array/ct/single_patch.hdf5','/eth');
eth = eth.r + 1j.*eth.i;
eph = h5read('results/aws-tess-array/ct/single_patch.hdf5','/eph');
eph = eph.r + 1j.*eph.i;
enorm_single = sqrt(abs(eth).^2 + abs(eph).^2);
%enorm_single = enorm_single./max(enorm_single,[],'all');

dlp = [];
dfp = [];
for ii = 1:(numel(dir('results/aws-tess-array/ct/patch'))-2)/2
    fc = h5readatt(sprintf('results/aws-tess-array/ct/patch/fc%d.hdf5',ii-1),'/','fc');
    l = 3e8/fc/2;
    l0 = 3e8/6e9;
    dlp(end+1) = l/l0;
    eth = h5read(sprintf('results/aws-tess-array/ct/patch/fc%d.hdf5',ii-1),'/eth');
    eth = eth.r + 1j.*eth.i;
    eph = h5read(sprintf('results/aws-tess-array/ct/patch/fc%d.hdf5',ii-1),'/eph');
    eph = eph.r + 1j.*eph.i;
    enorm = sqrt(abs(eth).^2 + abs(eph).^2);
    %enorm = enorm./max(enorm,[],'all');
    dfp(end+1) = norm(enorm-enorm_single,'fro');
end

plot(dl, 10.*log10(dfy),'b^-.'); 
hold on;
plot(dl,10.*log10(dfz),'ro-');
plot(dlp,10.*log10(dfp),'k<-.');
grid minor;
legend('Dipole Y','Dipole Z','Patch');
xlabel('$\frac{d}{\lambda_0}$','Interpreter','latex');
ylabel('$10log(|E_s - E_c|_F)$','Interpreter','latex');