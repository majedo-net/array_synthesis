clear; close all;
labels = ["Dipole Z","Dipole Y","Spiral","Patch"];
paths = ["z","y","spiral","patch"];
plot_markers = ["ko-","bx-.","r^:","c*-.","m->"];
d_over_ls = {};
errors = {};

for idx=1:numel(labels)
    % load single
    eth = h5read(sprintf('results/aws-tess-array/ct/single_%s.hdf5',paths(idx)),'/eth');
    eth = eth.r + 1j.*eth.i;
    eph = h5read(sprintf('results/aws-tess-array/ct/single_%s.hdf5',paths(idx)),'/eph');
    eph = eph.r + 1j.*eph.i;
    enorm_single = sqrt(abs(eth).^2 + abs(eph).^2);
    enorm_single = enorm_single./max(enorm_single,[],'all'); 
    d_over_ls{end+1} = [];
    errors{end+1} = [];
    
    for ii = 1:(numel(dir('results/aws-tess-array/ct/z'))-2)/2
        fc = h5readatt(sprintf('results/aws-tess-array/ct/%s/fc%d.hdf5',paths(idx),ii-1),'/','fc');
        l = 3e8/fc/2;
        l0 = 3e8/6e9;
        d_over_ls{end}(end+1) = l/l0;
        eth = h5read(sprintf('results/aws-tess-array/ct/%s/fc%d.hdf5',paths(idx),ii-1),'/eth');
        eth = eth.r + 1j.*eth.i;
        eph = h5read(sprintf('results/aws-tess-array/ct/%s/fc%d.hdf5',paths(idx),ii-1),'/eph');
        eph = eph.r + 1j.*eph.i;
        enorm = sqrt(abs(eth).^2 + abs(eph).^2);
        enorm = enorm./max(enorm,[],'all');
        errors{end}(end+1) = norm(enorm-enorm_single,"fro");
    end
    plot(d_over_ls{end}, log10(errors{end}),plot_markers(idx)); 
    hold on;
end
legend(labels);
grid minor;
xlabel('$\frac{d}{\lambda_0}$','Interpreter','latex');
ylabel('$\log(|E_s - E_c|_F)$','Interpreter','latex');
