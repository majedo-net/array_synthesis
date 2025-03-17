close all; clear;
labels = ["Dipole Z 37", "Dipole Y 37","Dipole Y 67","Dipole Z 67","Spiral 37"];
paths = ["dipole/z","dipole/y","dipole/y/NP5","dipole/z/NP5","spiral"];
elements=[37,37,67,67,37];
nearest_neighbors = {};
errors = {};
plot_markers = ["ko-","bx-.","r^:","c*-.","m->"];


for idx=1:numel(labels)
    if elements(idx) == 37
        nearest_neighbors{end+1} = [4,6,8,10,12,14,16,18,20,24,28,32,36];
    elseif elements(idx) == 67
        nearest_neighbors{end+1} = [4,6,8,10,12,14,16,18,20,24,28,32,36,40,44,48,52,58,62,66];
    end
    errors{end+1} = zeros(numel(nearest_neighbors{end}),1);
    for edx=0:elements(idx)-1
        epf = h5read(sprintf('results/aws-tess-array/%s/full/element%i.hdf5',paths(idx),edx),'/eph');
        epf = epf.r + 1j.*epf.i;
        etf = h5read(sprintf('results/aws-tess-array/%s/full/element%i.hdf5',paths(idx),edx),'/eth');
        etf = etf.r + 1j.*etf.i;
        e_full = abs(epf+etf)./max(abs(epf+etf),[],"all");
        for nid=1:numel(nearest_neighbors{end})
            epf = h5read(sprintf('results/aws-tess-array/%s/nr%i/element%i.hdf5',paths(idx),nearest_neighbors{end}(nid),edx),'/eph');
            epf = epf.r + 1j.*epf.i;
            etf = h5read(sprintf('results/aws-tess-array/%s/nr%i/element%i.hdf5',paths(idx),nearest_neighbors{end}(nid),edx),'/eth');
            etf = etf.r + 1j.*etf.i;
            e_tess = abs(epf+etf)./max(abs(epf+etf),[],"all");
            errors{end}(nid) = errors{end}(nid) + norm(abs(e_full-e_tess),'fro')./elements(idx);
        end
    end
    plot(nearest_neighbors{idx},errors{idx},plot_markers(idx)); hold on;
end
grid minor;
legend(labels);