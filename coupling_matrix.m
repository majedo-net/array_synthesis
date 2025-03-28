clear; close all;
%labels = ["Dipole Z 37", "Dipole Y 37","Dipole Y 67","Dipole Z 67","Spiral 37"];
labels = ["Z-Aligned Dipole","Y-Aligned Dipole"];
%paths = ["dipole/z","dipole/y","dipole/y/NP5","dipole/z/NP5","spiral"];
paths = ["dipole/z","dipole/y"];
%elements=[37,37,67,67,37];
elements=[37,37];
nearest_neighbors = {};
errors = {};
plot_markers = ["ko-","bx-.","r^:","c*-.","m->"];

times = [
    [0.43 0.43 0.43 0.43 0.57 0.57 0.71 0.71 0.71 0.86 1 0.86 1.14 1] % z dipole
    [0.33 0.33 0.33 0.44 0.56 0.67 0.56 0.78 0.78 0.78 1 0.89 1.11 1] % y dipole
    ];

for idx=1:numel(labels)
    fprintf('%s\n',labels(idx));
    if elements(idx) == 37
        nearest_neighbors{end+1} = [2,4,6,8,10,12,14,16,18,20,24,28,32,36];
    elseif elements(idx) == 67
        nearest_neighbors{end+1} = [4,6,8,10,12,14,16,18,20,24,28,32,36,40,44,48,52,58,62,66];
    end
    errors{end+1} = zeros(numel(nearest_neighbors{end}),1);

    smn = zeros([elements(idx),301,elements(idx)]);
    smnf = zeros([elements(idx),301,elements(idx)]);
    error = zeros(numel(nearest_neighbors{end}),1);
    s11_error = zeros(numel(nearest_neighbors{end}),1);
    for nid=1:numel(nearest_neighbors{end})
        nr = nearest_neighbors{end}(nid);
        for edx = linspace(0,elements(idx)-1,elements(idx))
            sp = h5read(sprintf('results/aws-tess-array/%s/full/element%i.hdf5',paths(idx),edx),'/smn');
            sps = sp.r + 1j.*sp.i;
            if any(abs(sps)>1,'all')
                fprintf('SP error FULL: %i \n',edx);
            end
            neighb_id = h5read(sprintf('results/aws-tess-array/%s/full/element%i.hdf5',paths(idx),edx),'/nids');
            smnf(edx+1,:,neighb_id+1) = sps;
            clear sps;
            sps = h5read(sprintf('results/aws-tess-array/%s/nr%i/element%i.hdf5',paths(idx),nr,edx),'/smn');
            sps = sps.r + 1j.*sps.i;
            smn_old = sps;
            if any(abs(sps)>1,'all')
                fprintf('SP error TESS: %i \n',edx);
            end
            neighb_id = h5read(sprintf('results/aws-tess-array/%s/nr%i/element%i.hdf5',paths(idx),nr,edx),'/nids');
            smn(edx+1,:,neighb_id+1) = sps;
            clear sps;
        end
        error(nid) = sum(abs(abs(smnf(:,151,:))-abs(smn(:,151,:))),'all')./sum(abs(smnf(:,151,:)),'all');
        error(nid) = error(nid) ./ elements(idx).^2;
        s11_error(nid) = sum(abs(abs(diag(squeeze(smnf(:,151,:)))) - abs(diag(squeeze(smn(:,151,:))))),'all')./sum(abs(diag(squeeze(smnf(:,151,:)))),'all');
        s11_error(nid) = s11_error(nid) ./ elements(idx);
        delta = abs(abs(smnf(:,151,:)) - abs(smn(:,151,:)));
        delta = squeeze(delta);
        %heatmap(delta);

    end
    %plot errors here
    disp('Plots');
    figure(idx);
    grid minor;
    xlabel('Nearest Neighbors N_r');
    yyaxis left;
    plot(nearest_neighbors{end},error,'k-o');
    hold on;
    %plot(nearest_neighbors{end},s11_error,'k-.s');
    
    ax = gca;
    ax.YColor = 'k';
    ylabel('Error \epsilon');
    yyaxis right;
    plot(nearest_neighbors{end},times(idx,:),'b-^');
    %legend('Full','Computation Time');
    ax = gca;
    ax.YColor = 'b';
    ylabel('Computation Time/t(N_r = 36)');
    title(labels{idx});
%     figure((2+(idx-1)*2));
%     title(labels(idx));
%     delta = abs(smnf(:,151,:)) - abs(smn(:,151,:));
%     delta = squeeze(delta);
%     heatmap(delta);

end
