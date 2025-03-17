clear; close all;
%labels = ["Dipole Z 37", "Dipole Y 37","Dipole Y 67","Dipole Z 67","Spiral 37"];
labels = ["Spiral 37"];
%paths = ["dipole/z","dipole/y","dipole/y/NP5","dipole/z/NP5","spiral"];
paths = ["spiral"];
%elements=[37,37,67,67,37];
elements=[37];
nearest_neighbors = {};
errors = {};
plot_markers = ["ko-","bx-.","r^:","c*-.","m->"];


for idx=1:numel(labels)
    fprintf('%s\n',labels(idx));
    if elements(idx) == 37
        nearest_neighbors{end+1} = [4,6,8,10,12,14,16,18,20,24,28,32,36];
    elseif elements(idx) == 67
        nearest_neighbors{end+1} = [4,6,8,10,12,14,16,18,20,24,28,32,36,40,44,48,52,58,62,66];
    end
    errors{end+1} = zeros(numel(nearest_neighbors{end}),1);

    smn = zeros([elements(idx),301,elements(idx)]);
    smnf = zeros([elements(idx),301,elements(idx)]);

    for nid=1:numel(nearest_neighbors{end})
        nr = nearest_neighbors{end}(nid);
        for edx = linspace(0,elements(idx)-1,elements(idx))
            sp = h5read(sprintf('results/aws-tess-array/%s/full/element%i.hdf5',paths(idx),edx),'/smn');
            sps = sp.r + 1j.*sp.i;
            smnf(edx+1,:,:) = sps;
            sps = h5read(sprintf('results/aws-tess-array/%s/nr%i/element%i.hdf5',paths(idx),nr,edx),'/smn');
            sps = sps.r + 1j.*sps.i;
            neighb_id = h5read(sprintf('results/aws-tess-array/%s/nr%i/element%i.hdf5',paths(idx),nr,edx),'/nids');
            smn(edx+1,:,neighb_id+1) = sps;
        end
    end
    % plot errors here
    disp('Plots');
    figure(idx);
    title(labels(idx));
    delta = abs(smnf(:,151,:)) - abs(smn(:,151,:));
    delta = squeeze(delta);
    heatmap(delta);

end
% %% Dipole Z
% smn = zeros([37,301,37]);
% smnf = zeros([37,301,37]);
% nrsd = [4,6,8,10,12,14,16,18,20,24,28,32,36];
% nfd = [];
% nfd_error = zeros([numel(nrsd) 37 37]);
% for nid=1:numel(nrsd)
%     nr = nrsd(nid);
%     for idx = linspace(0,36,37)
%         sp = h5read(sprintf('results/aws-tess-array/dipole/z/full/element%i.hdf5',idx),'/smn');
%         sps = sp.r + 1j.*sp.i;
%         smnf(idx+1,:,:) = sps;
%         sps = h5read(sprintf('results/aws-tess-array/dipole/z/nr%i/element%i.hdf5',nr,idx),'/smn');
%         sps = sps.r + 1j.*sps.i;
%         neighb_id = h5read(sprintf('results/aws-tess-array/dipole/z/nr%i/element%i.hdf5',nr,idx),'/nids');
%         smn(idx+1,:,neighb_id+1) = sps;
%     end
%     % sum errors across frequency
%     errors = abs(smnf-smn);
%     norm_smnf = abs(smnf);
%     
%     nfd(end+1) = norm(abs(smnf(:)-smn(:))./abs(max(smnf(:),[],2)),'fro')/(idx+1);
%     %nfd_error(nid,:,:) = abs(smnf(:)-smn)./abs(max(smnf,[],2));
% end
% 
% %% Dipole Y
% smn = zeros([37,37]);
% smnf = zeros([37,37]);
% nrsdy = [4,6,8,10,12,14,16,18,20,24,28,32,36];
% nfdy = [];
% nfdy_error = zeros([numel(nrsdy) 37 37]);
% for nid=1:numel(nrsdy)
%     nr = nrsdy(nid);
%     for idx = linspace(0,36,37)
%         sp = h5read(sprintf('results/aws-tess-array/dipole/y/full/element%i.hdf5',idx),'/smn');
%         sps = sp.r(151,:) + 1j.*sp.i(151,:);
%         smnf(idx+1,:) = sps;
%         sps = h5read(sprintf('results/aws-tess-array/dipole/y/nr%i/element%i.hdf5',nr,idx),'/smn');
%         sps = sps.r(151,:) + 1j.*sps.i(151,:);
%         neighb_id = h5read(sprintf('results/aws-tess-array/dipole/y/nr%i/element%i.hdf5',nr,idx),'/nids');
%         smn(idx+1,neighb_id+1) = sps;
%     end
%     nfdy(end+1) = norm(abs(smnf-smn)./abs(max(smnf,[],2)),'fro')/(idx+1);
%     nfdy_error(nid,:,:) = abs(smnf-smn)./abs(max(smnf,[],2));
% end
% 
% %% 67 z dipole
% smn = zeros([67,67]);
% smnf = zeros([67,67]);
% nrsd67 = [4,6,8,10,12,14,16,18,20,24,28,32,36,40,44,48,52,58,62,66];
% nfd67 = [];
% nfd67_error = zeros([numel(nrsd67) 67 67]);
% for nid=1:numel(nrsd67)
%     nr = nrsd67(nid);
%     for idx = linspace(0,66,67)
%         sp = h5read(sprintf('results/aws-tess-array/dipole/z/NP5/full/element%i.hdf5',idx),'/smn');
%         sps = sp.r(151,:) + 1j.*sp.i(151,:);
%         smnf(idx+1,:) = sps;
%         sps = h5read(sprintf('results/aws-tess-array/dipole/z/NP5/nr%i/element%i.hdf5',nr,idx),'/smn');
%         sps = sps.r(151,:) + 1j.*sps.i(151,:);
%         neighb_id = h5read(sprintf('results/aws-tess-array/dipole/z/NP5/nr%i/element%i.hdf5',nr,idx),'/nids');
%         smn(idx+1,neighb_id+1) = sps;
%     end
%     nfd67(end+1) = norm(abs(smnf-smn)./abs(max(smnf,[],2)),'fro')/(idx+1);
%     %nfd67_error(nid,:,:) = abs(smnf-smn)./abs(max(smnf,[],2));
% end
% 
% %% 67 y dipole
% smn = zeros([67,67]);
% smnf = zeros([67,67]);
% nrsdy67 = [4,6,8,10,12,14,16,20,24,28,32,36,40,44,48,52,58,62,66];
% nfdy67 = [];
% nfdy67_error = zeros([numel(nrsdy67) 67 67]);
% for nid=1:numel(nrsdy67)
%     nr = nrsdy67(nid);
%     for idx = linspace(0,66,67)
%         sp = h5read(sprintf('results/aws-tess-array/dipole/y/NP5/full/element%i.hdf5',idx),'/smn');
%         sps = sp.r(151,:) + 1j.*sp.i(151,:);
%         smnf(idx+1,:) = sps;
%         sps = h5read(sprintf('results/aws-tess-array/dipole/y/NP5/nr%i/element%i.hdf5',nr,idx),'/smn');
%         sps = sps.r(151,:) + 1j.*sps.i(151,:);
%         neighb_id = h5read(sprintf('results/aws-tess-array/dipole/y/NP5/nr%i/element%i.hdf5',nr,idx),'/nids');
%         smn(idx+1,neighb_id+1) = sps;
%     end
%     nfdy67(end+1) = norm(abs(smnf-smn)./abs(max(smnf,[],2)),'fro')/(idx+1);
%     nfdy67_error(nid,:,:) = abs(smnf-smn)./abs(max(smnf,[],2));
% end
% 
% %%  37 spiral
% smn = zeros([37,37]);
% smnf = zeros([37,37]);
% nrsd = [4,6,8,10,12,14,16,18,20,24,28,32,36];
% nfds = [];
% nfds_error = zeros([numel(nrsd) 37 37]);
% for nid=1:numel(nrsd)
%     nr = nrsd(nid);
%     for idx = linspace(0,36,37)
%         sp = h5read(sprintf('results/aws-tess-array/spiral/full/element%i.hdf5',idx),'/smn');
%         sps = sp.r(151,:) + 1j.*sp.i(151,:);
%         smnf(idx+1,:) = sps;
%         sps = h5read(sprintf('results/aws-tess-array/spiral/nr%i/element%i.hdf5',nr,idx),'/smn');
%         sps = sps.r(151,:) + 1j.*sps.i(151,:);
%         neighb_id = h5read(sprintf('results/aws-tess-array/spiral/nr%i/element%i.hdf5',nr,idx),'/nids');
%         smn(idx+1,neighb_id+1) = sps;
%     end
%     nfds(end+1) = norm(abs(smnf-smn)./abs(max(smnf,[],2)),'fro')/(idx+1);
%     nfds_error(nid,:,:) = abs(smnf-smn)./abs(max(smnf,[],2));
% end
% 
% %plot(nrsp,nfp/2,'bs-'); hold on;
% plot(nrsd,nfd.*100,'ko-'); hold on;
% plot(nrsdy,nfdy.*100,'bx-.');
% plot(nrsd67,nfd67.*100,'r^:');
% plot(nrsdy67,nfdy67.*100,'c*-.');
% %plot(nrsd,nfds.*100,'m->');
% grid minor;
% xlim([2 68]);
% %ylim([0 100]);
% xlabel('$N_R$ Nearest Neighbors','Interpreter','latex');
% ylabel('$100\times\left|\frac{|S_f - S_t|}{|S_f|}\right|_F$','Interpreter','latex');
% legend('Z Dipole 37','Y Dipole 37', 'Z Dipole 67', 'Y Dipole 67');%,'Spiral');