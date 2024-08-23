clear; close all;
%% Read the full array simulated coupling matrix
% smnf = fileread('results\fullarray_smn.txt');
% Sf = textscan(smnf,'(%f)',37*37);
% smnf = reshape(Sf{1},[37,37]);


%% Patch antenna 37 elements
smn = zeros([37,37]);
smnf = zeros([37 37]);
nrsp = [4,6,8,10];
nfp_error = zeros([numel(nrsp) 37 37]);
nfp = [];
for nid=1:numel(nrsp)
    nr = nrsp(nid);
    for idx=linspace(0,36,37)
        sp = fileread(sprintf('results/aws-tess-array/patch/full/fullarray_sm%i.txt',idx));
        S=textscan(sp,'(%f)',37);
        sps = S{1};
        smnf(idx+1,:) = sps;
        sp = fileread(sprintf('results/aws-tess-array/patch/nr%i/sn%i.txt',nr,idx));
        S = textscan(sp,'(%f)',2*nr+2);
        neighb_id = int32(S{1}(1:nr+1));
        sps = S{1}(nr+2:2*nr+2);
        smn(idx+1,neighb_id+1) = sps;
    end
    % percent error
    nfp(end+1) = norm(abs(smnf-smn)./abs(max(smnf,[],'all')),'fro')/(idx+1);
    nfp_error(nid,:,:) = abs(smnf-smn)./abs(max(smnf,[],'all'));
end

%% Dipole Z
smn = zeros([37,37]);
smnf = zeros([37,37]);
nrsd = [4,6,8,10,12,14,16,18,20,24,28,32,36];
nfd = [];
nfd_error = zeros([numel(nrsd) 37 37]);
for nid=1:numel(nrsd)
    nr = nrsd(nid);
    for idx = linspace(0,36,37)
        sp = fileread(sprintf('results/aws-tess-array/dipole/z/full/fullarray_sm%i.txt',idx));
        S=textscan(sp,'(%f)',37);
        sps = S{1};
        smnf(idx+1,:) = sps;
        sp = fileread(sprintf('results/aws-tess-array/dipole/z/tess/nr%i/sm%i.txt',nr,idx));
        S=textscan(sp,'(%f)',2*nr+2);
        neighb_id = int32(S{1}(1:nr+1));
        sps = S{1}(nr+2:2*nr+2);
        smn(idx+1,neighb_id+1) = sps;
    end
    nfd(end+1) = norm(abs(smnf-smn)./abs(max(smnf,[],'all')),'fro')/(idx+1);
    nfd_error(nid,:,:) = abs(smnf-smn)./abs(max(smnf,[],'all'));
end

%% Dipole Y
smn = zeros([37,37]);
smnf = zeros([37,37]);
nrsdy = [4,6,8,10,12,14,16,18,20,24,28,32,36];
nfdy = [];
nfdy_error = zeros([numel(nrsdy) 37 37]);
for nid=1:numel(nrsdy)
    nr = nrsdy(nid);
    for idx = linspace(0,36,37)
        sp = fileread(sprintf('results/aws-tess-array/dipole/y/full/fullarray_sm%i.txt',idx));
        S=textscan(sp,'(%f)',37);
        sps = S{1};
        smnf(idx+1,:) = sps;
        sp = fileread(sprintf('results/aws-tess-array/dipole/y/tess/nr%i/sm%i.txt',nr,idx));
        S=textscan(sp,'(%f)',2*nr+2);
        neighb_id = int32(S{1}(1:nr+1));
        sps = S{1}(nr+2:2*nr+2);
        smn(idx+1,neighb_id+1) = sps;
    end
    nfdy(end+1) = norm(abs(smnf-smn)./abs(max(smnf,[],'all')),'fro')/(idx+1);
    nfdy_error(nid,:,:) = abs(smnf-smn)./abs(max(smnf,[],'all'));
end

%% 67 z dipole
smn = zeros([67,67]);
smnf = zeros([67,67]);
nrsd67 = [4,6,8,10,12,14,16,18,20,24,28,32,36,40,44,48,52,58,62,66];
nfd67 = [];
nfd67_error = zeros([numel(nrsd67) 67 67]);
for nid=1:numel(nrsd67)
    nr = nrsd67(nid);
    for idx = linspace(0,66,67)
        sp = fileread(sprintf('results/aws-tess-array/dipole/z/NP5/full/fullarray_sm%i.txt',idx));
        S=textscan(sp,'(%f)',67);
        sps = S{1};
        smnf(idx+1,:) = sps;
        if nr > 18
            sps = h5read(sprintf('results/aws-tess-array/dipole/z/NP5/tess/nr%i/element%i.hdf5',nr,idx),'/smn');
            sps = sps.r(151,:) + 1j.*sps.i(151,:);
            neighb_id = h5read(sprintf('results/aws-tess-array/dipole/z/NP5/tess/nr%i/element%i.hdf5',nr,idx),'/nids');
            smn(idx+1,neighb_id+1) = sps;
        else
            sp = fileread(sprintf('results/aws-tess-array/dipole/z/NP5/tess/nr%i/sm%i.txt',nr,idx));
            S=textscan(sp,'(%f)',2*nr+2);
            neighb_id = int32(S{1}(1:nr+1));
            sps = S{1}(nr+2:2*nr+2);
            smn(idx+1,neighb_id+1) = sps;
        end
    end
    nfd67(end+1) = norm(abs(smnf-smn)./abs(max(smnf,[],'all')),'fro')/(idx+1);
    nfd67_error(nid,:,:) = abs(smnf-smn)./abs(max(smnf,[],'all'));
end

%% 67 y dipole
smn = zeros([67,67]);
smnf = zeros([67,67]);
nrsdy67 = [4,6,8,10,12,14,16,20,24,28,32,36,40,44,48,52,58,62,66];
nfdy67 = [];
nfdy67_error = zeros([numel(nrsdy67) 67 67]);
for nid=1:numel(nrsdy67)
    nr = nrsdy67(nid);
    for idx = linspace(0,66,67)
        sp = fileread(sprintf('results/aws-tess-array/dipole/y/tess/NP5/full/fullarray_sm%i.txt',idx));
        S=textscan(sp,'(%f)',67);
        sps = S{1};
        smnf(idx+1,:) = sps;
        sps = h5read(sprintf('results/aws-tess-array/dipole/y/tess/NP5/nr%i/element%i.hdf5',nr,idx),'/smn');
        sps = sps.r(151,:) + 1j.*sps.i(151,:);
        neighb_id = h5read(sprintf('results/aws-tess-array/dipole/y/tess/NP5/nr%i/element%i.hdf5',nr,idx),'/nids');
        smn(idx+1,neighb_id+1) = sps;
    end
    nfdy67(end+1) = norm(abs(smnf-smn)./abs(max(smnf,[],'all')),'fro')/(idx+1);
    nfdy67_error(nid,:,:) = abs(smnf-smn)./abs(max(smnf,[],'all'));
end

%%  37 spiral
smn = zeros([37,37]);
smnf = zeros([37,37]);
nrsd = [4,6,8,10,12,14,16,18,20,24,28,32,36];
nfds = [];
nfds_error = zeros([numel(nrsd) 37 37]);
for nid=1:numel(nrsd)
    nr = nrsd(nid);
    for idx = linspace(0,36,37)
        sp = h5read(sprintf('results/aws-tess-array/spiral/full/element%i.hdf5',idx),'/smn');
        sps = sp.r(151,:) + 1j.*sp.i(151,:);
        smnf(idx+1,:) = sps;
        sps = h5read(sprintf('results/aws-tess-array/spiral/nr%i/element%i.hdf5',nr,idx),'/smn');
        sps = sps.r(151,:) + 1j.*sps.i(151,:);
        neighb_id = h5read(sprintf('results/aws-tess-array/spiral/nr%i/element%i.hdf5',nr,idx),'/nids');
        smn(idx+1,neighb_id+1) = sps;
    end
    nfds(end+1) = norm(abs(smnf-smn)./abs(max(smnf,[],'all')),'fro')/(idx+1);
    nfds_error(nid,:,:) = abs(smnf-smn)./abs(max(smnf,[],'all'));
end

%plot(nrsp,nfp/2,'bs-'); hold on;
plot(nrsd,nfd.*100,'ko-'); hold on;
plot(nrsdy,nfdy.*100,'bx-.');
plot(nrsd67,nfd67.*100,'r^:');
plot(nrsdy67,nfdy67.*100,'c*-.');
%plot(nrsd,nfds.*100,'m->');
grid minor;
xlim([2 68]);
%ylim([0 100]);
xlabel('$N_R$ Nearest Neighbors','Interpreter','latex');
ylabel('$100\times\left|\frac{|S_f - S_t|}{|S_f|}\right|_F$','Interpreter','latex');
legend('Z Dipole 37','Y Dipole 37', 'Z Dipole 67', 'Y Dipole 67','Spiral');