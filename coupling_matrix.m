clear; close all;
%% Read the full array simulated coupling matrix
smnf = fileread('results\fullarray_smn.txt');
Sf = textscan(smnf,'(%f)',37*37);
smnf = reshape(Sf{1},[37,37]);
%% Construct the tessellated array coupling matrix
smn = zeros([37,37]);
nrsp = [4,6,8,10];
nfp = [];
for nid=1:numel(nrsp)
    nr = nrsp(nid);
    for idx=linspace(0,36,37)
        sp = fileread(sprintf('results/aws-tess-array/patch/nr%i/sn%i.txt',nr,idx));
        S = textscan(sp,'(%f)',2*nr+2);
        neighb_id = int32(S{1}(1:nr+1));
        sps = S{1}(nr+2:2*nr+2);
        smn(idx+1,neighb_id+1) = sps;
    end
    % frobenius norm of difference
    nfp(end+1) = sqrt(norm(smnf-smn,'fro'));
end

%% Dipole Z
smn = zeros([37,37]);
smnf = zeros([37,37]);
nrsd = [4,6,8,10,12,14,16,18];
nfd = [];
for nid=1:numel(nrsd)
    nr = nrsd(nid);
    for idx = linspace(0,36,37)
        sp = fileread(sprintf('results/aws-tess-array/dipole/full/fullarray_sm%i.txt',idx));
        S=textscan(sp,'(%f)',37);
        sps = S{1};
        smnf(idx+1,:) = sps;
        sp = fileread(sprintf('results/aws-tess-array/dipole/tess/nr%i/sm%i.txt',nr,idx));
        S=textscan(sp,'(%f)',2*nr+2);
        neighb_id = int32(S{1}(1:nr+1));
        sps = S{1}(nr+2:2*nr+2);
        smn(idx+1,neighb_id+1) = sps;
    end
    nfd(end+1) = sqrt(norm(smnf-smn,'fro'));
end

%% 67 dipole
smn = zeros([67,67]);
smnf = zeros([67,67]);
nrsd = [4,6,8,10,12,14,16,18];
nfd67 = [];
for nid=1:numel(nrsd)
    nr = nrsd(nid);
    for idx = linspace(0,66,67)
        sp = fileread(sprintf('results/aws-tess-array/dipole/NP5/full/fullarray_sm%i.txt',idx));
        S=textscan(sp,'(%f)',67);
        sps = S{1};
        smnf(idx+1,:) = sps;
        sp = fileread(sprintf('results/aws-tess-array/dipole/NP5/tess/nr%i/sm%i.txt',nr,idx));
        S=textscan(sp,'(%f)',2*nr+2);
        neighb_id = int32(S{1}(1:nr+1));
        sps = S{1}(nr+2:2*nr+2);
        smn(idx+1,neighb_id+1) = sps;
    end
    nfd67(end+1) = sqrt(norm(smnf-smn,'fro'));
end

%plot(nrsp,nfp,'bs-'); hold on;
plot(nrsd,nfd,'ko-'); hold on;
plot(nrsd,nfd67,'r^-');
grid minor;
xlim([2 20]);
xlabel('$N_R$ Nearest Neighbors','Interpreter','latex');
ylabel('$|S_f - S_t|_F$','Interpreter','latex');
legend('Dipole 37 Elements','Dipole 67 Elements');