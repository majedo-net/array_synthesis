clear; close all;
%% Read the full array simulated coupling matrix
smnf = fileread('results\fullarray_smn.txt');
Sf = textscan(smnf,'(%f)',37*37);
smnf = reshape(Sf{1},[37,37]);
%% Construct the tessellated array coupling matrix
smn = zeros([37,37]);
nrs = [6,8,10];
nf = [];
for nid=1:numel(nrs)
    nr = nrs(nid);
    for idx=linspace(0,36,37)
        sp = fileread(sprintf('results/aws-tess-array/nr%i/sn%i.txt',nr,idx));
        S = textscan(sp,'(%f)',2*nr+2);
        neighb_id = int32(S{1}(1:nr+1));
        sps = S{1}(nr+2:2*nr+2);
        smn(idx+1,neighb_id+1) = sps;
    end
    % frobenius norm of difference
    nf(end+1) = norm(smnf-smn,'fro');
end

plot(nrs,nf,'bs-');
grid minor;
xlim([2 14]);
xlabel('k-Nearest Neighbors');
ylabel('$|S_f - S_t|_F$','Interpreter','latex');