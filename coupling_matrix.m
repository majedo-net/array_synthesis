%% Read the full array simulated coupling matrix
smnf = fileread('results\fullarray_smn.txt');
Sf = textscan(smnf,'(%f)',37*37);
smnf = reshape(Sf{1},[37,37]);
%% Construct the tessellated array coupling matrix
smn = zeros([37,37]);
for idx=linspace(0,36,37)
    sp = fileread(sprintf('results/aces_tess_array_revised/sn%i.txt',idx));
    S = textscan(sp,'(%f)',14);
    neighb_id = int32(S{1}(1:7));
    sps = S{1}(8:14);
    smn(idx+1,neighb_id+1) = sps;
end
%smnf = smnf./max(smnf,[],'all');
%smn = smn./max(smn,[],'all');
% frobenius norm of difference
nf = norm(smnf-smn,'fro');
disp(nf);