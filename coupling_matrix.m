%% Read the full array simulated coupling matrix

%% Construct the tessellated array coupling matrix
smn = zeros([37,37]);
for idx=linspace(0,36,37)
    sp = fileread(sprintf('results/aces_tess_array_revised/sn%i.txt',idx));
    S = textscan(sp,'(%f)',14);
    neighb_id = int32(S{1}(1:7));
    sps = S{1}(8:14);
    smn(idx+1,neighb_id+1) = sps;
end