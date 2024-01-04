clear; close all;

max_eid = 36;
figure;
hold on;
for id=0:max_eid
    dat = readmatrix(sprintf('results/s11_%d.txt',id));
    freq = dat(1,:);
    s11 = dat(2,:);
    plot(freq,s11);
    xlim([4e9 9e9]);
    if any(s11>0)
        fprintf('Bad eid: %d \n',id);
    end
end
legend;