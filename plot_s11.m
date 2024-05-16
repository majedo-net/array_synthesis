clear; close all;

max_eid = 6;
figure;
hold on;
for id=0:max_eid
    dat = fileread(sprintf('results/aws-tess-array/dipole/full/full_s11_%d.txt',id));
    dat = textscan(dat,'(%f)',602);
    dat = dat{1};
    freq = dat(1:301)./1e9;
    s11 = 10.*log10(abs(dat(302:602)));
    plot(freq,s11,'-.','LineWidth',1);
    xlim([5.5 6.5]);
    ylim([-35 0]);
    if any(s11>0)
        fprintf('Bad eid: %d \n',id);
    end
end
% dat=readmatrix('results/aces_abstract_run/s11_single_patch.txt');
% plot(dat(1,:)./1e9,dat(2,:),'k-','LineWidth',3);
% qw{1} = plot(nan,'k-','LineWidth',3);
% qw{2} = plot(nan,'k-.','LineWidth',1);
% legend([qw{:}],{'Isolated Element','Elements With Mutual Coupling'},'location','best');
grid minor;
xlabel('Frequency (GHz)');
ylabel('S_{11} dB');