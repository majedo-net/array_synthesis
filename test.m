clear;
close all;
fmax = 7.5e9;
fplot = 2e9;
d = RaisedPowerSeries(fmax,0.85,4);

[xs,ys,td]=CalcRectArrayFactor(fplot,d,60,45,PlotsOn=1,Quantize=0);
%td = td ./ max(td,[],'all'); %normalize time delays
% td_bits = dec2bin(td./(25e-12)); % compute 7 bit delay with 25ps lsb
% tdmap = cat(2,xs,ys,td_bits);
% [X,Y] = meshgrid(unique(xs),unique(ys));
% tdz = zeros(length(X),length(Y));
% for idx=1:length(tdmap)
%     [i,~] = find(Y==tdmap(idx,2));
%     [~,j] = find(X==tdmap(idx,1));
%     i = i(1);
%     j= j(1);  
%     tdz(i,j)=tdmap(idx,3);
% end
% xs = sort(unique(xs));
% ys = sort(unique(ys));
% fprintf('Max Time Delay: %g ps \n Min Time Delay: %g ps\n',10^12.*max(tdz,[],'all'),10^12.*min(tdz(tdz>0),[],"all"));
% heatmap(xs,ys,tdz,'Colormap',turbo);