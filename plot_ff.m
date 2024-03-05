clear; close all;

max_eid = 25;
theta=linspace(0,2*pi,361);
epax = polaraxes;
hpax = polaraxes;
subplot(1,2,1,epax);
subplot(1,2,2,hpax);



for id=0:max_eid
    dat = abs(readmatrix(sprintf('results/ff_%d.txt',id)));
    if dat(1,1) > dat(end,1)
        Eplane = [dat(:,1)' , flip(dat(2:end,181)')];
        Hplane = [dat(:,91)' , flip(dat(2:end,271)')];
    else
        Eplane = [flip(dat(:,1)') , dat(2:end,181)'];
        Hplane = [flip(dat(:,91)') , dat(2:end,271)'];
    end
    Eplane = Eplane./max(Eplane,[],"all");
    Hplane = Hplane./max(Hplane,[],"all");
    Eplanedb = 20000*log10(Eplane);
    Eplanedb(Eplanedb < -20) = -20; 
    Eplanedb = movmean(Eplanedb,12);
    Hplanedb = 20000*log10(Hplane);
    Hplanedb(Hplanedb< -20) = -20;
    Hplanedb = movmean(Hplanedb,12);
    polaraxes(epax);
    polarplot(theta,Eplanedb,'-.','LineWidth',1);
    hold on;
    polaraxes(hpax);
    polarplot(theta,Hplanedb,'-.','LineWidth',1);
    hold on;
end
dat = readmatrix("results/ff_single_patch.txt");

Eplane = [dat(:,1)' , flip(dat(2:end,181)')];
Eplane = Eplane./max(Eplane,[],"all");
Eplanedb = 20000*log10(Eplane);

Hplane = [dat(:,91)' , flip(dat(2:end,271)')];
Hplane = Hplane./max(Hplane,[],"all");
Hplanedb = 20000*log10(Hplane);

polaraxes(epax);
polarplot(theta,Eplanedb,'k-','LineWidth',3);
hold on;
polaraxes(hpax);
polarplot(theta,Hplanedb,'r-','LineWidth',3);
hold on;

epax.ThetaZeroLocation = 'top';
epax.RLim = [-20 0];
hpax.ThetaZeroLocation = 'top';
hpax.RLim = [-20 0];
subtitle('H-Plane');
hpax.ThetaTickLabel = {'0','30','60','90','120','150','180',...
    '150','120','90','60','30'};
text(pi/6,6.5,'\leftarrow\phi = 90','HorizontalAlignment','center','FontSize',12);
text(-pi/6,6.5,'\phi = 270\rightarrow','HorizontalAlignment','center','FontSize',12);
epax.ThetaTickLabel = {'0','30','60','90','120','150','180',...
    '150','120','90','60','30'};
polaraxes(epax);
text(pi/6,6.5,'\leftarrow\phi = 0','HorizontalAlignment','center','FontSize',12);
text(-pi/6,6.5,'\phi = 180\rightarrow','HorizontalAlignment','center','FontSize',12);
subtitle('E-Plane');