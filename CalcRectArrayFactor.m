function [xs,ys,time_delays]=CalcRectArrayFactor(freq,d,theta_0,phi_0,options)
%%
% @file CalcArrayFactor.m
%
% @brief Synthesize an array with symmetric, but not necessarily uniform
% spacing, and calculate the required phase shifts/time delays for a given
% steering angle. Finally, plot the array factor and optionally total
% pattern if provided an element pattern.
% 
% input:
%   freq: Frequency to evaluate array factor
%   
%   d: vector of element spacings for one quadrant of the planar array,
%   starting from center element
%   
%   theta_0: Desired steering angle in elevation, broadside is zero
%   
%   phi_0: Desired steering angle in azimuth, -y is zero
%   
%   element_pattern_file: path to antenna element pattern file with
%   electric fiel/gain as function of theta and phi. Interposer function
%   required to import simulated pattern
%
% @copyright Copyright (c) 2022 Augustus Aerospace Company, all rights reserved.
%
arguments
    freq
    d
    theta_0
    phi_0
    options.element_pattern_file =0
    options.PlotsOn (1,1) {mustBeNumeric} = 1
    options.Quantize (1,1) {mustBeNumeric} = 0
end


set(0,'DefaultAxesFontSize',18);
set(0,'DefaultLineLinewidth',2);
set(0,'DefaultTextFontSize',18);

if nargin < 4
    disp('Missing arguments');
    return;
end
if freq<1
    disp('Invalid Frequency');
    return;
end
if numel(d)<2
    disp('Invalid spacing array, must be at least 2 elements');
    return;
end
if (theta_0 < 0) || (theta_0>180)
    disp('Theta_0 must be between 0 and 180');
    return;
end
if (phi_0 < 0) || (phi_0 > 360)
    disp('Phi_0 must be between 0 and 360');
    return;
end



% Load Element Pattern from File
if options.element_pattern_file
    element_pattern=load_element(options.element_pattern_file);
    %normalize element pattern
    element_pattern = element_pattern - max(element_pattern,[],'all');
    element_pattern(element_pattern< -20) = -20;
end


N = length(d)-1;
theta_0 = deg2rad(theta_0);
phi_0 = deg2rad(phi_0);
lambda = physconst('Lightspeed')./freq;
k = 2.*pi./lambda;
k = reshape(k,1,1,[]);
theta = linspace(0,pi,181);
theta_yz = cat(2,flip(-theta,2),theta);
theta_xz = theta_yz;
phi = transpose(linspace(0,2*pi,361));
[el,az] = meshgrid(theta,phi);

% Generate array pattern
xs = zeros((2*N+1)^2,1);
ys = zeros((2*N+1)^2,1);
id = 1;
for nd = 0:(size(d)-1)
    yd = nd;
    xd = yd;
    while yd >= 0
        xs(id) = xd.*(d(nd+1));
        ys(id) = yd.*(d(nd+1));
        id=id+1;
        if (xd~=0)
            xs(id)= -xd*(d(nd+1));
            ys(id)= yd*(d(nd+1));
            id = id +1;
        end
        if (yd~=0)
            xs(id)= xd*(d(nd+1));
            ys(id)= -yd*(d(nd+1));
            id = id +1;
        end
        if (xd~=0)&(yd~=0)
            xs(id)= -xd*(d(nd+1));
            ys(id)= -yd*(d(nd+1));
            id = id +1;
        end
        yd = yd-1;
    end
    if xd>0
        xd = xd-1;
        while xd>=0
            yd = nd;
            xs(id)=xd*(d(nd+1));
            ys(id)=yd*(d(nd+1));
            id = id+1;
            if (yd~=0)&(xd~=yd)
                xs(id)=xd*(d(nd+1));
                ys(id)=-yd*(d(nd+1));
                id = id+1;
            end
            if (xd~=0)&(xd~=yd)
                xs(id)=-xd*(d(nd+1));
                ys(id)=yd*(d(nd+1));
                id = id+1;
                xs(id)=-xd*(d(nd+1));
                ys(id)=-yd*(d(nd+1));
                id = id+1;
            end
            xd = xd-1;
        end
    end
end

%xs = sort(xs);
%ys = sort(ys);

% Plot element arrangement
figure();
plot(xs,ys,'ko');
grid minor;
title('Array Element Locations');

% Calculate array factor
ang_pos = atan(ys./xs);


path_len_delta = (xs.*sin(theta_0).*cos(phi_0)+ys.*sin(theta_0).*sin(phi_0));
time_delays = (path_len_delta-min(path_len_delta))./ 3e8;
if options.Quantize
    partitions = 25e-12.*linspace(1,127,127);
    index = quantiz(time_delays,partitions);
    codebook = 25e-12.*linspace(0,127,128);
    time_delays = codebook(index+1).';
    path_len_delta = time_delays .* 3e8;
end
steering_vectors = exp(-1j.*k.*path_len_delta);
incident_waves = zeros(361,181,size(xs,1));
AF = 0;

for n = 1:size(xs)
    if isnan(ang_pos(n))
        ang_pos(n) = 0;
    end
    xs3 = reshape(xs(n),1,1,[]);
    ys3 = reshape(ys(n),1,1,[]);
    incident_waves(:,:,n) = exp(1j.*k.*(xs3.*(sin(theta).*cos(phi))+ys3.*(sin(theta).*sin(phi))));
    AF = AF + (steering_vectors(n).*incident_waves(:,:,n));   
end



% Fix null values at theta/phi = 0
AF(isnan(AF)) = 0;
% Normalize  and convert to db array factor
AF = 10.*log10(abs(AF)./max(abs(AF),[],'all'));
AF(AF< -20) = -20;

%% Plots
if options.PlotsOn

    if options.element_pattern_file
        % Plot element pattern
        figure();
        ax1=subplot(1,3,1,polaraxes);
        el_xy = transpose(element_pattern(:,90));
        polarplot(ax1,phi,el_xy);
        rlim([-20 0]);
        title(ax1,'Element Pattern XY-Plane');
        ax1.ThetaZeroLocation = 'right';
        subtitle('\leftarrow \phi \rightarrow');
        
        ax2=subplot(1,3,2,polaraxes);
        el_yz = cat(2,flip(element_pattern(270,:),2),element_pattern(90,:));
        polarplot(ax2,theta_yz,el_yz);
        rlim([-20 0]);
        title(ax2,'Element Pattern YZ-Plane');
        ax2.ThetaTickLabels = {'0';'30';'60';'90';'120';'150';'180';'150';'120';'90';'60';'30'};
        ax2.ThetaZeroLocation = 'top';
        subtitle('\leftarrow \theta \rightarrow');
        
        ax3=subplot(1,3,3,polaraxes);
        el_xz = cat(2,flip(element_pattern(180,:),2),element_pattern(1,:));
        polarplot(ax3,theta_xz,el_xz);
        rlim([-20 0]);
        title(ax3,'Element Pattern XZ-Plane');
        ax3.ThetaTickLabels = {'0';'30';'60';'90';'120';'150';'180';'150';'120';'90';'60';'30'};
        ax3.ThetaZeroLocation = 'top';
        subtitle('\leftarrow \theta \rightarrow');  
        
        % Element 3d plot
        figure;
        [X,Y,Z] = sph2cart(az,(el-pi/2),(10.^(element_pattern./10)));
        s=surf(X,Y,Z,'CData',(10.^(element_pattern./10)),'FaceColor','interp');
        zlim([-1 1]);
        s.EdgeColor = 'None';
        title('Element Pattern');
        xlabel('X');
        ylabel('Y');
        zlabel('Z');
        colorbar;
    end
    
    % Plot Array Factor------------------------------
    figure();
    
    ax1=subplot(1,3,1,polaraxes);
    sgtitle(['Freq = ',num2str(freq,3),' Hz','\newline', ...
        '\theta_0 =',num2str(rad2deg(theta_0)),'\phi_0 =',num2str(rad2deg(phi_0)),'\newline', ...
        'd space =',strjoin(string(d))]);
    
    
    
    AF_xy = transpose(AF(:,90));
    polarplot(ax1,phi,AF_xy);
    rlim([-20 0]);
    title(ax1,'Array Factor XY-Plane');
    ax1.ThetaZeroLocation = 'right';
    subtitle('\leftarrow \phi \rightarrow');
    
    ax2=subplot(1,3,2,polaraxes);
    AF_yz = cat(2,flip(AF(270,:),2),AF(90,:));
    polarplot(ax2,theta_yz,AF_yz);
    rlim([-20 0]);
    title(ax2,'Array Factor YZ-Plane');
    ax2.ThetaTickLabels = {'0';'30';'60';'90';'120';'150';'180';'150';'120';'90';'60';'30'};
    ax2.ThetaZeroLocation = 'top';
    subtitle('\leftarrow \theta \rightarrow');
    
    ax3=subplot(1,3,3,polaraxes);
    AF_xz = cat(2,flip(AF(180,:),2),AF(1,:));
    polarplot(ax3,theta_xz,AF_xz);
    rlim([-20 0]);
    title(ax3,'Array Factor XZ-Plane');
    ax3.ThetaTickLabels = {'0';'30';'60';'90';'120';'150';'180';'150';'120';'90';'60';'30'};
    ax3.ThetaZeroLocation = 'top';
    subtitle('\leftarrow \theta \rightarrow');
    
    % AF 3d plot
    figure();
    [az,el] = meshgrid(phi,theta);
    [X,Y,Z] = sph2cart(az,(el-pi/2),transpose(10.^(AF(:,:)./10)));
    s=surf(X,Y,Z,'CData',10.^(transpose(AF(:,:))./10),'FaceColor','interp');
    shading interp;
    s.EdgeColor = 'None';
    title('Array Factor');
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    colorbar;
    
    if options.element_pattern_file
        % Compute the total pattern
        total = AF + element_pattern;
        total = total - max(total,[],'all');
        total(total<-20) = -20;
        
        %------------------------------
        % Plot the total pattern
        figure();
        ax1=subplot(1,3,1,polaraxes);
        sgtitle(['Freq = ',num2str(freq,3),' Hz']);
        total_xy = transpose(total(:,90));
        polarplot(ax1,phi,total_xy);
        rlim([-20 0]);
        title(ax1,'Total XY-Plane');
        ax1.ThetaZeroLocation = 'right';
        subtitle('\leftarrow \phi \rightarrow');
        
        ax2=subplot(1,3,2,polaraxes);
        total_yz = cat(2,flip(total(270,:),2),total(90,:));
        polarplot(ax2,theta_yz,total_yz);
        rlim([-20 0]);
        title(ax2,'Total YZ-Plane');
        ax2.ThetaTickLabels = {'0';'30';'60';'90';'120';'150';'180';'150';'120';'90';'60';'30'};
        ax2.ThetaZeroLocation = 'top';
        subtitle('\leftarrow \theta \rightarrow');
        
        ax3=subplot(1,3,3,polaraxes);
        total_xz = cat(2,flip(total(180,:),2),total(1,:));
        polarplot(ax3,theta_xz,total_xz);
        rlim([-20 0]);
        title(ax3,'Total XZ-Plane');
        ax3.ThetaTickLabels = {'0';'30';'60';'90';'120';'150';'180';'150';'120';'90';'60';'30'};
        ax3.ThetaZeroLocation = 'top';
        subtitle('\leftarrow \theta \rightarrow');
        
        % Total 3d plot
        figure;
        [az,el] = meshgrid(phi,theta);
        [X,Y,Z] = sph2cart(az,(el-pi/2),transpose(10.^(total./10)));
        s=surf(X,Y,Z,'CData',10.^(transpose(total)./10),'FaceColor','interp');
        shading interp
        s.EdgeColor = 'None';
        xlabel('X');
        ylabel('Y');
        zlabel('Z');
        title(['Total Pattern, ','Freq = ',num2str(freq,3),' Hz']);
        colorbar;
    end % if element pattern file
end % if plots on

end