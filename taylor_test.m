clear; close all;

f = 15e9;
lambda = 3e8/f;
k = 2.*pi./lambda;

f2 = 5e9;
lambda2 = 3e8/f2;
k2 = 2.*pi./lambda2;

N = 51;
z = linspace(-lambda*(N-1)/4,lambda*(N-1)/4,N);
L = lambda*N;

th0 = deg2rad(0);
scanw = exp(1j.*k.*z.*cos(th0));

%% Taylor one parameter distribution
%   B           SLL(dB)     Efficiency
%   0           13.26       1
%   0.7386      20          0.933
%   1.0229      25          0.995
%   1.7415      40          0.709
%%
Bs =[0, 0.7386,1.0229,1.7415];

for ii = 1:numel(Bs)
    B = Bs(ii);
    Jarg = 1j.*pi.*B.*sqrt(1-(2.*z./L).^2);
    Iz = besselj(0,Jarg);
    figure(1);
    plot(z,Iz./(max(Iz))); hold on;

    Iznorm = Iz./max(Iz);
    dn = (1./Iznorm).*lambda/2;
    idx_middle = (N-1)/2 +1;
    zn(idx_middle) = 0;
    for idx=idx_middle+1:N
        zn(idx) = zn(idx-1) + dn(idx);
        opp_idx = 2*idx_middle - idx;
        zn(opp_idx) = zn(opp_idx+1) - dn(opp_idx);
    end
    scatter(zn',ii./4,'filled');
    theta = linspace(0,pi/2,501);
    
    AF = Iz*exp(-1j.*(k.*zn.*sin(theta'-th0')))';
    AF = 20.*log10(abs(AF./max(AF)));
    
    AF2 = Iz*exp(-1j.*(k2.*zn.*sin(theta'-th0')))';
    AF2 = 20.*log10(abs(AF2./max(AF2)));

    figure(2);
    plot(rad2deg(theta),AF2); hold on;

    figure(3);
    plot(rad2deg(theta),AF); hold on;
end

figure(1);
xlabel('Radial Position (m)');
ylabel('Amplitude Coefficient');
grid minor;
legend('SLL = 13.26dB','SLL=20dB','SLL=25dB','SLL=40dB');

figure(2);
title(sprintf('Frequency = %f GHz',f2/1e9));
xlabel('Theta');
ylabel('Normalized Array Factor (dB)');
grid minor;
ylim([-50,0]);
legend('SLL = 13.26dB','SLL=20dB','SLL=25dB','SLL=40dB');

figure(3);
title(sprintf('Frequency = %f GHz',f/1e9));
xlabel('Theta');
ylabel('Normalized Array Factor (dB)');
grid minor;
ylim([-50,0]);
legend('SLL = 13.26dB','SLL=20dB','SLL=25dB','SLL=40dB','Location','south');
