clear; close all;

L = 0.6;
z = linspace(-L/2,L/2,201);
f = 25e9;
lambda = 3e8/f;
k = 2.*pi./lambda;

th0 = deg2rad(60);
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
    
    theta = linspace(0,pi,501);
    
    AF = Iz*exp(-1j.*(k.*z.*sin(theta'-th0')))';
    AF = 20.*log10(abs(AF./max(AF)));
    figure(2);
    plot(rad2deg(theta),AF); hold on;
end

figure(1);
xlabel('Radial Position (m)');
ylabel('Amplitude Coefficient');
grid minor;
legend('SLL = 13.26dB','SLL=20dB','SLL=25dB','SLL=40dB','Location','south');

figure(2);
xlabel('Theta');
ylabel('Normalized Array Factor (dB)');
grid minor;
legend('SLL = 13.26dB','SLL=20dB','SLL=25dB','SLL=40dB','Location','south');
