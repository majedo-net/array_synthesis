close all; clear;
N = 20;


theta = optimvar('theta',N,'LowerBound',0,'UpperBound',2.*pi);
thetaprob = optimproblem;
energy = optimexpr(1);

rw = linspace(-10,10,501);
w = taylorwin(501);
rs = rejectionSample(w,rw,N)';
ws = interp1(rw,w,rs);

init.theta = rand(N,1).*2.*pi;

x = rs.*cos(theta);
y = rs.*sin(theta);
for ii = 1:(N-1)
    jj = (ii+1):N;
    tempe = (x(ii)-x(jj)).^2 + (y(ii)-y(jj)).^2;
    energy = energy + sum(tempe.^(-1/2));
end
thetaprob.Objective = energy;

x0 = rs.*cos(init.theta);
y0 = rs.*sin(init.theta);

[sol,fval,exitlfag,output] = solve(thetaprob,init);

xf = rs.*cos(sol.theta);
yf = rs.*sin(sol.theta);
dx = xf -x0;
dy = yf -y0;
scatter(x0,y0,'filled'); hold on;
scatter(xf,yf,'filled');
legend('Start','Finish');
quiver(x0,y0,dx,dy,'off');

theta = linspace(0,pi,181)';
phi = linspace(0,2*pi,361);

AF = zeros(numel(theta),numel(phi));
for id=1:N
    AF =AF + ws(id).*exp(1j.*2.*pi.*(xf(id).*sin(theta).*cos(phi)+yf(id).*sin(theta).*sin(phi)));
end

[TH,PH] = ndgrid(theta,phi);
U = sin(TH).*cos(PH);
V = sin(TH).*sin(PH);
figure;
surf(U,V,10.*log10(abs(AF)),'EdgeColor','interp');
colorbar;
view(2);

figure;
for phid=1:10:numel(phi)
    plot(theta,10.*log10(abs(AF(:,phid)))); hold on;
end