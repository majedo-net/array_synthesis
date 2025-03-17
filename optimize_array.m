close all; clear;
N = 50;

theta = optimvar('theta',N,'LowerBound',0,'UpperBound',2.*pi);
thetaprob = optimproblem;
energy = optimexpr(1);

rs = abs(randn(N,1));
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
