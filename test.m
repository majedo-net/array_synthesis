fmax = 6e9;

d = RaisedPowerSeries(fmax,0.8,3);

[xs,ys,td]=CalcArrayFactor(0.9e9,d,30,30,PlotsOn=0);
[X,Y] = meshgrid(xs,ys);
X = sort(X);
Y = sort(Y);
Z = 
contourf(X,Y,Z);
colorbar;