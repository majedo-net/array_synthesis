function [CSX,mesh] = make_spiral(CSX,mesh,center,radii,alpha,h,hs, ...
    Nidx,excite)

%% Calculate top spiral
r0 = radii(1);
rmax = radii(2);
N = log(rmax/r0)/(2*pi*alpha);
a0 = pi/4;
a1 = linspace(a0,N*2*pi,100*N);
a2 = a1 + pi/2;
a3 = a2 + pi/2;
a4 = a3 + pi/2;
r = r0 .* exp(alpha.*a1);
tp = zeros([2 2*numel(r)+52]);
bp = tp;
tp(1,:) = [0,r,r(end).*ones([1 50]),flip(r),0];
tp(2,:) = [0,a1,linspace(a1(end),a2(end),50),flip(a2),0];
bp(1,:) = [0,r,r(end).*ones([1 50]),flip(r),0];
bp(2,:) = [0,a3,linspace(a3(end),a4(end),50),flip(a4),0];

% move coordinates to center if necessary
tp(1,:) = center(1) + tp(1,:);
tp(2,:) = center(2) + tp(2,:);

%create fixed lines for the simulation box and port
mesh.x = [mesh.x center(1)-r0 center(1)-1 center(1)+1 center(1)+r0];
mesh.y = [mesh.y center(2)-r0 center(2)-1 center(2)+1 center(2)+r0];
tp_cart(1,:) = tp(1,:) .* cos(tp(2,:));
tp_cart(2,:) = tp(1,:) .* sin(tp(2,:));
CSX = AddMetal(CSX, sprintf('top_spiral%d',Nidx));
CSX = AddLinPoly( CSX, sprintf('top_spiral%d',Nidx), ...
    10,'z',h+hs,tp_cart,0,'CoordSystem',0);
bp_cart(1,:) = bp(1,:) .* cos(bp(2,:));
bp_cart(2,:) = bp(1,:) .* sin(bp(2,:));
CSX = AddMetal(CSX, sprintf('bot_spiral%d',Nidx));
CSX = AddLinPoly( CSX, sprintf('bot_spiral%d',Nidx), ...
    10,'z',h-hs,bp_cart,0,'CoordSystem',0);

CSX = AddMaterial(CSX,sprintf('substrate%d',Nidx));
CSX = SetMaterialProperty(CSX,sprintf('substrate%d',Nidx), ...
    'Epsilon',4.2);
start = [center(1)-rmax center(2)-rmax h-hs];
stop =  [center(1)+rmax center(2)+rmax h+hs];
CSX = AddBox(CSX, sprintf('substrate%d',Nidx), ...
    0, start, stop,'CoordSystem',0);


%% apply the excitation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

start=[center(1)-1 center(2)-1 h-hs];
stop =[center(1)+1 center(2)+1 h+hs];
[CSX] = AddLumpedPort(CSX,5,Nidx,50,start,stop,[0 0 1],excite);
