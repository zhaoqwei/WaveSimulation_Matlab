clc,clear
% 3D - Acoustic wave simulation using finite difference. with staggered-grid.
% By zhaoqingwei
% Chengdu University of Technology (CDUT), 2021-2025
%%%%%%%%%%%%%%%%%%%%%%%%%
%         P--Vx--P
%       /       /
%     Vy --- Vy
%   /      /
% P--Vx--P
% |      |
% Vz-----vz
% |      |
% P--Vx--P
%%%%%%%%%%%%%%%%%%%%%%%%%

FM=20;
DT=0.001;dt=DT;
T=2;
nx=200;
ny=200;
nz=200;
DH=10;dx=DH;dz=DH;dy=DH;

s=wavelet(FM,DT,T);
vx0=zeros(nz,nx,ny);vx1=vx0;
vz0=vx0;vz1=vx0;
vy0=vx0;vy1=vx0;
p0=vx0;p1=vx0;
v=ones(nz,nx,ny)*2000;
z0=round(nz/2);x0=round(nx/2);y0=round(ny/2);

nn=3;
dxd=FDcoeffDx(nn);ddz0(:,1,1)=dxd';ddz1(:,1,1)=[dxd 0]';
ddx0=permute(ddz0,[2 1 3]);ddx1=permute(ddz1,[2 1 3]);ddy0=permute(ddz0,[3 2 1]);ddy1=permute(ddz1,[3 2 1]);
TX=dt/dx;TZ=dt/dz;TY=dt/dy;

for	t=DT:DT:T
    k=round(t/DT);
    vz1=vz0+TZ*imfilter(p0,ddz0);
    vx1=vx0+TX*imfilter(p0,ddx0);
    vy1=vy0+TY*imfilter(p0,ddy0);
    p1=p0+v.*v.*TZ.*imfilter(vz1,ddz1)+v.*v.*TX.*imfilter(vx1,ddx1)+v.*v.*TY.*imfilter(vy1,ddy1);
    
    p1(z0,x0,y0)=p1(z0,x0,y0)+s(k);
    p0=p1;
    vz0=vz1;vx0=vx1;vy0=vy1;

    if mod(t,0.5)==0
        figure();
        imagesc(p1(:,:,y0)); 
    end
end
