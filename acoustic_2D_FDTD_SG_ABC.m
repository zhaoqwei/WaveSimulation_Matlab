clc,clear
% Acoustic wave simulation using finite difference. 
% with staggered-grid; with sponge absorbing boundary condition;
% By zhaoqingwei
% Chengdu University of Technology (CDUT), 2021-2025
%%%%%%%%%%%%%%%%
% P--Vx--P
% |      |
% Vz-----vz
% |      |
% P--Vx--P
%%%%%%%%%%%%%%%

FM=20;
DT=0.001;dt=DT;
T=4;
nx=600;
nz=600;
DH=10;dx=DH;dz=DH;

s=wavelet(FM,DT,T);
v=ones(nz,nx)*2000;
z0=round(nz/2);x0=round(nx/2);

vx0=zeros(nz,nx);vx1=vx0;
vz0=vx0;vz1=vx0;
p0=vx0;p1=vx0;

[aaz,bbz,aax,bbx] = cpml(1e-10,[nz nx],20,max(v(:)),dt,dz);

nn=3;
dxd=FDcoeffDx(nn);ddz0=dxd';ddz1=[dxd 0]';ddx0=ddz0';ddx1=ddz1';
TX=dt/dx;TZ=dt/dz;
for	t=DT:DT:T
    k=round(t/DT);
    vz1=vz0+TZ*imfilter(p0,ddz0);
    vx1=vx0+TX*imfilter(p0,ddx0);
    vz1=vz1.*bbz;
    vx1=vx1.*bbx;
    p1=p0+v.*v.*TZ.*imfilter(vz1,ddz1)+v.*v.*TX.*imfilter(vx1,ddx1);
    p1=p1.*bbz.*bbx;
    
    p1(z0,x0)=p1(z0,x0)+s(k);
    p0=p1;
    vz0=vz1;vx0=vx1;
    
    if mod(t,1)==0
        figure();
        imagesc(p1);
    end
end
