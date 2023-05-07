clc,clear
% 3D - Acoustic wave simulation using finite difference. 
% with staggered-grid; with Convolutional Perfectly Matched Layer;
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
v=ones(nz,nx,ny)*2000;
z0=round(nz/2);x0=round(nx/2);y0=round(ny/2);
vx0=zeros(nz,nx,ny);vx1=vx0;
vz0=vx0;vz1=vx0;
vy0=vx0;vy1=vx0;
p0=vx0;p1=vx0;
EXX0=vx0;EXX1=vx0;
EYY0=vx0;EYY1=vx0;
EZZ0=vx0;EZZ1=vx0;
HXX0=vx0;HXX1=vx0;
HYY0=vx0;HYY1=vx0;
HZZ0=vx0;HZZ1=vx0;
[aaz,bbz,aax,bbx,aay,bby] = cpml(1e-10,[nz nx ny],20,max(v(:)),dt,dz);

nn=3;
dxd=FDcoeffDx(nn);ddz0(:,1,1)=dxd';ddz1(:,1,1)=[dxd 0]';
ddx0=permute(ddz0,[2 1 3]);ddx1=permute(ddz1,[2 1 3]);ddy0=permute(ddz0,[3 2 1]);ddy1=permute(ddz1,[3 2 1]);
TX=dt/dx;TZ=dt/dz;TY=dt/dy;
for	t=DT:DT:T
    disp(t);
    k=round(t/DT);
    EXX1=bbx.*EXX0+aax.*TX.*imfilter(p0,ddx0)/dt;
    EYY1=bby.*EYY0+aay.*TY.*imfilter(p0,ddy0)/dt;
    EZZ1=bbz.*EZZ0+aaz.*TZ.*imfilter(p0,ddz0)/dt;
    vz1=vz0+TZ*imfilter(p0,ddz0)+EZZ1*dt;
    vx1=vx0+TX*imfilter(p0,ddx0)+EXX1*dt;
    vy1=vy0+TY*imfilter(p0,ddy0)+EYY1*dt;
    
    HXX1=bbx.*HXX0+aax.*TX.*imfilter(vx1,ddx1)/dt;
    HYY1=bby.*HYY0+aay.*TY.*imfilter(vy1,ddy1)/dt;
    HZZ1=bbz.*HZZ0+aaz.*TZ.*imfilter(vz1,ddz1)/dt;
    p1=p0+v.*v.*TZ.*imfilter(vz1,ddz1)+v.*v.*TX.*imfilter(vx1,ddx1)+v.*v.*TY.*imfilter(vy1,ddy1)+v.*v.*(HXX1*dt+HZZ1*dt+HYY1*dt);
    
    p1(z0,x0,y0)=p1(z0,x0,y0)+s(k);
    p0=p1;
    vz0=vz1;vx0=vx1;vy0=vy1;
    EXX0=EXX1;EYY0=EYY1;EZZ0=EZZ1;
    HXX0=HXX1;HYY0=HYY1;HZZ0=HZZ1;
    
    if mod(t,0.5)==0
        figure();
        imagesc(p1(:,:,y0)); 
    end
end
