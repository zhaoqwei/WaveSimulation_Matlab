clc,clear
% ViscoAcoustic wave simulation using finite difference. with staggered-grid.
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
T=2;
nx=600;
nz=600;
DH=10;dx=DH;dz=DH;

s=wavelet(FM,DT,T)*1e5;
v=ones(nz,nx)*4000;Q=ones(nz,nx)*100;
z0=round(nz/2);x0=round(nx/2);

vx0=zeros(nz,nx);vx1=vx0;
vz0=vx0;vz1=vx0;
p0=vx0;p1=vx0;

EXX0=vx0;EXX1=vx0;
EZZ0=vx0;EZZ1=vx0;
HXX0=vx0;HXX1=vx0;
HZZ0=vx0;HZZ1=vx0;
[aaz,bbz,aax,bbx] = cpml(1e-6,[nz nx],20,max(v(:)),dt,dz);

nn=3;
dxd=FDcoeffDx(nn);ddz0=dxd';ddz1=[dxd 0]';ddx0=ddz0';ddx1=ddz1';
TX=dt/dx;TZ=dt/dz;
for	t=DT:DT:T
    disp(t);
    k=round(t/DT);
    EXX1=bbx.*EXX0+aax.*TX.*imfilter(p0,ddx0)/dt;
    EZZ1=bbz.*EZZ0+aaz.*TZ.*imfilter(p0,ddz0)/dt;
    vz1=vz0+TZ*imfilter(p0,ddz0)+EZZ1*dt;
    vx1=vx0+TX*imfilter(p0,ddx0)+EXX1*dt;
    
    HXX1=bbx.*HXX0+aax.*TX.*imfilter(vx1,ddx1)/dt;
    HZZ1=bbz.*HZZ0+aaz.*TZ.*imfilter(vz1,ddz1)/dt;
    p1=p0+(v.*v+1i*v.*v./Q).*(TZ.*imfilter(vz1,ddz1)+TX.*imfilter(vx1,ddx1)+(HXX1*dt+HZZ1*dt));
    
    EXX0=EXX1;EZZ0=EZZ1;
    HXX0=HXX1;HZZ0=HZZ1;
    p1(z0,x0)=p1(z0,x0)+s(k);
    p0=p1;
    vz0=vz1;vx0=vx1;
    
    if mod(t,0.5)==0
        figure();
        imagesc(real(p1));
    end
end
