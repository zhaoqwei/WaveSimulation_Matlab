clc,clear
% 2D BIOT Elastic wave simulation using finite difference in poroelastic media
% with staggered-grid; with Convolutional Perfectly Matched Layer;
% By zhaoqingwei
% Chengdu University of Technology (CDUT), 2021-2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%txx(tzz)--Vx--txx(tzz)
% |        |    |
% Vz------txz---vz
% |        |    |
%txx(tzz)--Vx--txx(tzz)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FM=20;
DT=0.001;dt=DT;
T=2;
nx=600;
nz=600;
DH=10;dx=DH;dz=DH;

s=wavelet(FM,DT,T);
vp=ones(nz,nx)*4000;vs=vp/2;rou=ones(600,600)*2;
lamda=rou.*(vp.*vp-vs.*vs*2);mu=rou.*vs.*vs;
c11=rou.*vp.*vp;c13=rou.*(vp.*vp-vs.*vs*2);c33=c11;c55=rou.*vs.*vs;%各向同性
% Epsilon=0.25*ones(600,600);Delta=0.1*ones(600,600);
% c33=vp.*vp.*rou;c55=vs.*vs.*rou;c11=2*c33.*Epsilon+c33;c13=rou.*sqrt(((1+2*Delta).*vp.*vp-vs.*vs).*(vp.*vp-vs.*vs))-rou.*vs.*vs;%VTI

%  BIOT
A=12.72*1e9;N=6.84*1e9; Ro11=2.17*1e3;
R=0.331*1e9;Ro22=0.191*1e3;
Q=0.953*1e9;Ro12=-0.083*1e3;
b=3.0*1e-2;
Ro=(Ro11*Ro22-Ro12*Ro12);
R11=Ro11/Ro;R12=Ro12/Ro;R22=Ro22/Ro;
B1=b/(R11+R12);B2=b/(R12+R22);
c11=(A+2*N)*ones(nz,nx);c33=c11;
c13=(A)*ones(nz,nx);
c55=(N)*ones(nz,nx);
R11=R11*ones(nz,nx);
R12=R12*ones(nz,nx);
R22=R22*ones(nz,nx);
B1=B1*ones(nz,nx);B2=B2*ones(nz,nx);


vx0=zeros(nz,nx);vx1=vx0;
vz0=vx0;vz1=vx0;
taoxx0=vx0;taoxx1=vx0;
taozz0=vx0;taozz1=vx0;
taoxz0=vx0;taoxz1=vx0;
VX0=vx0;VX1=vx0;
VZ0=vx0;VZ1=vx0;
sf0=vx0;sf1=vx0;
z0=round(nz/2);x0=round(nx/2);

EXX0=vx0;EXX1=vx0;
EXZ0=vx0;EXZ1=vx0;
EZX0=vx0;EZX1=vx0;
EZZ0=vx0;EZZ1=vx0;
HXX0=vx0;HXX1=vx0;
HXZ0=vx0;HXZ1=vx0;
HZX0=vx0;HZX1=vx0;
HZZ0=vx0;HZZ1=vx0;
FX0=vx0;FX1=vx0;
FZ0=vx0;FZ1=vx0;
GX0=vx0;GX1=vx0;
GZ0=vx0;GZ1=vx0;

[aaz,bbz,aax,bbx] = cpml(1e-5,[nz nx],20,max(vp(:)),dt,dz);

nn=3;
dxd=FDcoeffDx(nn);ddz0=dxd';ddz1=[dxd 0]';ddx0=ddz0';ddx1=ddz1';
TX=dt/dx;TZ=dt/dz;
for	t=DT:DT:T
    disp(t);
    k=round(t/DT);
    EXX1=bbx.*EXX0+aax.*TX.*imfilter(taoxx0,ddx0)/dt;
    EXZ1=bbz.*EXZ0+aaz.*TZ.*imfilter(taoxz0,ddz1)/dt;
    EZX1=bbx.*EZX0+aax.*TX.*imfilter(taoxz0,ddx1)/dt;
    EZZ1=bbz.*EZZ0+aaz.*TZ.*imfilter(taozz0,ddz0)/dt;
    FX1=bbx.*FX0+aax.*TX.*imfilter(sf0,ddx0)/dt;
    FZ1=bbz.*FZ0+aaz.*TZ.*imfilter(sf0,ddz0)/dt;
    vx1=vx0+R22.*(TZ.*imfilter(taoxz0,ddz1)+TX.*imfilter(taoxx0,ddx0)+EXZ1*dt+EXX1*dt)...
        -R12.*(TX.*imfilter(sf0,ddx0)+FX1*dt)+B2.*(VX0-vx0)*dt;
    vz1=vz0+R22.*(TZ.*imfilter(taozz0,ddz0)+TX.*imfilter(taoxz0,ddx1)+EZZ1*dt+EZX1*dt)...
        -R12.*(TZ.*imfilter(sf0,ddz0)+FZ1*dt)+B2.*(VZ0-vz0)*dt;
    VX1=VX0-R12.*(TZ.*imfilter(taoxz0,ddz1)+TX.*imfilter(taoxx0,ddx0)+EXZ1*dt+EXX1*dt)...
        +R11.*(TX.*imfilter(sf0,ddx0)+FX1*dt)-B1.*(VX0-vx0)*dt;
    VZ1=VZ0-R12.*(TZ.*imfilter(taozz0,ddz0)+TX.*imfilter(taoxz0,ddx1)+EZZ1*dt+EZX1*dt)...
        +R11.*(TZ.*imfilter(sf0,ddz0)+FZ1*dt)-B1.*(VZ0-vz0)*dt;
    
    vz1(z0,x0)=vz1(z0,x0)+s(k);
    
    HXX1=bbx.*HXX0+aax.*TX.*imfilter(vx1,ddx1)/dt;
    HXZ1=bbz.*HXZ0+aaz.*TZ.*imfilter(vx1,ddz0)/dt;
    HZX1=bbx.*HZX0+aax.*TX.*imfilter(vz1,ddx0)/dt;
    HZZ1=bbz.*HZZ0+aaz.*TZ.*imfilter(vz1,ddz1)/dt;
    GX1=bbx.*GX0+aax.*TX.*imfilter(VX1,ddx1)/dt;
    GZ1=bbz.*GZ0+aaz.*TZ.*imfilter(VZ1,ddz1)/dt;
    taoxx1=taoxx0+c13.*TZ.*imfilter(vz1,ddz1)+c11.*TX.*imfilter(vx1,ddx1)+c13.*HZZ1*dt+c11.*HXX1*dt+...
        Q.*(TZ.*imfilter(VZ1,ddz1)+TX.*imfilter(VX1,ddx1)+GX1*dt+GZ1*dt);
    taozz1=taozz0+c33.*TZ.*imfilter(vz1,ddz1)+c13.*TX.*imfilter(vx1,ddx1)+c33.*HZZ1*dt+c13.*HXX1*dt+...
        Q.*(TZ.*imfilter(VZ1,ddz1)+TX.*imfilter(VX1,ddx1)+GX1*dt+GZ1*dt);
    taoxz1=taoxz0+c55.*TX.*imfilter(vx1,ddz0)+c55.*TZ.*imfilter(vz1,ddx0)+c55.*HXZ1*dt+c55.*HZX1*dt;
    sf1=sf0+Q.*(TZ.*imfilter(vz1,ddz1)+TX.*imfilter(vx1,ddx1)+HZZ1*dt+HXX1*dt)+R.*(TZ.*imfilter(VZ1,ddz1)+TX.*imfilter(VX1,ddx1)+GX1*dt+GZ1*dt);
    
%     taoxx1(z0,x0)=taoxx1(z0,x0)+s(k);
%     taozz1(z0,x0)=taozz1(z0,x0)+s(k);
    
    EXX0=EXX1;EXZ0=EXZ1;EZX0=EZX1;EZZ0=EZZ1;
    HXX0=HXX1;HXZ0=HXZ1;HZX0=HZX1;HZZ0=HZZ1;
    FX0=FX1;FZ0=FZ1;GX0=GX1;GZ0=GZ1;
    taoxx0=taoxx1;taozz0=taozz1;taoxz0=taoxz1;
    vz0=vz1;vx0=vx1;
    VZ0=VZ1;VX0=VX1;sf0=sf1;
    
    if mod(t,0.5)==0
        figure();
        imagesc(VZ1);
    end
end
