clc,clear
% 2D Elastic wave simulation using finite difference. 
% Up layer with free boundary generate surfce wave.
% with staggered-grid; with Convolutional Perfectly Matched Layer;with free Boundary;
% By zhaoqingwei
% Chengdu University of Technology (CDUT), 2021-2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      vz---txz---vz
%      |     |    |
%  txx（tzz）vx   ...
%      |     |    |
%      vz --txz-- ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FM=20;
DT=0.001;dt=DT;
T=2;
nx=600;
nz=600;
DH=10;dx=DH;dz=DH;

s=wavelet(FM,DT,T);
vp=ones(nz,nx)*2000;vs=vp/2;rou=ones(600,600)*2;
lamda=rou.*(vp.*vp-vs.*vs*2);mu=rou.*vs.*vs;
c11=rou.*vp.*vp;c13=rou.*(vp.*vp-vs.*vs*2);c33=c11;c55=rou.*vs.*vs;%各向同性
% Epsilon=0.25*ones(600,600);Delta=0.1*ones(600,600);
% c33=vp.*vp.*rou;c55=vs.*vs.*rou;c11=2*c33.*Epsilon+c33;c13=rou.*sqrt(((1+2*Delta).*vp.*vp-vs.*vs).*(vp.*vp-vs.*vs))-rou.*vs.*vs;%VTI

vx0=zeros(nz,nx);vx1=vx0;
vz0=vx0;vz1=vx0;
taoxx0=vx0;taoxx1=vx0;
taozz0=vx0;taozz1=vx0;
taoxz0=vx0;taoxz1=vx0;
z0=round(3);x0=round(nx/2);

EXX0=vx0;EXX1=vx0;
EXZ0=vx0;EXZ1=vx0;
EZX0=vx0;EZX1=vx0;
EZZ0=vx0;EZZ1=vx0;
HXX0=vx0;HXX1=vx0;
HXZ0=vx0;HXZ1=vx0;
HZX0=vx0;HZX1=vx0;
HZZ0=vx0;HZZ1=vx0;
cpmln=20;[aaz,bbz,aax,bbx] = cpml(1e-10,[nz nx],cpmln,max(vp(:)),dt,dz); aaz(1:cpmln,:)=0;bbz(1:cpmln,:)=1;

nn=3;
dxd=FDcoeffDx(nn);ddz0=dxd';ddz1=[dxd 0]';ddx0=ddz0';ddx1=ddz1';
TX=dt/dx;TZ=dt/dz;
for	t=DT:DT:T
    disp(t);
    k=round(t/DT);
    EXX1=bbx.*EXX0+aax.*TX.*imfilter(taoxx0,ddx0)/dt;
    EXZ1=bbz.*EXZ0+aaz.*TZ.*imfilter(taoxz0,ddz0)/dt;
    EZX1=bbx.*EZX0+aax.*TX.*imfilter(taoxz0,ddx1)/dt;
    EZZ1=bbz.*EZZ0+aaz.*TZ.*imfilter(taozz0,ddz1)/dt;
    vx1=vx0+(TZ.*imfilter(taoxz0,ddz0)+TX.*imfilter(taoxx0,ddx0)+EXZ1*dt+EXX1*dt)./rou;
    vz1=vz0+(TZ.*imfilter(taozz0,ddz1)+TX.*imfilter(taoxz0,ddx1)+EZZ1*dt+EZX1*dt)./rou;
    
    vz1(z0,x0)=vz1(z0,x0)+s(k);
    vx1(1:2,:)=0;
    vz1(1:2,:)=0;
    vx1(1,:)=vx1(3,:)+2*imfilter(vz1(3,:),ddx0);
    vz1(1,:)=vz1(3,:)+2*(1-2*(c55(3,:)./c11(3,:))).*imfilter(vx1(1,:),ddx1);
    vx1(2,:)=vx1(1,:)+2*imfilter(vz1(3,:),ddx0);
    vz1(2,:)=vz1(1,:)+2*(1-2*(c55(3,:)./c11(3,:))).*imfilter(vx1(2,:),ddx1);
    
    HXX1=bbx.*HXX0+aax.*TX.*imfilter(vx1,ddx1)/dt;
    HXZ1=bbz.*HXZ0+aaz.*TZ.*imfilter(vx1,ddz1)/dt;
    HZX1=bbx.*HZX0+aax.*TX.*imfilter(vz1,ddx0)/dt;
    HZZ1=bbz.*HZZ0+aaz.*TZ.*imfilter(vz1,ddz0)/dt;
    taoxx1=taoxx0+c13.*TZ.*imfilter(vz1,ddz0)+c11.*TX.*imfilter(vx1,ddx1)+c13.*HZZ1*dt+c11.*HXX1*dt;
    taozz1=taozz0+c33.*TZ.*imfilter(vz1,ddz0)+c13.*TX.*imfilter(vx1,ddx1)+c33.*HZZ1*dt+c13.*HXX1*dt;
    taoxz1=taoxz0+c55.*TZ.*imfilter(vx1,ddz1)+c55.*TX.*imfilter(vz1,ddx0)+c55.*HXZ1*dt+c55.*HZX1*dt;
    
    taoxx1(1:2,:)=0;
    taozz1(2,:)=-taozz1(3,:);taozz1(1,:)=-taozz1(4,:);
    taoxz1(1:3,:)=0;

    
    taoxx0=taoxx1;taozz0=taozz1;taoxz0=taoxz1;
    vz0=vz1;vx0=vx1;
    
    if mod(t,0.5)==0
        figure();
        imagesc(vz1);
    end
end
