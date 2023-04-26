clc,clear
% 3D - Elastic wave simulation using finite difference. 
% with staggered-grid; with Convolutional Perfectly Matched Layer;
% By zhaoqingwei
% Chengdu University of Technology (CDUT), 2021-2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%grid point
%               |vx  |txx      |vx 
%          |txy   |vy     |txy  
%     vx  -- txx --  vx    |   |txz
%      |       | /tyz|     | /
%     txz  --  vz -- txz  /    |v
%      |       |     |     |txy
%     vx  --  txx  -- vx  / 
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


FM=20;
DT=0.001;dt=DT;
T=2;
nx=200;
ny=200;
nz=200;
DH=10;dx=DH;dz=DH;dy=DH;

s=wavelet(FM,DT,T);z0=round(nz/2);x0=round(nx/2);y0=round(ny/2);
vp=ones(nz,nx,ny)*2000;vs=vp/2;rou=ones(nz,nx,ny)*2;
Epsilon=0.1*ones(nz,nx,ny);Delta=0.1*ones(nz,nx,ny);Gamma=0.2*ones(nz,nx,ny);
c11=rou.*vp.*vp;c22=c11;c33=c11;c55=rou.*vs.*vs;c44=c55;c66=c55;c13=rou.*(vp.*vp-vs.*vs*2);c12=c13;c23=c13;%各向同性
% c33=vp.*vp.*rou;c11=2*c33.*Epsilon+c33;c22=c11;
% c55=vs.*vs.*rou;c44=c55;c66=2*c44.*Gamma+c44;
% c13=rou.*sqrt(((1+2*Delta).*vp.*vp-vs.*vs).*(vp.*vp-vs.*vs))-rou.*vs.*vs;c23=c13;c12=c11-2*c66;%VTI

vx0=zeros(nz,nx,ny);vx1=vx0;
vz0=vx0;vz1=vx0;
vy0=vx0;vy1=vx0;
taoxx0=vx0;taoxx1=vx0;
taoyy0=vx0;taoyy1=vx0;
taozz0=vx0;taozz1=vx0;
taoxz0=vx0;taoxz1=vx0;
taoyz0=vx0;taoyz1=vx0;
taoxy0=vx0;taoxy1=vx0;

EXX0=vx0;EXX1=vx0;
EXY0=vx0;EXY1=vx0;
EXZ0=vx0;EXZ1=vx0;
EYX0=vx0;EYX1=vx0;
EYY0=vx0;EYY1=vx0;
EYZ0=vx0;EYZ1=vx0;
EZX0=vx0;EZX1=vx0;
EZY0=vx0;EZY1=vx0;
EZZ0=vx0;EZZ1=vx0;

HXX0=vx0;HXX1=vx0;
HXY0=vx0;HXY1=vx0;
HXZ0=vx0;HXZ1=vx0;
HYX0=vx0;HYX1=vx0;
HYY0=vx0;HYY1=vx0;
HYZ0=vx0;HYZ1=vx0;
HZX0=vx0;HZX1=vx0;
HZY0=vx0;HZY1=vx0;
HZZ0=vx0;HZZ1=vx0;

[aaz,bbz,aax,bbx,aay,bby] = cpml(1e-10,[nz nx ny],20,max(vp(:)),dt,dz);

nn=2;
dxd=FDcoeffDx(nn);ddz0(:,1,1)=dxd';ddz1(:,1,1)=[dxd 0]';
ddx0=permute(ddz0,[2 1 3]);ddx1=permute(ddz1,[2 1 3]);ddy0=permute(ddz0,[3 2 1]);ddy1=permute(ddz1,[3 2 1]);
TX=dt/dx;TZ=dt/dz;TY=dt/dy;
for	t=DT:DT:T
    disp(t);
    k=round(t/DT);
    EXX1=bbx.*EXX0+aax.*TX.*imfilter(taoxx0,ddx1)/dt;
    EXY1=bby.*EXY0+aay.*TY.*imfilter(taoxy0,ddy1)/dt;
    EXZ1=bbz.*EXY0+aaz.*TZ.*imfilter(taoxz0,ddz1)/dt;
    EYX1=bbx.*EYX0+aax.*TX.*imfilter(taoxy0,ddx0)/dt;
    EYY1=bby.*EYY0+aay.*TY.*imfilter(taoyy0,ddy0)/dt;
    EYZ1=bbz.*EYY0+aaz.*TZ.*imfilter(taoyz0,ddz1)/dt;
    EZX1=bbx.*EZX0+aax.*TX.*imfilter(taoxz0,ddx0)/dt;
    EZY1=bby.*EZY0+aay.*TY.*imfilter(taoyz0,ddy1)/dt;
    EZZ1=bbz.*EZY0+aaz.*TZ.*imfilter(taozz0,ddz0)/dt;
    
    vx1=vx0+(TX.*imfilter(taoxx0,ddx1)+TY.*imfilter(taoxy0,ddy1)+TZ.*imfilter(taoxz0,ddz1)+dt*EXX1+dt*EXY1+dt*EXZ1)./rou;
    vy1=vy0+(TX.*imfilter(taoxy0,ddx0)+TY.*imfilter(taoyy0,ddy0)+TZ.*imfilter(taoyz0,ddz1)+dt*EYX1+dt*EYY1+dt*EYZ1)./rou;
    vz1=vz0+(TX.*imfilter(taoxz0,ddx0)+TY.*imfilter(taoyz0,ddy1)+TZ.*imfilter(taozz0,ddz0)+dt*EZX1+dt*EZY1+dt*EZZ1)./rou;
    
    vz1(z0,x0,y0)=vz1(z0,x0,y0)+s(k);

    HXX1=bbx.*HXX0+aax.*TX.*imfilter(vx1,ddx0)/dt;
    HXY1=bby.*HXY0+aay.*TY.*imfilter(vx1,ddy0)/dt;
    HXZ1=bbz.*HXY0+aaz.*TZ.*imfilter(vx1,ddz0)/dt;
    HYX1=bbx.*HYX0+aax.*TX.*imfilter(vy1,ddx1)/dt;
    HYY1=bby.*HYY0+aay.*TY.*imfilter(vy1,ddy1)/dt;
    HYZ1=bbz.*HYY0+aaz.*TZ.*imfilter(vy1,ddz0)/dt;
    HZX1=bbx.*HZX0+aax.*TX.*imfilter(vz1,ddx1)/dt;
    HZY1=bby.*HZY0+aay.*TY.*imfilter(vz1,ddy0)/dt;
    HZZ1=bbz.*HZY0+aaz.*TZ.*imfilter(vz1,ddz1)/dt;
    
    taoxx1=taoxx0+c11.*TX.*imfilter(vx1,ddx0)+c12.*TY.*imfilter(vy1,ddy1)+c13.*TZ.*imfilter(vz1,ddz1)+c11.*HXX1*dt+c12.*HYY1*dt+c13.*HZZ1*dt;
    taoyy1=taoyy0+c12.*TX.*imfilter(vx1,ddx0)+c22.*TY.*imfilter(vy1,ddy1)+c23.*TZ.*imfilter(vz1,ddz1)+c12.*HXX1*dt+c22.*HYY1*dt+c23.*HZZ1*dt;
    taozz1=taozz0+c13.*TX.*imfilter(vx1,ddx0)+c23.*TY.*imfilter(vy1,ddy1)+c33.*TZ.*imfilter(vz1,ddz1)+c13.*HXX1*dt+c23.*HYY1*dt+c33.*HZZ1*dt;
    taoyz1=taoyz0                            +c44.*TZ.*imfilter(vy1,ddz0)+c44.*TY.*imfilter(vz1,ddy0)+c44.*HYZ1*dt+c44.*HZY1*dt;
    taoxz1=taoxz0+c55.*TZ.*imfilter(vx1,ddz0)                             +c55.*TX.*imfilter(vz1,ddx1)+c55.*HXZ1*dt+c55.*HZX1*dt;
    taoxy1=taoxy0+c66.*TY.*imfilter(vx1,ddy0)+c66.*TX.*imfilter(vy1,ddx1)                          +c66.*HXY1*dt+c66.*HYX1*dt;
    
    EXX0=EXX1;EXY0=EXY1;EXZ0=EXZ1;
    EYX0=EYX1;EYY0=EYY1;EYZ0=EYZ1;
    EZX0=EZX1;EZY0=EZY1;EZZ0=EZZ1;

    HXX0=HXX1;HXY0=HXY1;HXZ0=HXZ1;
    HYX0=HYX1;HYY0=HYY1;HYZ0=HYZ1;
    HZX0=HZX1;HZY0=HZY1;HZZ0=HZZ1;
    vx0=vx1;
    vy0=vy1;
    vz0=vz1;
    taoxx0=taoxx1;
    taoyy0=taoyy1;
    taozz0=taozz1;
    taoxz0=taoxz1;
    taoxy0=taoxy1;
    taoyz0=taoyz1;

    if mod(t,0.5)==0
        figure();
        imagesc(vz1(:,:,y0)); 
    end
end
