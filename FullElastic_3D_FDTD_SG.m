clc,clear
% 3D - full elastic wave simulation using finite difference. with staggered-grid.
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
c33=vp.*vp.*rou;c11=2*c33.*Epsilon+c33;c22=c11;
c55=vs.*vs.*rou;c44=c55;c66=2*c44.*Gamma+c44;
c13=rou.*sqrt(((1+2*Delta).*vp.*vp-vs.*vs).*(vp.*vp-vs.*vs))-rou.*vs.*vs;c23=c13;c12=c11-2*c66;%VTI

c_old=[c11(1) c12(1) c13(1) 0 0 0;
   c12(1) c22(1) c23(1) 0 0 0;
   c13(1) c23(1) c33(1) 0 0 0;
    0  0  0  c44(1) 0 0;
    0  0  0  0  c55(1) 0;
    0  0  0  0  0  c66(1);];
M=bond([0 0 pi/4]);
c_new= M*c_old*M';
c11=ones(nz,nx,ny)*c_new(1,1);c12=ones(nz,nx,ny)*c_new(1,2);c13=ones(nz,nx,ny)*c_new(1,3);
c14=ones(nz,nx,ny)*c_new(1,4);c15=ones(nz,nx,ny)*c_new(1,5);c16=ones(nz,nx,ny)*c_new(1,6);
c22=ones(nz,nx,ny)*c_new(2,2);c23=ones(nz,nx,ny)*c_new(2,3);
c24=ones(nz,nx,ny)*c_new(2,4);c25=ones(nz,nx,ny)*c_new(2,5);c26=ones(nz,nx,ny)*c_new(2,6);
c33=ones(nz,nx,ny)*c_new(3,3);c34=ones(nz,nx,ny)*c_new(3,4);c35=ones(nz,nx,ny)*c_new(3,5);c36=ones(nz,nx,ny)*c_new(3,6);
c44=ones(nz,nx,ny)*c_new(4,4);c45=ones(nz,nx,ny)*c_new(4,5);c46=ones(nz,nx,ny)*c_new(4,6);
c55=ones(nz,nx,ny)*c_new(5,5);c56=ones(nz,nx,ny)*c_new(5,6);c66=ones(nz,nx,ny)*c_new(6,6);


vx0=zeros(nz,nx,ny);vx1=vx0;
vz0=vx0;vz1=vx0;
vy0=vx0;vy1=vx0;
taoxx0=vx0;taoxx1=vx0;
taoyy0=vx0;taoyy1=vx0;
taozz0=vx0;taozz1=vx0;
taoxz0=vx0;taoxz1=vx0;
taoyz0=vx0;taoyz1=vx0;
taoxy0=vx0;taoxy1=vx0;

nn=3;
dxd=FDcoeffDx(nn);ddz0(:,1,1)=dxd';ddz1(:,1,1)=[dxd 0]';
ddx0=permute(ddz0,[2 1 3]);ddx1=permute(ddz1,[2 1 3]);ddy0=permute(ddz0,[3 2 1]);ddy1=permute(ddz1,[3 2 1]);
Xdx0=[-0.25 0 0.25;-0.25 0 0.25];Xdx1=[Xdx0;zeros(1,3)];Xdz0=Xdx0';Xdz1=Xdx1';Xdy0=permute(Xdx0,[1 3 2]);Xdy1=permute(Xdx1,[1 3 2]);
TX=dt/dx;TZ=dt/dz;TY=dt/dy;

for	t=DT:DT:T
    disp(t);
    k=round(t/DT);

    vx1=vx0+(TX.*imfilter(taoxx0,ddx1)+TY.*imfilter(taoxy0,ddy1)+TZ.*imfilter(taoxz0,ddz1))./rou;
    vy1=vy0+(TX.*imfilter(taoxy0,ddx0)+TY.*imfilter(taoyy0,ddy0)+TZ.*imfilter(taoyz0,ddz1))./rou;
    vz1=vz0+(TX.*imfilter(taoxz0,ddx0)+TY.*imfilter(taoyz0,ddy1)+TZ.*imfilter(taozz0,ddz0))./rou;
       
    vz1(z0,x0,y0)=vz1(z0,x0,y0)+s(k);

    taoxx1=taoxx0+c11.*TX.*imfilter(vx1,ddx0)+c16.*TY.*imfilter(vx1,Xdy1)+c15.*TZ.*imfilter(vx1,Xdz1)+...
        c12.*TY.*imfilter(vy1,ddy1)+c16.*TX.*imfilter(vy1,Xdx1)+c14.*TZ.*imfilter(vy1,Xdz1)+...
        c13.*TZ.*imfilter(vz1,ddz1)+c15.*TX.*imfilter(vz1,Xdx1)+c14.*TZ.*imfilter(vz1,Xdy1);
    taoyy1=taoyy0+c12.*TX.*imfilter(vx1,ddx0)+c26.*TY.*imfilter(vx1,Xdy1)+c25.*TZ.*imfilter(vx1,Xdz1)+...
        c22.*TY.*imfilter(vy1,ddy1)+c26.*TX.*imfilter(vy1,Xdx1)+c24.*TZ.*imfilter(vy1,Xdz1)+...
        c23.*TZ.*imfilter(vz1,ddz1)+c25.*TX.*imfilter(vz1,Xdx1)+c24.*TZ.*imfilter(vz1,Xdy1);
    taozz1=taozz0+c13.*TX.*imfilter(vx1,ddx0)+c36.*TY.*imfilter(vx1,Xdy1)+c35.*TZ.*imfilter(vx1,Xdz1)+...
        c23.*TY.*imfilter(vy1,ddy1)+c36.*TX.*imfilter(vy1,Xdx1)+c34.*TZ.*imfilter(vy1,Xdz1)+...
        c33.*TZ.*imfilter(vz1,ddz1)+c35.*TX.*imfilter(vz1,Xdx1)+c34.*TZ.*imfilter(vz1,Xdy1);
    taoyz1=taoyz0+c14.*TX.*imfilter(vx1,Xdx1)+c46.*TY.*imfilter(vx1,Xdy1)+c45.*TZ.*imfilter(vx1,Xdz1)+...
        c44.*TZ.*imfilter(vy1,ddz0)+c46.*TY.*imfilter(vy1,Xdx1)+c24.*TZ.*imfilter(vy1,Xdy1)+...
        c44.*TY.*imfilter(vz1,ddy0)+c45.*TY.*imfilter(vz1,Xdx1)+c34.*TZ.*imfilter(vz1,Xdz1);
    taoxz1=taoxz0+c55.*TZ.*imfilter(vx1,ddz0)+c15.*TY.*imfilter(vx1,Xdx1)+c56.*TZ.*imfilter(vx1,Xdy1)+...
        c56.*TX.*imfilter(vy1,Xdx1)+c25.*TY.*imfilter(vy1,Xdy1)+c45.*TZ.*imfilter(vy1,Xdz1)+...
        c55.*TX.*imfilter(vz1,ddx1)+c45.*TY.*imfilter(vz1,Xdy1)+c35.*TZ.*imfilter(vz1,Xdz1);
    taoxy1=taoxy0+c66.*TY.*imfilter(vx1,ddy0)+c16.*TY.*imfilter(vx1,Xdx1)+c56.*TZ.*imfilter(vx1,Xdz1)+...
        +c66.*TX.*imfilter(vy1,ddx1)+c26.*TY.*imfilter(vy1,Xdy1)+c46.*TZ.*imfilter(vy1,Xdz1)+...
        c56.*TX.*imfilter(vz1,Xdx1)+c46.*TY.*imfilter(vz1,Xdy1)+c36.*TZ.*imfilter(vz1,Xdz1);
    
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
