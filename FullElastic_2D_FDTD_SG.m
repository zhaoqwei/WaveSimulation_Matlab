clc,clear
% 2D - full elastic  wave simulation using finite difference. with staggered-grid.
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
vp=ones(nz,nx)*4000;vs=vp/2;rou=ones(nz,nx)*2;
lamda=rou.*(vp.*vp-vs.*vs*2);mu=rou.*vs.*vs;
Epsilon=0.2*ones(600,600);Delta=0.1*ones(600,600);
c33=vp.*vp.*rou;c55=vs.*vs.*rou;c11=2*c33.*Epsilon+c33;c13=rou.*sqrt(((1+2*Delta).*vp.*vp-vs.*vs).*(vp.*vp-vs.*vs))-rou.*vs.*vs;c15=zeros(nz,nx);c35=zeros(nz,nx);%VTI
c_old=[c11(1) c13(1) c15(1);
   c13(1) c33(1) c35(1);
   c15(1) c35(1) c55(1);];
M=bond([0 0 pi/4]);M1=M(1:2:5,1:2:5);
c_new= M1*c_old*M1';
c11=ones(nz,nx)*c_new(1,1);c13=ones(nz,nx)*c_new(1,2);c15=ones(nz,nx)*c_new(1,3);
c33=ones(nz,nx)*c_new(2,2);c35=ones(nz,nx)*c_new(2,3);c55=ones(nz,nx)*c_new(3,3);

vx0=zeros(nz,nx);vx1=vx0;
vz0=vx0;vz1=vx0;
taoxx0=vx0;taoxx1=vx0;
taozz0=vx0;taozz1=vx0;
taoxz0=vx0;taoxz1=vx0;
z0=round(nz/2);x0=round(nx/2);

nn=3;
dxd=FDcoeffDx(nn);ddz0=dxd';ddz1=[dxd 0]';ddx0=ddz0';ddx1=ddz1';Xdx0=[-0.25 0 0.25;-0.25 0 0.25];Xdx1=[Xdx0;zeros(1,3)];Xdz0=Xdx0';Xdz1=Xdx1';
TX=dt/dx;TZ=dt/dz;
for	t=DT:DT:T
     disp(t);
    k=round(t/DT);
    vx1=vx0+(TZ.*imfilter(taoxz0,ddz1)+TX.*imfilter(taoxx0,ddx0))./rou;%
    vz1=vz0+(TZ.*imfilter(taozz0,ddz0)+TX.*imfilter(taoxz0,ddx1))./rou;
    
    vz1(z0,x0)=vz1(z0,x0)+s(k);
    
    taoxx1=taoxx0+c13.*TZ.*imfilter(vz1,ddz1)+c15.*TX.*imfilter(vz1,Xdx1)+c11.*TX.*imfilter(vx1,ddx1)+c15.*TZ.*imfilter(vx1,Xdz1);
    taozz1=taozz0+c33.*TZ.*imfilter(vz1,ddz1)+c35.*TX.*imfilter(vz1,Xdx1)+c13.*TX.*imfilter(vx1,ddx1)+c35.*TZ.*imfilter(vx1,Xdz1);
    taoxz1=taoxz0+c55.*TZ.*imfilter(vx1,ddz0)+c15.*TX.*imfilter(vx1,Xdx0)+c55.*TX.*imfilter(vz1,ddx0)+c35.*TZ.*imfilter(vz1,Xdz0);
    
    taoxx0=taoxx1;taozz0=taozz1;taoxz0=taoxz1;
    vz0=vz1;vx0=vx1;
    
    if mod(t,0.5)==0
        figure();
        imagesc(vz1);
    end
end
