clc,clear
% 2D - full elastic wave simulation using finite difference.
% By zhaoqingwei
% Chengdu University of Technology (CDUT), 2021-2025
FM=20;
DT=0.001;dt=DT;
T=2;
nx=600;
nz=600;
DH=10;dx=DH;dz=DH;

s=wavelet(FM,DT,T);
vp=ones(nz,nx)*4000;vs=vp/2;rou=ones(nz,nx)*2;
lamda=rou.*(vp.*vp-vs.*vs*2);mu=rou.*vs.*vs;
c11=rou.*vp.*vp;c13=rou.*(vp.*vp-vs.*vs*2);c33=c11;c55=rou.*vs.*vs;
c15=zeros(nz,nx);c35=c15;%各向同性
% Epsilon=0.25*ones(600,600);Delta=0.1*ones(600,600);
% c33=vp.*vp.*rou;c55=vs.*vs.*rou;c11=2*c33.*Epsilon+c33;c13=rou.*sqrt(((1+2*Delta).*vp.*vp-vs.*vs).*(vp.*vp-vs.*vs))-rou.*vs.*vs;%VTI
% [c11 c13 c15;
%  c13 c33 c35;
%  c15 c35 c55;]%full stiffness


vx0=zeros(nz,nx);vx1=vx0;vx2=vx0;
vz0=vx0;vz1=vx0;vz2=vx0;
z0=round(nz/2);x0=round(nx/2);

nn=2;
dxx=FDcoeffDxx(nn);dxz=FDcoeffDxz(nn);dzz=dxx';
TX=dt/dx;TZ=dt/dz;TTXX=TX*TX;
for	t=DT:DT:T
    disp(t);
    k=round(t/DT);
	VxDxx=imfilter(vx1,dxx);
    VxDzz=imfilter(vx1,dzz);
    VxDxz=imfilter(vx1,dxz);
    VzDxx=imfilter(vz1,dxx);
    VzDzz=imfilter(vz1,dzz);
    VzDxz=imfilter(vz1,dxz);
	
	vx2=2*vx1-vx0+(TX*TX*c11.*VxDxx+2*TX*TZ*c15.*VxDxz+TZ*TZ*c55.*VxDzz +  TX*TX*c15.*VzDxx+TX*TZ*(c55+c13).*VzDxz+TZ*TZ*c35.*VzDzz)./rou;
    vz2=2*vz1-vz0+(TX*TX*c15.*VxDxx+TX*TZ*(c13+c55).*VxDxz+TZ*TZ*c35.*VxDzz +  TX*TX*c55.*VzDxx+2*TX*TZ*c35.*VzDxz+TZ*TZ*c33.*VzDzz)./rou;
    vz2(z0,x0)=vz2(z0,x0)+s(k);
	vx0=vx1;vx1=vx2;
	vz0=vz1;vz1=vz2;
    
    
    if mod(t,0.5)==0
        figure();
        imagesc(vz2);
    end
end
