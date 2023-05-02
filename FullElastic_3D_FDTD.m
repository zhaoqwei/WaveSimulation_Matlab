clc,clear
% 3D - full elastic wave simulation using finite difference.
% By zhaoqingwei
% Chengdu University of Technology (CDUT), 2021-2025
FM=20;
DT=0.001;dt=DT;
T=2;
nx=200;
ny=200;
nz=200;
DH=10;dx=DH;dz=DH;dy=DH;


c11=zeros(nz,nx,ny);c12=c11;c13=c11;c14=c11;c15=c11;c16=c11;
c22=c11;c23=c11;c24=c11;c25=c11;c26=c11;
c33=c11;c34=c11;c35=c11;c36=c11;
c44=c11;c45=c11;c46=c11;c55=c11;c56=c11;c66=c11;
s=wavelet(FM,DT,T);z0=round(nz/2);x0=round(nx/2);y0=round(ny/2);
vp=ones(nz,nx,ny)*2000;vs=vp/2;rou=ones(nz,nx,ny)*2;
c11=rou.*vp.*vp;c22=c11;c33=c11;c55=rou.*vs.*vs;c44=c55;c66=c55;c13=rou.*(vp.*vp-vs.*vs*2);c12=c13;c23=c13;%各向同性
% Epsilon=0.1*ones(nz,nx,ny);Delta=0.1*ones(nz,nx,ny);Gamma=0.2*ones(nz,nx,ny);
% c33=vp.*vp.*rou;c11=2*c33.*Epsilon+c33;c22=c11;
% c55=vs.*vs.*rou;c44=c55;c66=2*c44.*Gamma+c44;
% c13=rou.*sqrt(((1+2*Delta).*vp.*vp-vs.*vs).*(vp.*vp-vs.*vs))-rou.*vs.*vs;c23=c13;c12=c11-2*c66;%VTI
%     [c11 c12 c13 c14 c15 c16;
%     c12 c22 c23 c24 c25 c26;
%     c13 c23 c33 c34 c35 c36;
%     c14 c24 c34 c44 c45 c46;
%     c15 c25 c35 c45 c55 c56;
%     c16 c26 c36 c46 c56 c66;];%full stiffness


vx0=zeros(nz,nx,ny);vx1=vx0;vx2=vx0;
vz0=vx0;vz1=vx0;vz2=vx0;
vy0=vx0;vy1=vx0;vy2=vx0;

nn=2;
dxx=FDcoeffDxx(nn);dzz=dxx';dyy=permute(dzz,[3 2 1]);
dxz=FDcoeffDxz(nn);dxy=permute(dxz,[3 2 1]);dyz=permute(dxz,[1 3 2]);
TX=dt/dx;TZ=dt/dz;TY=dt/dy;TTXX=TX*TX;
for	t=DT:DT:T
    disp(t);
    k=round(t/DT);
	VxDxx=imfilter(vx1,dxx);VxDxy=imfilter(vx1,dxy);VxDxz=imfilter(vx1,dxz);
    VxDyy=imfilter(vx1,dyy);VxDyz=imfilter(vx1,dyz);VxDzz=imfilter(vx1,dzz);
	VyDxx=imfilter(vy1,dxx);VyDxy=imfilter(vy1,dxy);VyDxz=imfilter(vy1,dxz);
    VyDyy=imfilter(vy1,dyy);VyDyz=imfilter(vy1,dyz);VyDzz=imfilter(vy1,dzz);
    VzDxx=imfilter(vz1,dxx);VzDxy=imfilter(vz1,dxy);VzDxz=imfilter(vz1,dxz);
    VzDyy=imfilter(vz1,dyy);VzDyz=imfilter(vz1,dyz);VzDzz=imfilter(vz1,dzz);
	
	vx2=2*vx1-vx0+TTXX*(c11.*VxDxx+2*c16.*VxDxy+2*c15.*VxDxz+2*c56.*VxDyz+c66.*VxDyy+c55.*VxDzz+...
        c16.*VyDxx+(c66+c12).*VyDxy+(c56+c14).*VyDxz+(c25+c46).*VyDyz+c26.*VyDyy+c45.*VyDzz+...
        c15.*VzDxx+(c56+c14).*VzDxy+(c55+c13).*VzDxz+(c45+c36).*VzDyz+c46.*VzDyy+c35.*VzDzz ...
        )./rou;
    vy2=2*vy1-vy0+TTXX*(c16.*VxDxx+(c12+c66).*VxDxy+(c14+c56).*VxDxz+(c46+c25).*VxDyz+c26.*VxDyy+c45.*VxDzz+...
        c66.*VyDxx+(c26+c26).*VyDxy+(c46+c46).*VyDxz+(c24+c24).*VyDyz+c22.*VyDyy+c44.*VyDzz+...
        c56.*VzDxx+(c25+c46).*VzDxy+(c45+c36).*VzDxz+(c44+c23).*VzDyz+c24.*VzDyy+c34.*VzDzz ...
        )./rou;
    vz2=2*vz1-vz0+TTXX*(c15.*VxDxx+(c14+c56).*VxDxy+(c13+c55).*VxDxz+(c36+c45).*VxDyz+c46.*VxDyy+c35.*VxDzz+...
        c56.*VyDxx+(c46+c24).*VyDxy+(c36+c45).*VyDxz+(c23+c44).*VyDyz+c24.*VyDyy+c34.*VyDzz+...
        c55.*VzDxx+(c45+c45).*VzDxy+(c35+c35).*VzDxz+(c34+c34).*VzDyz+c44.*VzDyy+c33.*VzDzz ...
        )./rou;
        
    vz2(z0,x0,y0)=vz2(z0,x0,y0)+s(k);
	vx0=vx1;vx1=vx2;
    vy0=vy1;vy1=vx2;
	vz0=vz1;vz1=vz2;
    
    
    if mod(t,0.5)==0
        figure();
        imagesc(vz1(:,:,y0)); 
    end
end
