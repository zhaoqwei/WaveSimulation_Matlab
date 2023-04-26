clc,clear
% 3D - Acoustic wave simulation using finite difference.
% By zhaoqingwei
% Chengdu University of Technology (CDUT), 2021-2025
FM=20;
DT=0.001;
T=2;
nx=200;
ny=200;
nz=200;
DH=10;

s=wavelet(FM,DT,T);
u1=zeros(nz,nx,ny);u2=u1;u3=u1;
v=ones(nz,nx,ny)*2000;
z0=round(nz/2);x0=round(nx/2);y0=round(ny/2);
f=zeros(nz,nx,ny);f(z0,x0,x0)=1;

nn=3;
dxd=FDcoeffDxx(nn);ddz=zeros(2*nn-1,2*nn-1,2*nn-1);ddz(:,nn,nn)=dxd;ddx = permute(ddz,[2 1 3]);ddy = permute(ddz,[3 2 1]);dd=ddz+ddx+ddy;
for	t=DT:DT:T
    k=round(t/DT);
	const1=v.*v*DT*DT/DH/DH;
	UU=imfilter(u2,dd);
	u3=2*u2-u1+const1.*UU+s(k)*f;

	u1=u2;
	u2=u3;
    
    if mod(t,0.5)==0
        figure();
        imagesc(u3(:,:,y0)); 
    end
end
