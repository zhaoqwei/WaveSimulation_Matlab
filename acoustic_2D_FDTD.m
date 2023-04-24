clc,clear
% Acoustic wave simulation using finite difference.
% By zhaoqingwei
% Chengdu University of Technology (CDUT), 2021-2025
FM=20;
DT=0.001;
T=4;
nx=600;
nz=600;
DH=10;

s=wavelet(FM,DT,T);
u1=zeros(nz,nx);u2=u1;u3=u1;
v=ones(nz,nx)*2000;
z0=round(nz/2);x0=round(nx/2);
f=zeros(nx,nz);f(z0,x0)=1;

nn=3;
dxd=FDcoeffDxx(nn);dd=zeros(2*nn-1,2*nn-1);dd(nn,:)=dd(nn,:)+dxd;dd(:,nn)=dd(:,nn)+dxd';
for	t=DT:DT:T
    k=round(t/DT);
	const1=v.*v*DT*DT/DH/DH;
	UU=imfilter(u2,dd);
	u3=2*u2-u1+const1.*UU+s(k)*f;

	u1=u2;
	u2=u3;
    
    if mod(t,1)==0
        figure();
        imagesc(u3);
    end
end
