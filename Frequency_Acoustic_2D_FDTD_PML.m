clc,clear
% Acoustic wave simulation using finite difference depend  Frequency. 
% with Perfectly Matched Layer;
% By zhaoqingwei
% Chengdu University of Technology (CDUT), 2021-2025
%%
%获得雷克子波
FM=20;
DT=0.001;dt=DT;
T=2;
nx=600;
nz=600;
DH=10;dx=DH;dz=DH;
pmln=20;

s=wavelet(FM,DT,T);
hs=hilbert(s);
fs=fft(hs);
nt=length(fs);

df=1/T;
w_pi=[0:1/nt:1/2 (-1/2+1/nt):1/nt:-1/nt]*2*pi;
w=w_pi/dt;
f=w/2/pi;

fmax=(70);fmin=(0);
fn=find(f>=fmax, 1 );%index in f for fmax.
f1=find(f>fmin, 1);%index in f for fmin.

%%
% init model
v=ones(nz,nx)*2000;
z0=round(nz/2);x0=round(nx/2);

fs_index=(x0-1)*nz+z0;
u_line=zeros(nz*nx,nt);
for i=f1:fn
    disp(i);
    seismat=get2d_9seismat(v,w(i),dx,dz,pmln,pmln);

    Fi=zeros(nz*nx,1);Fi(fs_index)=fs(i);
    u_1line=seismat\Fi;
%     [L,U]=lu(seismat);
%     u_1line=U\(L\Fi);
%     u_1line=lsqr(seismat,Fi,1e-10,300);
    u_line(:,i)=u_1line;
%     U=reshape(u_1line,nx,nz).';
%     figure;
%     imagesc(real(U));
end
u_l=ifft(u_line,[],2);
u=reshape(u_l,nz,nx,nt);
figure;imagesc(real(u(:,:,500)));


