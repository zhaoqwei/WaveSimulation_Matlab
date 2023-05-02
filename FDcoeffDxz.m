function fdcoeff=FDcoeffDxz(n)
%FDcoeffDxz£º  Finite difference coefficient
% By zhaoqingwei
% Chengdu University of Technology (CDUT), 2021-2025
% fdcoeff=zeros(3,3);
dxP5=FDcoeffDx(n);
dx=imfilter(dxP5,[0.5 0.5],'full');
fdcoeff=imfilter(dx,dx','full');
% dx=[-1 1];
% fdcoeff=imfilter(dx,dx','full');