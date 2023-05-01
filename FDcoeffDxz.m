function fdcoeff=FDcoeffDxz()
%FDcoeffDxz£º  Finite difference coefficient
% By zhaoqingwei
% Chengdu University of Technology (CDUT), 2021-2025
% fdcoeff=zeros(3,3);
dx=[-0.5 0 0.5];
fdcoeff=imfilter(dx,dx','full');