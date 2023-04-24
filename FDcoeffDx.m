function fdcoeff=FDcoeffDx(n)
%FDcoeffDx£º 1st Finite difference coefficient
% By zhaoqingwei
% Chengdu University of Technology (CDUT), 2021-2025
B=zeros(n,1);B(1)=1;
Aa=zeros(1,n);
for k=1:n
    Aa(k)=2*k-1;
end
AA=ones(n,1)*Aa;
A=AA.^(AA');
x=A\B;
fdcoeff=[-fliplr(x') x']; 