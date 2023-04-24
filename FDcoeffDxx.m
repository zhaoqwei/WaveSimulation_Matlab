function fdcoeff=FDcoeffDxx(n)
%FDcoeffDxx£º 2st Finite difference coefficient
% By zhaoqingwei
% Chengdu University of Technology (CDUT), 2021-2025
a=zeros(n,n);b=zeros(n,1);fdcoeff=zeros(1,2*n-1);
for i=0:n-1
    for j=0:n-1
        a(i+1,j+1)=2*j^(2*i)/factorial(2*i);
    end
end
b(2)=1;
x=a\b;
fdcoeff(1:n)=fliplr(x');
fdcoeff(n:2*n-1)=fdcoeff(n:2*n-1)+x';