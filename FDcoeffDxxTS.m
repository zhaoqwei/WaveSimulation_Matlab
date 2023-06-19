function a = FDcoeffDxxTS(M,Ang,r)
%FDcoeffDxxTS£º 2st Finite difference coefficient Time-space
% By zhaoqingwei
% Chengdu University of Technology (CDUT), 2021-2025
a=zeros(2*M+1,1);
A=zeros(M,M);
B=zeros(M,1);

for j=1:M
     for m=1:M
            A(j,m)=m^(2*j)*((cos(Ang))^(2*j)+ (sin(Ang))^(2*j));
     
     end
    B(j)=r^(2*j-2);
end

size(A);
size(B);

c=A\B;

a(M+1)=-2*sum(c);
a(1:M)=flip(c);
a(M+2:end)=c;

end

