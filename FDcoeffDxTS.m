function a = FDcoeffDxTS(M,Ang,r)
%FDcoeffDxTS£º 1st Finite difference coefficient Time-space
% By zhaoqingwei
% Chengdu University of Technology (CDUT), 2021-2025
a=zeros(2*M,1);
A=zeros(M,M);
B=zeros(M,1);

for j=1:M
     for m=1:M
         
            A(j,m)=(2*m-1)^(2*j-1)*((cos(Ang))^(2*j-1)+ (sin(Ang))^(2*j-1));
     
     end
    B(j)=r^(2*j-2);
end

size(A);
size(B);

c=A\B;


a(1:M)=-flip(c);
a(M+1:end)=c;

end

