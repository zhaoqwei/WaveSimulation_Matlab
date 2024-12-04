function [ M ] = bond( angel )
%bond£º  calculates the bond transformation of the elastic stiffness
% By zhaoqingwei
% Chengdu University of Technology (CDUT), 2021-2025
    a1=angel(1);a2=angel(2);a3=angel(3);
    A1=[cos(a1) sin(a1) 0;-sin(a1) cos(a1) 0;0 0 1;];
    A2=[1 0 0;0 cos(a2) sin(a2);0 -sin(a2) cos(a2);];
    A3=[cos(a3) 0 -sin(a3);0 1 0;sin(a3) 0 cos(a3);];
    A=A1*A2*A3;
    M1=[A(1,1)*A(1,1) A(1,2)*A(1,2) A(1,3)*A(1,3);
        A(2,1)*A(2,1) A(2,2)*A(2,2) A(2,3)*A(2,3);
        A(3,1)*A(3,1) A(3,2)*A(3,2) A(3,3)*A(3,3);];
    M2=[A(1,2)*A(1,3) A(1,3)*A(1,1) A(1,1)*A(1,2);
        A(2,2)*A(2,3) A(2,3)*A(2,1) A(2,1)*A(2,2);
        A(3,2)*A(3,3) A(3,3)*A(3,1) A(3,1)*A(3,2);];
    M3=[A(2,1)*A(3,1) A(2,2)*A(3,2) A(2,3)*A(3,3);
        A(3,1)*A(1,1) A(3,2)*A(1,2) A(3,3)*A(1,3);
        A(1,1)*A(2,1) A(1,2)*A(2,2) A(1,3)*A(2,3);];
    M4=[A(2,2)*A(3,3)+A(2,3)*A(3,2) A(2,1)*A(3,3)+A(2,3)*A(3,1) A(2,2)*A(3,1)+A(2,1)*A(3,2);
        A(1,2)*A(3,3)+A(1,3)*A(3,2) A(1,1)*A(3,3)+A(1,3)*A(3,1) A(1,1)*A(3,2)+A(1,2)*A(3,1);
        A(2,2)*A(1,3)+A(1,2)*A(2,3) A(1,1)*A(2,3)+A(1,3)*A(2,1) A(2,2)*A(1,1)+A(1,2)*A(2,1);];
    M = [M1, 2*M2; M3, M4];

end

