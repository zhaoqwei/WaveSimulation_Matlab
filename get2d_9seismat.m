function seismat=get2d_9seismat(v,w,dx,dz,pmlx,pmlz)
% nine points rotate grid different
%
[z,x]=size(v);
a=0.5461;c=0.6248;d=0.09381;  
%1996year jo a=0.5461,c=0.6248,d=0.09381;

%% 
dt=1;% Fake, no-uses 
[aaz,bbz,aax,bbx] = cpml(1e-10,[z,x],pmlx,max(v(:)),dt,dz);
ddz=-log(bbz)/dt;ddx=-log(bbx)/dt;

ex=1-1i*ddx/w;
ez=1-1i*ddz/w;
%%
V_ma=reshape(v.',z*x,1);
EX_ma=reshape(ex.',z*x,1);
EZ_ma=reshape(ez.',z*x,1);

%%
med=2*(a./(EZ_ma.*EZ_ma*dz*dz))+2*(a./(EX_ma.*EX_ma*dx*dx))+2*((1-a)./(EZ_ma.*EZ_ma*(dx*dx+dz*dz)))...
    +2*((1-a)./(EX_ma.*EX_ma*(dx*dx+dz*dz)))-c*w*w./(V_ma.*V_ma);
ax=2*(1-a)./(2*EX_ma.*EX_ma*(dx*dx+dz*dz))+(1-c-4*d)*w*w./(4*V_ma.*V_ma);
az=2*(1-a)./(2*EZ_ma.*EZ_ma*(dx*dx+dz*dz))+(1-c-4*d)*w*w./(4*V_ma.*V_ma);
bx=2*a./(2*(EX_ma.*EX_ma*dx*dx))+d*w*w./(V_ma.*V_ma);
bz=2*a./(2*(EZ_ma.*EZ_ma*dz*dz))+d*w*w./(V_ma.*V_ma);
seismat=-spdiags(med,0,x*z,x*z)+spdiags([0;bx(1:end-1)],1,x*z,x*z)+spdiags([bx(2:end);0],-1,x*z,x*z)+spdiags([zeros(x,1);bz(1:end-x)],x,x*z,x*z)+...
    spdiags([bz(1+x:end);zeros(x,1)],-x,x*z,x*z)+spdiags([zeros(x,1);ax(2:end-x);0],(x-1),x*z,x*z)+spdiags([0;ax(1+x:(end-1));zeros(x,1)],-(x-1),x*z,x*z)+...
    spdiags([zeros(x+1,1);az(1:(end-1-x))],(x+1),x*z,x*z)+spdiags([az(2+x:(end));zeros(x+1,1)],-(x+1),x*z,x*z);
for k=(x+1):x:z*x
    seismat(k-1,k)=0;seismat(k,k-1)=0;
end
for k=(2*x+1):x:z*x
    seismat(k-x-1,k)=0;seismat(k,k-x-1)=0;seismat(k-x,k-1)=0;seismat(k-1,k-x)=0;
end
