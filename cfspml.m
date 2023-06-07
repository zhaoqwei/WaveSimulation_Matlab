function [bbz,ccz,dddz,bbx,ccx,dddx ] = cfspml( R,Size,pmln,Vmax,dt,dz,Gmax,Amax)
%CFS-PML: make CFS Perfectly Matched Layer
% By zhaoqingwei
% Chengdu University of Technology (CDUT), 2021-2025
    if nargin < 6
        R=1e-10;
        pmln=20;
        Size=[600 600];
        Vmax=2000;
        dt=0.001;
        dz=10;
    end
    nD=length(Size);
    if nD<2
        error('dimension size is error ');
    end
    
    if nD==2
        z=Size(1);x=Size(2);
        ddz=zeros(z,x);ddx=zeros(x,z);
        ggz=ones(z,x);ggx=ones(x,z);
        alphz=zeros(z,x);alphx=zeros(x,z);
        I=pmln;
        for iz=1:z
            if iz<=I
                ddz(iz,:)=(I-iz+1)*(I-iz+1)/I/I*ones(1,x);
                ggz(iz,:)=1+(Gmax-1)*(((I-iz+1)/I)^2)*ones(1,x);
                alphz(iz,:)=Amax*(iz)/I*ones(1,x);
            end
            if iz>z-I
                ddz(iz,:)=(-z+I+iz)*(-z+I+iz)/I/I*ones(1,x);
                ggz(iz,:)=1+(Gmax-1)*(((-z+I+iz)/I)^2)*ones(1,x);
                alphz(iz,:)=Amax*(z-iz+1)/I*ones(1,x);
            end
        end
        for ix=1:x
            if ix<=I
                ddx(ix,:)=(I-ix+1)*(I-ix+1)/I/I*ones(1,z);
                ggx(ix,:)=1+(Gmax-1)*(((I-ix+1)/I)^2)*ones(1,z);
                alphx(ix,:)=Amax*(ix)/I*ones(1,z);
            end
            if ix>x-I
                ddx(ix,:)=(-x+I+ix)*(-x+I+ix)/I/I*ones(1,z);
                ggx(ix,:)=1+(Gmax-1)*(((-x+I+ix)/I)^2)*ones(1,z);
                alphx(ix,:)=Amax*(x-ix+1)/I*ones(1,z);
            end
        end
        AA=-log(R)*Vmax*3/2/(pmln*dz);
        ddz=ddz*AA;
        bbz=exp(-(alphz+ddz./ggz)*dt);
        ccz=ddz./(ddz.*ggz+ggz.*ggz.*alphz).*(bbz-1);
        ccz(isnan(ccz))=0;
        dddz=(1-ggz)./ggz;


        AA=-log(R)*Vmax*3/2/(pmln*dz);
        ddx=ddx*AA;
        bbx1=exp(-(alphx+ddx./ggx)*dt);
        ccx1=ddx./(ddx.*ggx+ggx.*ggx.*alphx).*(bbx1-1);
        ccx1(isnan(ccx1))=0;
        dddx1=(1-ggx)./ggx;
        
        dddx = permute(dddx1,[2 1 3]);
        ccx = permute(ccx1,[2 1 3]);
        bbx = permute(bbx1,[2 1 3]);
    end
end

