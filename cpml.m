function [ aaz,bbz,aax,bbx,aay,bby ] = cpml( R,Size,pmln,Vmax,dt,dz)
%CPML: make Conv Perfectly Matched Layer
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
        I=pmln;
        for iz=1:z
            if iz<=I
                ddz(iz,:)=(I-iz+1)*(I-iz+1)/I/I*ones(1,x);
            end
            if iz>=z-I
                ddz(iz,:)=(-z+I+iz)*(-z+I+iz)/I/I*ones(1,x);
            end
        end
        for ix=1:x
            if ix<=I
                ddx(ix,:)=(I-ix+1)*(I-ix+1)/I/I*ones(1,z);
            end
            if ix>=x-I
                ddx(ix,:)=(-x+I+ix)*(-x+I+ix)/I/I*ones(1,z);
            end
        end
        
        AA=-log(R)*Vmax*3/2/(pmln*dz);
        ddz=ddz*AA;
        bbz=exp(-ddz*dt);
        aaz=bbz-1;
        ddx=ddx*AA;
        bbx1=exp(-ddx*dt);
        aax1=bbx1-1;
        
        bbx = permute(bbx1,[2 1 3]);
        aax = permute(aax1,[2 1 3]);
        bby=[];
        aay=[];
    end
    if nD==3
        z=Size(1);x=Size(2);y=Size(3);
        ddz=zeros(z,x,y);ddx=zeros(x,z,y);ddy=zeros(y,x,z);
        I=pmln;
        for iz=1:z
            if iz<=I
                ddz(iz,:,:)=(I-iz+1)*(I-iz+1)/I/I*ones(x,y);
            end
            if iz>=z-I
                ddz(iz,:,:)=(-z+I+iz)*(-z+I+iz)/I/I*ones(x,y);
            end
        end
        for ix=1:x
            if ix<=I
                ddx(ix,:,:)=(I-ix+1)*(I-ix+1)/I/I*ones(z,y);
            end
            if ix>=x-I
                ddx(ix,:,:)=(-x+I+ix)*(-x+I+ix)/I/I*ones(z,y);
            end
        end
        for iy=1:y
            if iy<=I
                ddy(iy,:,:)=(I-iy+1)*(I-iy+1)/I/I*ones(x,z);
            end
            if iy>=y-I
                ddy(iy,:,:)=(-y+I+iy)*(-y+I+iy)/I/I*ones(x,z);
            end
        end
        AA=-log(R)*Vmax*3/2/(pmln*dz);
        ddz=ddz*AA;ddx=ddx*AA;ddy=ddy*AA;
        bbz=exp(-ddz*dt);aaz=bbz-1;
        bbx1=exp(-ddx*dt);aax1=bbx1-1;
        bby1=exp(-ddy*dt);aay1=bby1-1;
        bbx = permute(bbx1,[2 1 3]);aax = permute(aax1,[2 1 3]);
        bby = permute(bby1,[3 2 1]);aay = permute(aay1,[3 2 1]);
    end
end

