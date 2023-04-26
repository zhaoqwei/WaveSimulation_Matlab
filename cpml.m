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
        ddz=zeros(z,x);
        I=pmln;
        for iz=1:z
            if iz<=I
                ddz(iz,:)=(I-iz+1)*(I-iz+1)/I/I*ones(1,x);
            end
            if iz>=z-I
                ddz(iz,:)=(-z+I+iz)*(-z+I+iz)/I/I*ones(1,x);
            end
        end
        AA=-log(R)*Vmax*3/2/(pmln*dz);
        ddz=ddz*AA;
        bbz=exp(-ddz*dt);
        aaz=bbz-1;
        bbx = permute(bbz,[2 1 3]);
        aax = permute(aaz,[2 1 3]);
        bby=[];
        aay=[];
    end
    if nD==3
        z=Size(1);x=Size(2);y=Size(3);
        ddz=zeros(z,x,y);
        I=pmln;
        for iz=1:z
            if iz<=I
                ddz(iz,:,:)=(I-iz+1)*(I-iz+1)/I/I*ones(x,y);
            end
            if iz>=z-I
                ddz(iz,:,:)=(-z+I+iz)*(-z+I+iz)/I/I*ones(x,y);
            end
        end
        AA=-log(R)*Vmax*3/2/(pmln*dz);
        ddz=ddz*AA;
        bbz=exp(-ddz*dt);aaz=bbz-1;
        bbx = permute(bbz,[2 1 3]);aax = permute(aaz,[2 1 3]);
        bby = permute(bbz,[3 2 1]);aay = permute(aaz,[3 2 1]);
    end
end

