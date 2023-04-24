function s = wavelet(f,dt,T)
%WAVELET make wavlet.
% By zhaoqingwei
% Chengdu University of Technology (CDUT), 2021-2025
    if nargin < 3
        f=20;dt=0.001;T=4;
    end
    fm=f;
    dtt=dt;
    t=1/f;
    tt=-t:dtt:t;
    s=zeros(1,round(T/dtt));
    wave=(1-2*(pi*fm*tt).^2).*exp(-(pi*fm*tt).^2);
    
    s(1:length(wave))=wave;
end

