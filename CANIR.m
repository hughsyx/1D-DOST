function ANIRW=CANIR(vb,va,npt)
% calculate the ANIR of each window pair
% va is the virtual source
% npt is the number used for running window average

len1=length(va);len2=length(vb);
len=len1+len2-1;
fva=fft(va,len);
fvb=fft(vb,len);

eps=1e-9;
%running window average
% fva2=zeros(len+2*npt,1);
% if size(fva,1)~=1
%     fva2(1:npt)=fliplr(fva(1:npt)')';
%     fva2(npt+len+1:end)=fliplr(fva(len-npt+1:len)')';
% else
%     fva2(1:npt)=fliplr(fva(1:npt));
%     fva2(npt+len+1:end)=fliplr(fva(len-npt+1:len));
% end
% fva2(npt+1:npt+len)=fva;
% fva2=abs(fva2);

fva2=RWA2(fva,npt); fvb2=RWA2(fvb,npt);
p1=ones(size(fva))*eps; p2=ones(size(fvb))*eps;
for i=1:len
    p1(i)=sum(fva2(i:i+npt-1))/npt;
    p2(i)=sum(fvb2(i:i+npt-1))/npt;
end

ANIRWF=fvb.*conj(fva)./(p1.*p2);% p2?
ANIRW=real(fftshift(ifft(ANIRWF)));

%function v2 = RWA2(v,npt)
%    len = length(v);
%    v2 = zeros(len+npt-1,1);
%    if size(v,1)~=1
%        v2(1:npt)=fliplr(v(1:npt)')';
%        v2(npt+len+1:end)=fliplr(v(len-npt+1:len)')';
%    else
%        v2(1:npt)=fliplr(v(1:npt));
%        v2(npt+len+1:end)=fliplr(v(len-npt+1:len));
%    end
%    v2(npt+1:npt+len)=v;
%    v2=abs(v2);
%end
function v2 = RWA2(v,npt)
    temp = zeros(1,npt-1);
    len = length(v);
    v2 = zeros(len+npt-1,1);
    v2(1:len) = v;
    temp=v(end-npt+1:end-1);
    v2(len+1:end)=fliplr(temp);
    v2=abs(v2);
end
end
