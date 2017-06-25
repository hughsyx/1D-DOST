s = [1/5.6,1/6.0,1/6.4,1/6.8];
f = [0:4096-1]/4096/0.05;
f = f-10;

d05 = distance(latlon(5,1),latlon(5,2),latlon(4:15,1),latlon(4:15,2),almanac('earth','ellipsoid'));
d15 = distance(latlon(15,1),latlon(15,2),latlon(4:15,1),latlon(4:15,2),almanac('earth','ellipsoid'));
for i=1:length(s)
    for j=1:length(d05)
        delta05(i,j) = sqrt(d05(j)^2+50^2)*s(i);
    end
    delta05(i,:) = delta05(i,:) - delta05(i,2);
end
delta_Ind05 = fix(delta05/0.05);

for i=1:length(s)
    for j=1:length(d15)
        delta15(i,j) = sqrt(d15(j)^2+50^2)*s(i);
    end
    delta15(i,:) = delta15(i,:) - delta15(i,end);
end
delta_Ind15 = fix(delta15/0.05);



cc_05_15_filt = zeros(1,7999);
cc_05_15 = zeros(1,7999);
tic
for iwin = 1:96
    test = squeeze(C1_A01_filt(18000+12000:18000+16000-1,iwin,4:15));
    S05 = test(:,2); S15 = test(:,12);
    if sum(S05) ==0 | sum(S15) ==0
        continue
    end
    
    nsta = size(test,2);
    test_pad = zeros(4096,nsta);
    test_pad(49:end-48,:) = test;
    test_dost = dost(test_pad);
    rearr_dost_coeff = zeros(4096,4096,nsta);
    test_phi = zeros(4096,4096,nsta);
    for i=1:nsta
        if sum(test_dost(:,i)) ~= 0
            rearr_dost_coeff(:,:,i) = rearrange_dost(test_dost(:,i));
            test_phi(:,:,i) = rearr_dost_coeff(:,:,i)./(abs(rearr_dost_coeff(:,:,i))+eps);
        end
    end
    
    c1 = zeros(size(rearr_dost_coeff(:,:,1)));
    for i=max(delta_Ind05(1,:)+1):4096-max(delta_Ind05(1,:)+1)
        for j = 2049:2049+300
            tmp = 0;
            tmp = tmp+test_phi(j,i-delta_Ind05(1,1),1)*exp(1i*2*pi*f(j)*delta05(1,1));
            tmp = tmp+test_phi(j,i+delta_Ind05(1,2),2)*exp(-1i*2*pi*f(j)*delta05(1,2));
            tmp = tmp+test_phi(j,i+delta_Ind05(1,3),3)*exp(-1i*2*pi*f(j)*delta05(1,3));
            tmp = tmp+test_phi(j,i+delta_Ind05(1,4),4)*exp(-1i*2*pi*f(j)*delta05(1,4));
            tmp = tmp+test_phi(j,i+delta_Ind05(1,5),5)*exp(-1i*2*pi*f(j)*delta05(1,5));
            tmp = tmp+test_phi(j,i+delta_Ind05(1,6),6)*exp(-1i*2*pi*f(j)*delta05(1,6));
            tmp = tmp+test_phi(j,i+delta_Ind05(1,7),7)*exp(-1i*2*pi*f(j)*delta05(1,7));
            tmp = tmp+test_phi(j,i+delta_Ind05(1,8),8)*exp(-1i*2*pi*f(j)*delta05(1,8));
            tmp = tmp+test_phi(j,i+delta_Ind05(1,9),9)*exp(-1i*2*pi*f(j)*delta05(1,9));
            tmp = tmp+test_phi(j,i+delta_Ind05(1,10),10)*exp(-1i*2*pi*f(j)*delta05(1,10));
            tmp = tmp+test_phi(j,i+delta_Ind05(1,11),11)*exp(-1i*2*pi*f(j)*delta05(1,11));
            tmp = tmp+test_phi(j,i+delta_Ind05(1,12),12)*exp(-1i*2*pi*f(j)*delta05(1,12));
            c1(j,i) = tmp/12;
        end
    end
    c1 = abs(c1);
    
    c2 = zeros(size(c1));
    for i=max(delta_Ind05(2,:)+1):4096-max(delta_Ind05(2,:)+1)
        for j = 2049:2049+300
            tmp = 0;
            tmp = tmp+test_phi(j,i-delta_Ind05(2,1),1)*exp(1i*2*pi*f(j)*delta05(2,1));
            tmp = tmp+test_phi(j,i+delta_Ind05(2,2),2)*exp(-1i*2*pi*f(j)*delta05(2,2));
            tmp = tmp+test_phi(j,i+delta_Ind05(2,3),3)*exp(-1i*2*pi*f(j)*delta05(2,3));
            tmp = tmp+test_phi(j,i+delta_Ind05(2,4),4)*exp(-1i*2*pi*f(j)*delta05(2,4));
            tmp = tmp+test_phi(j,i+delta_Ind05(2,5),5)*exp(-1i*2*pi*f(j)*delta05(2,5));
            tmp = tmp+test_phi(j,i+delta_Ind05(2,6),6)*exp(-1i*2*pi*f(j)*delta05(2,6));
            tmp = tmp+test_phi(j,i+delta_Ind05(2,7),7)*exp(-1i*2*pi*f(j)*delta05(2,7));
            tmp = tmp+test_phi(j,i+delta_Ind05(2,8),8)*exp(-1i*2*pi*f(j)*delta05(2,8));
            tmp = tmp+test_phi(j,i+delta_Ind05(2,9),9)*exp(-1i*2*pi*f(j)*delta05(2,9));
            tmp = tmp+test_phi(j,i+delta_Ind05(2,10),10)*exp(-1i*2*pi*f(j)*delta05(2,10));
            tmp = tmp+test_phi(j,i+delta_Ind05(2,11),11)*exp(-1i*2*pi*f(j)*delta05(2,11));
            tmp = tmp+test_phi(j,i+delta_Ind05(2,12),12)*exp(-1i*2*pi*f(j)*delta05(2,12));
        c2(j,i) = max(abs(tmp)/12,c1(j,i));
        end
    end
    
    c3 = zeros(size(c1));
    for i=max(delta_Ind05(3,:)+1):4096-max(delta_Ind05(3,:)+1)
        for j = 2049:2049+300
            tmp = 0;
            tmp = tmp+test_phi(j,i-delta_Ind05(3,1),1)*exp(1i*2*pi*f(j)*delta05(3,1));
            tmp = tmp+test_phi(j,i+delta_Ind05(3,2),2)*exp(-1i*2*pi*f(j)*delta05(3,2));
            tmp = tmp+test_phi(j,i+delta_Ind05(3,3),3)*exp(-1i*2*pi*f(j)*delta05(3,3));
            tmp = tmp+test_phi(j,i+delta_Ind05(3,4),4)*exp(-1i*2*pi*f(j)*delta05(3,4));
            tmp = tmp+test_phi(j,i+delta_Ind05(3,5),5)*exp(-1i*2*pi*f(j)*delta05(3,5));
            tmp = tmp+test_phi(j,i+delta_Ind05(3,6),6)*exp(-1i*2*pi*f(j)*delta05(3,6));
            tmp = tmp+test_phi(j,i+delta_Ind05(3,7),7)*exp(-1i*2*pi*f(j)*delta05(3,7));
            tmp = tmp+test_phi(j,i+delta_Ind05(3,8),8)*exp(-1i*2*pi*f(j)*delta05(3,8));
            tmp = tmp+test_phi(j,i+delta_Ind05(3,9),9)*exp(-1i*2*pi*f(j)*delta05(3,9));
            tmp = tmp+test_phi(j,i+delta_Ind05(3,10),10)*exp(-1i*2*pi*f(j)*delta05(3,10));
            tmp = tmp+test_phi(j,i+delta_Ind05(3,11),11)*exp(-1i*2*pi*f(j)*delta05(3,11));
            tmp = tmp+test_phi(j,i+delta_Ind05(3,12),12)*exp(-1i*2*pi*f(j)*delta05(3,12));
        c3(j,i) = max(abs(tmp)/12,c2(j,i));
        end
    end
    
    c4 = zeros(size(c1));
    for i=max(delta_Ind05(4,:)+1):4096-max(delta_Ind05(4,:)+1)
        for j = 2049:2049+300
            tmp = 0;
            tmp = tmp+test_phi(j,i-delta_Ind05(4,1),1)*exp(1i*2*pi*f(j)*delta05(4,1));
            tmp = tmp+test_phi(j,i+delta_Ind05(4,2),2)*exp(-1i*2*pi*f(j)*delta05(4,2));
            tmp = tmp+test_phi(j,i+delta_Ind05(4,3),3)*exp(-1i*2*pi*f(j)*delta05(4,3));
            tmp = tmp+test_phi(j,i+delta_Ind05(4,4),4)*exp(-1i*2*pi*f(j)*delta05(4,4));
            tmp = tmp+test_phi(j,i+delta_Ind05(4,5),5)*exp(-1i*2*pi*f(j)*delta05(4,5));
            tmp = tmp+test_phi(j,i+delta_Ind05(4,6),6)*exp(-1i*2*pi*f(j)*delta05(4,6));
            tmp = tmp+test_phi(j,i+delta_Ind05(4,7),7)*exp(-1i*2*pi*f(j)*delta05(4,7));
            tmp = tmp+test_phi(j,i+delta_Ind05(4,8),8)*exp(-1i*2*pi*f(j)*delta05(4,8));
            tmp = tmp+test_phi(j,i+delta_Ind05(4,9),9)*exp(-1i*2*pi*f(j)*delta05(4,9));
            tmp = tmp+test_phi(j,i+delta_Ind05(4,10),10)*exp(-1i*2*pi*f(j)*delta05(4,10));
            tmp = tmp+test_phi(j,i+delta_Ind05(4,11),11)*exp(-1i*2*pi*f(j)*delta05(4,11));
            tmp = tmp+test_phi(j,i+delta_Ind05(4,12),12)*exp(-1i*2*pi*f(j)*delta05(4,12));
        c4(j,i) = max(abs(tmp)/12,c3(j,i));
        end
    end 
    
    c4 = c4/(max(max(abs(c4)))+eps);
    c4(2:2048,:) = flipud(c4(2050:4096,:));
    c4_2 = c4.^2;
    
    S05_filt_rearr_dost_coeff = rearr_dost_coeff(:,:,2).*c4_2;
    S05_filt_dost_coeff = rearrange_dost_back(S05_filt_rearr_dost_coeff);
    S05_filt_pad = real(idost(S05_filt_dost_coeff));
    S05_filt = S05_filt_pad(49:end-48);
    
    clear c1 c2 c3 c4 c4_2 tmp
    
    
    
    
    c1 = zeros(size(rearr_dost_coeff(:,:,1)));
    for i=max(delta_Ind15(1,:)+1):4096-max(delta_Ind15(1,:)+1)
        for j = 2049:2049+300
            tmp = 0;
            tmp = tmp+test_phi(j,i-delta_Ind15(1,1),1)*exp(1i*2*pi*f(j)*delta15(1,1));
            tmp = tmp+test_phi(j,i-delta_Ind15(1,2),2)*exp(1i*2*pi*f(j)*delta15(1,2));
            tmp = tmp+test_phi(j,i-delta_Ind15(1,3),3)*exp(1i*2*pi*f(j)*delta15(1,3));
            tmp = tmp+test_phi(j,i-delta_Ind15(1,4),4)*exp(1i*2*pi*f(j)*delta15(1,4));
            tmp = tmp+test_phi(j,i-delta_Ind15(1,5),5)*exp(1i*2*pi*f(j)*delta15(1,5));
            tmp = tmp+test_phi(j,i-delta_Ind15(1,6),6)*exp(1i*2*pi*f(j)*delta15(1,6));
            tmp = tmp+test_phi(j,i-delta_Ind15(1,7),7)*exp(1i*2*pi*f(j)*delta15(1,7));
            tmp = tmp+test_phi(j,i-delta_Ind15(1,8),8)*exp(1i*2*pi*f(j)*delta15(1,8));
            tmp = tmp+test_phi(j,i-delta_Ind15(1,9),9)*exp(1i*2*pi*f(j)*delta15(1,9));
            tmp = tmp+test_phi(j,i-delta_Ind15(1,10),10)*exp(1i*2*pi*f(j)*delta15(1,10));
            tmp = tmp+test_phi(j,i-delta_Ind15(1,11),11)*exp(1i*2*pi*f(j)*delta15(1,11));
            tmp = tmp+test_phi(j,i+delta_Ind15(1,12),12)*exp(-1i*2*pi*f(j)*delta15(1,12));
            c1(j,i) = tmp/12;
        end
    end
    c1 = abs(c1);
    
    c2 = zeros(size(c1));
    for i=max(delta_Ind15(2,:)+1):4096-max(delta_Ind15(2,:)+1)
        for j = 2049:2049+300
            tmp = 0;
            tmp = tmp+test_phi(j,i-delta_Ind15(2,1),1)*exp(1i*2*pi*f(j)*delta15(2,1));
            tmp = tmp+test_phi(j,i-delta_Ind15(2,2),2)*exp(1i*2*pi*f(j)*delta15(2,2));
            tmp = tmp+test_phi(j,i-delta_Ind15(2,3),3)*exp(1i*2*pi*f(j)*delta15(2,3));
            tmp = tmp+test_phi(j,i-delta_Ind15(2,4),4)*exp(1i*2*pi*f(j)*delta15(2,4));
            tmp = tmp+test_phi(j,i-delta_Ind15(2,5),5)*exp(1i*2*pi*f(j)*delta15(2,5));
            tmp = tmp+test_phi(j,i-delta_Ind15(2,6),6)*exp(1i*2*pi*f(j)*delta15(2,6));
            tmp = tmp+test_phi(j,i-delta_Ind15(2,7),7)*exp(1i*2*pi*f(j)*delta15(2,7));
            tmp = tmp+test_phi(j,i-delta_Ind15(2,8),8)*exp(1i*2*pi*f(j)*delta15(2,8));
            tmp = tmp+test_phi(j,i-delta_Ind15(2,9),9)*exp(1i*2*pi*f(j)*delta15(2,9));
            tmp = tmp+test_phi(j,i-delta_Ind15(2,10),10)*exp(1i*2*pi*f(j)*delta15(2,10));
            tmp = tmp+test_phi(j,i-delta_Ind15(2,11),11)*exp(1i*2*pi*f(j)*delta15(2,11));
            tmp = tmp+test_phi(j,i+delta_Ind15(2,12),12)*exp(-1i*2*pi*f(j)*delta15(2,12));
        c2(j,i) = max(abs(tmp)/12,c1(j,i));
        end
    end
    
    c3 = zeros(size(c1));
    for i=max(delta_Ind15(3,:)+1):4096-max(delta_Ind15(3,:)+1)
        for j = 2049:2049+300
            tmp = 0;
            tmp = tmp+test_phi(j,i-delta_Ind15(3,1),1)*exp(1i*2*pi*f(j)*delta15(3,1));
            tmp = tmp+test_phi(j,i-delta_Ind15(3,2),2)*exp(1i*2*pi*f(j)*delta15(3,2));
            tmp = tmp+test_phi(j,i-delta_Ind15(3,3),3)*exp(1i*2*pi*f(j)*delta15(3,3));
            tmp = tmp+test_phi(j,i-delta_Ind15(3,4),4)*exp(1i*2*pi*f(j)*delta15(3,4));
            tmp = tmp+test_phi(j,i-delta_Ind15(3,5),5)*exp(1i*2*pi*f(j)*delta15(3,5));
            tmp = tmp+test_phi(j,i-delta_Ind15(3,6),6)*exp(1i*2*pi*f(j)*delta15(3,6));
            tmp = tmp+test_phi(j,i-delta_Ind15(3,7),7)*exp(1i*2*pi*f(j)*delta15(3,7));
            tmp = tmp+test_phi(j,i-delta_Ind15(3,8),8)*exp(1i*2*pi*f(j)*delta15(3,8));
            tmp = tmp+test_phi(j,i-delta_Ind15(3,9),9)*exp(1i*2*pi*f(j)*delta15(3,9));
            tmp = tmp+test_phi(j,i-delta_Ind15(3,10),10)*exp(1i*2*pi*f(j)*delta15(3,10));
            tmp = tmp+test_phi(j,i-delta_Ind15(3,11),11)*exp(1i*2*pi*f(j)*delta15(3,11));
            tmp = tmp+test_phi(j,i+delta_Ind15(3,12),12)*exp(-1i*2*pi*f(j)*delta15(3,12));
        c3(j,i) = max(abs(tmp)/12,c2(j,i));
        end
    end
    
    c4 = zeros(size(c1));
    for i=max(delta_Ind15(4,:)+1):4096-max(delta_Ind15(4,:)+1)
        for j = 2049:2049+300
            tmp = 0;
            tmp = tmp+test_phi(j,i-delta_Ind15(4,1),1)*exp(1i*2*pi*f(j)*delta15(4,1));
            tmp = tmp+test_phi(j,i-delta_Ind15(4,2),2)*exp(1i*2*pi*f(j)*delta15(4,2));
            tmp = tmp+test_phi(j,i-delta_Ind15(4,3),3)*exp(1i*2*pi*f(j)*delta15(4,3));
            tmp = tmp+test_phi(j,i-delta_Ind15(4,4),4)*exp(1i*2*pi*f(j)*delta15(4,4));
            tmp = tmp+test_phi(j,i-delta_Ind15(4,5),5)*exp(1i*2*pi*f(j)*delta15(4,5));
            tmp = tmp+test_phi(j,i-delta_Ind15(4,6),6)*exp(1i*2*pi*f(j)*delta15(4,6));
            tmp = tmp+test_phi(j,i-delta_Ind15(4,7),7)*exp(1i*2*pi*f(j)*delta15(4,7));
            tmp = tmp+test_phi(j,i-delta_Ind15(4,8),8)*exp(1i*2*pi*f(j)*delta15(4,8));
            tmp = tmp+test_phi(j,i-delta_Ind15(4,9),9)*exp(1i*2*pi*f(j)*delta15(4,9));
            tmp = tmp+test_phi(j,i-delta_Ind15(4,10),10)*exp(1i*2*pi*f(j)*delta15(4,10));
            tmp = tmp+test_phi(j,i-delta_Ind15(4,11),11)*exp(1i*2*pi*f(j)*delta15(4,11));
            tmp = tmp+test_phi(j,i+delta_Ind15(4,12),12)*exp(-1i*2*pi*f(j)*delta15(4,12));
        c4(j,i) = max(abs(tmp)/12,c3(j,i));
        end
    end 
    
    c4 = c4/(max(max(abs(c4)))+eps);
    c4(2:2048,:) = flipud(c4(2050:4096,:));
    c4_2 = c4.^2;
    
    S15_filt_rearr_dost_coeff = rearr_dost_coeff(:,:,end).*c4_2;
    S15_filt_dost_coeff = rearrange_dost_back(S15_filt_rearr_dost_coeff);
    S15_filt_pad = real(idost(S15_filt_dost_coeff));
    S15_filt = S15_filt_pad(49:end-48);
    
    clear c1 c2 c3 c4 c4_2 tmp test*
    cc_05_15_filt = cc_05_15_filt + CANIR(S15_filt,S05_filt,10)';
    cc_05_15= cc_05_15 + CANIR(S15,S05,10)';
    %plot(filtfilt(B,A,cc_05_15_filt),'k')
    iwin
end
toc
    
    