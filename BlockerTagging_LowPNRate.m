%
%  07/18/2014
%  Using lower PN signal rates to watch what happened.


clear;clc;clf;close all
warning off

% determined whether to apply estimation of power and detection
do_est=0;
do_det=0;
plot_intermediate=1;

symbols=5e4;
upsam=500;


N=400;
M=10;

%  tune power of each blockers
blocker_num=3;
int_power_pool=[1 2 4 8 16 32]/1024';
%int_power_pool=0.5';

blocker_pow_pool=zeros(length(int_power_pool),blocker_num);
blocker_pow_pool(:,1)=1/1024;
blocker_pow_pool(:,2)=int_power_pool;
blocker_pow_pool(:,3)=int_power_pool;


mark=zeros(1,length(int_power_pool));   

L=8e5;

for ii=1:3
    clear data
    hmod = modem.pskmod('M', 4, 'InputType', 'integer');
    hdesign  = fdesign.pulseshaping(upsam,'Square Root Raised Cosine');
    hpulse = design(hdesign);
    data = randi(4,symbols,1)-1;
    data = modulate(hmod, data);
    data = upsample(data,upsam);
    sig_bb(:,ii) = conv(data,hpulse.Numerator);
end

%% standard
%  We define Blocker with certain power, let's say 10, as standard.
%  It can be regarded as training signal to determine TH
blocker(1,:)=real(sig_bb(1:L,1))./sqrt(mean(real(sig_bb(1:L,1)).^2))*blocker_pow_pool(1,1);
cal0=PNgenerator_v1(M*N*4);

    % Equivalently lower rate of PN signal
Dsampling=4;
    for ii=1:3
        % reuse PN sequences of length M*N
        cal(ii,:)=kron(ones(L/M/N/Dsampling,1),cal0(M*N*ii+1:M*N*(ii+1)));
        % Rate of PN is lower than Fs, which is equivalent to upsamling PN sequence
        cal_temp(ii,:)=kron(cal(ii,:),ones(1,Dsampling));
    end
cal=cal_temp;
clearvars cal_temp
    temp3=0.5*blocker(1,:).^2.*cal(1,:);
    temp4=blocker(1,:);
raw_standard=temp3+temp4;
for ii=1:fix(L/N/2)
    Nave_standard(ii)=abs(sum(raw_standard((ii-1)*N+1:ii*N)/N));
end
l_need=fix(length(Nave_standard)/M);    
MNave_standard=mean(reshape(Nave_standard(1:M*l_need),M,l_need));
sorted_standard=sort(MNave_standard);
TH_standard=sorted_standard(fix(length(sorted_standard)*0.2));   


%% h0,h1 sweep with different interference blocker power
for int_power_index=1:length(int_power_pool)
    
    %  QPSK Blocker generator with power scaled    
    for ii=2:3
        blocker(ii,:)=real(sig_bb(1:L,ii))./sqrt(mean(real(sig_bb(1:L,ii)).^2))*blocker_pow_pool(int_power_index,ii);
    end

    %  calculate statistics (including intermediate)
    temp1=0.5*(blocker(2,:).^2+blocker(3,:).^2).*cal(1,:);
    temp2=blocker(3,:).*cal(1,:).*cal(3,:)+blocker(2,:).*cal(1,:).*cal(2,:);
    
    raw_h0=temp1+temp2;
    raw_h1=raw_standard+raw_h0; 
    
    for ii=1:fix(L/N/2)
        temp1_sum(ii)=sum(temp1((ii-1)*N+1:ii*N)/N);
        temp2_sum(ii)=sum(temp2((ii-1)*N+1:ii*N)/N);
        temp3_sum(ii)=sum(temp3((ii-1)*N+1:ii*N)/N);
        temp4_sum(ii)=sum(temp4((ii-1)*N+1:ii*N)/N);
        Nave_h0(ii)=abs(sum(raw_h0((ii-1)*N+1:ii*N)/N));
        Nave_h1(ii)=abs(sum(raw_h1((ii-1)*N+1:ii*N)/N));
    end
    
    %% plot for intermediate variable
    if plot_intermediate
    figure
    subplot(2,2,1)
    stem(temp1_sum(1e2:1e2+100));
    subplot(2,2,2)
    stem(temp2_sum(1e2:1e2+100));
    subplot(2,2,3)
    stem(temp3_sum(1e2:1e2+100));
    subplot(2,2,4)
    stem(temp4_sum(1e2:1e2+100));
    end
%     clearvars temp1 temp2 temp3 temp4
%     clearvars temp1_sum temp2_sum temp3_sum temp4_sum
    %% M averaging part

    % statistics after N,M averaging
    MNave_h0(:,int_power_index)=mean(reshape(Nave_h0(1:M*l_need),M,l_need));
    MNave_h1(:,int_power_index)=mean(reshape(Nave_h1(1:M*l_need),M,l_need));

end


%% watch histogram and TH_standard
for ii=1:length(int_power_pool)
figure(99)
subplot(2,3,ii)
stem(sort(MNave_h0(:,ii)));hold on
title(['h0, ',num2str(int_power_pool(ii))])

figure(100)
subplot(2,3,ii)
stem(sort(MNave_h1(:,ii)));hold on
title(['h1, ',num2str(int_power_pool(ii))])
end

