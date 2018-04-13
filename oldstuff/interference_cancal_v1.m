%  03/31/2016
% let s_k be the strongest blocker
%  use digital bandpass filters to find s_k * c_i e^{j(k-i)\omega t}
% it works well when only one block exists

clear;clc;%clf;close all


rand('seed',3)
blk_num = 4;
M = 1e2;
N = 63;
SampleNumberAve = 32;
stat_num = 2e2;
P = M*N;
L = M*SampleNumberAve*stat_num;


other_pow = 1e-4;
sig_pow0 = ones(blk_num,1)*other_pow;
sig_pow1 = ones(blk_num,1)*other_pow;

sig_pow1([1,4]) = [0.1/2,1];
sig_pow0([1,4]) = [0,1];

for ii = 1:blk_num
%     sig_r(:,ii) = kron(randi(2,M*stat_num,1)*2-3,ones(SampleNumberAve,1));
%     sig_i(:,ii) = kron(randi(2,M*stat_num,1)*2-3,ones(SampleNumberAve,1));
%     sig(:,ii) = (sig_r(:,ii)+sig_i(:,ii))/sqrt(2);
    upsam = SampleNumberAve;
    symbols = fix(L*2/upsam);
    
    clearvars data temp_data
    hmod = modem.qammod('M', 4, 'InputType', 'integer');
    hdesign  = fdesign.pulseshaping(upsam,'Square Root Raised Cosine');
    hpulse = design(hdesign);
    data = randi(4,symbols,1)-1;
    data = modulate(hmod, data);
    data = upsample(data,upsam);
    temp_data = conv(data,hpulse.Numerator);
%     unnormalized = real(temp_data(end-L+1:end));
%     sig_r(:,ii) = unnormalized./sqrt(sum(unnormalized.^2)/L);
%     unnormalized = imag(temp_data(end-L+1:end));
%     sig_i(:,ii) = unnormalized./sqrt(sum(unnormalized.^2)/L);
%     sig_i(:,ii) = zeros(size(sig_r(:,ii)));
%     sig(:,ii) = (sig_r(:,ii)+sig_i(:,ii))/sqrt(2);
    sig(:,ii) = temp_data(end-L+1:end)./sqrt(temp_data(end-L+1:end)'*temp_data(end-L+1:end)/L);

end
%% calibration sequences
cal0=PNgenerator_v5(N,N,1);
CAL_noperm=LowerRate_v2(cal0,P);
CAL=zeros(P,N);
for mm=1:M
    CAL((mm-1)*N+1:mm*N,:)=CAL_noperm((mm-1)*N+1:mm*N,randperm(N));
end
%% tagging

CAL=zeros(P,N);
for mm=1:M
    CAL((mm-1)*N+1:mm*N,:)=CAL_noperm((mm-1)*N+1:mm*N,randperm(N));
end

%CAL=CAL(:,randperm(N));
start_point = randi(L-P-1,blk_num,M);
for ii=1:blk_num
    sig_cache(:,ii) = sig(start_point(ii,mm):start_point(ii,mm)+P-1,randi(blk_num));
end

sig_mix0_part1 = zeros(P,1);
sig_mix1_part1 = zeros(P,1);

sig_mix0_part2 = zeros(P,1);
sig_mix1_part2 = zeros(P,1);

sig_mix0 = zeros(P,1);
sig_mix1 = zeros(P,1);
PhaseShift = exp(1j*rand(1,blk_num));

for ii=1:blk_num
    sig_mix0_part1 = sig_mix0_part1+sig_cache(:,ii).*CAL(:,ii)*PhaseShift(ii)*sqrt(sig_pow0(ii));
    sig_mix1_part1 = sig_mix1_part1+sig_cache(:,ii).*CAL(:,ii)*PhaseShift(ii)*sqrt(sig_pow1(ii));
end

for ii=1:blk_num
    for jj=1:blk_num
        if ii ~= jj
        carrier = exp(1j*((ii-jj)/32*2*pi*(1:P)'))*PhaseShift(ii);
        sig_mix0_part2 = sig_mix0_part2 + sqrt(sig_pow0(ii))*CAL(:,jj).*...
            real(sig_cache(:,ii).*carrier);
        sig_mix1_part2 = sig_mix1_part2 + sqrt(sig_pow1(ii))*CAL(:,jj).*...
            real(sig_cache(:,ii).*carrier);
        end
    end
end
sig_mix0 = sig_mix0_part1+sig_mix0_part2;
sig_mix1 = sig_mix1_part1+sig_mix1_part2;
%%
% noise = (randn(P,1)+1j*randn(P,1))*sqrt(10);
% temp = (sig_mix1_part1(1:4950)+noise(1:4950)).*CAL(1:4950,1);
% temp2 = abs(fft(temp));
% lengthtemp = length(temp2);
% figure;
% plot([temp2(lengthtemp/2:lengthtemp);temp2(1:lengthtemp/2-1)]);hold on
% %
% noise_filtered = filter(bhi,1,noise.*CAL(1:N*M,2));
% noise_new = noise(1:4950)-noise_filtered(51:5000).*CAL(1:4950,2);
% 
% %
% temp = (sig_mix1_part1(1:4950)+noise_new(1:4950)).*CAL(1:4950,1);
% temp2 = abs(fft(temp));
% lengthtemp = length(temp2);
% plot([temp2(lengthtemp/2:lengthtemp);temp2(1:lengthtemp/2-1)],'r');hold on
%%
% for ii=1:blk_num   
%     temp = filter(bhi,1,sig_mix1.*CAL(1:N*M,ii));
% %     temp = filter(bhi,1,sig_mix1_part2.*CAL(1:N*M,ii));
%     cancle_term(:,ii) = temp(51:5050).*CAL(1:5000,ii);
% end
% sig_nocancal = sig_mix1_part1(1:5e3)+sig_mix1_part2(1:5e3);
% sig_cancal = sig_nocancal-mean(cancle_term(:,2:end),2);
% 
% temp1 = abs(fft(sig_nocancal.*CAL(1:5e3,1)));
% temp2 = abs(fft(sig_cancal.*CAL(1:5e3,1)));
% figure(10)
% plot([temp1(2500:end);temp1(1:2500-1)]);hold on
% figure(20)
% plot([temp2(2500:end);temp2(1:2500-1)]);hold on
% 
% for ii=1:blk_num 
%     temp = filter(bhi,1,sig_mix0.*CAL(1:N*M,ii));
% %     temp = filter(bhi,1,sig_mix0_part2.*CAL(1:N*M,ii));
%     cancle_term(:,ii) = temp(51:5050).*CAL(1:5000,ii);
% end
% sig_nocancal = sig_mix0_part1(1:5e3)+sig_mix0_part2(1:5e3);
% sig_cancal = sig_nocancal-mean(cancle_term(:,2:end),2);
% 
% temp1 = abs(fft(sig_nocancal.*CAL(1:5e3,1)));
% temp2 = abs(fft(sig_cancal.*CAL(1:5e3,1)));
% figure(10)
% plot([temp1(2500:end);temp1(1:2500-1)]);hold on
% figure(20)
% plot([temp2(2500:end);temp2(1:2500-1)]);hold on
% 
% 
% %%
% ii=1;
% 
% ydata0 = abs(fft(sig_mix0(1:5e3).*CAL(1:5e3,1)));
% ydata1 = abs(fft(sig_mix1(1:5e3).*CAL(1:5e3,1)));
% l_ydata = length(ydata1);
% xdata = linspace(-1,1,l_ydata);
% figure;
% plot(xdata,[ydata1(l_ydata/2:end);ydata1(1:l_ydata/2-1)],'k');hold on
% plot(xdata,[ydata0(l_ydata/2:end);ydata0(1:l_ydata/2-1)],'g');hold on
% DLPF = abs(xdata)<1/16;
% plot(xdata,800*DLPF,'k--');
% ylim([0,1000])
% xlabel('Frequency (\pi)')
% ylabel('relative amplitude')
% 
% 


%% figure

% figure
% plot(abs(fft(sig_cache(:,1))));hold on
% xlabel('Number of Frames (M)')
% ylabel('Identification Probability')
% %%
% 
% ydata = abs(fft(sig_mix0(1:N*M).*CAL(1:N*M,ii)));
% l_ydata = length(ydata)
% figure;
% plot(ydata)
% figure;
% plot([ydata(l_ydata/2:end);ydata(1:l_ydata/2)]); 
%%
bhi(1,:) = fir1(100,[5/32-1/64,7/32+1/64]);
% freqz(bhi(1,:),1)
%
bhi(2,:) = fir1(100,[3/32-1/64,5/32+1/64]);
% freqz(bhi(2,:),1)
%
bhi(3,:) = fir1(100,[1/32-1/64,3/32+1/64]);
% freqz(bhi(3,:),1)
%
bhi(4,:) = fir1(100,1/32+1/64,'low');
% freqz(bhi(4,:),1)
%
for ii=1:4
    yyy = filter(bhi(ii,:),1,sig_mix1.*CAL(1:M*N,ii));
    comp_term(:,ii) = yyy(51:5050).*CAL(1:5e3,ii);
end
sig_nocancal = sig_mix1(1:5e3);
sig_cancal = sig_nocancal-sum(comp_term(:,1:4),2);
temp1 = abs(fft(sig_nocancal.*CAL(1:5e3,1)));
temp2 = abs(fft(sig_cancal.*CAL(1:5e3,1)));
xdata = linspace(-1,1,5e3);
figure(10)
subplot(211)
plot(xdata,[temp1(2500:end);temp1(1:2500-1)]);hold on
figure(10)
subplot(212)
plot(xdata,[temp2(2500:end);temp2(1:2500-1)]);hold on



for ii=1:4
    yyy = filter(bhi(ii,:),1,sig_mix0.*CAL(1:M*N,ii));
    comp_term(:,ii) = yyy(51:5050).*CAL(1:5e3,ii);
end
sig_nocancal = sig_mix0(1:5e3);
sig_cancal = sig_nocancal-sum(comp_term(:,1:4),2);
temp1 = abs(fft(sig_nocancal.*CAL(1:5e3,1)));
temp2 = abs(fft(sig_cancal.*CAL(1:5e3,1)));
figure(10)
subplot(211)
plot(xdata,[temp1(2500:end);temp1(1:2500-1)]);hold on
figure(10)
subplot(212)
plot(xdata,[temp2(2500:end);temp2(1:2500-1)]);hold on

%%
xxx = sig_mix0(1:5e3);
yyy = filter(bhi(4,:),1,sig_mix0.*CAL(1:M*N,4));
zzz = yyy(51:5050).*CAL(1:5e3,4);

%
% 
figure;
subplot(211)
plot(abs(fft(xxx.*CAL(1:5e3,1))));hold on
subplot(212)
plot(abs(fft(zzz.*CAL(1:5e3,1))))
