% 03/25/2014
%
% This code is to generate QSKP signal by comm_toolbox and verify 
% certian CLT approximation valid or not

% 03/27/2014
% PN sequences with characterize polynomial have been put into use

clear;clc;clf;close all
warning off
symbols=3e4;
upsam=8;
for jj=1:10

for i=1:3
clear data
hmod = modem.pskmod('M', 4, 'InputType', 'integer');
hdesign  = fdesign.pulseshaping(upsam,'Square Root Raised Cosine');
hpulse = design(hdesign);
data = randi(4,symbols,1)-1;
data = modulate(hmod, data);
data = upsample(data,upsam);
data = conv(data,hpulse.Numerator);

sig_bb(:,i)=data;
sig_bb(:,i)=sig_bb(:,i)/(sqrt(mean(abs(sig_bb(:,i)).^2)));
end





%% distribution of signal_r^2*cal
N=4e4;
M=30;
clear ii i
sig_use=real(sig_bb(1:N,1)).^2;
sig_power=mean(real(sig_bb(:,1)).^2);
cal=PNgenerator_v1(N);
sig_pn=sig_use.*cal;

% % simulate "uniformly random starting point"
% for ii=1:N-M
%     rslt1(ii)=mean(sig_pn(ii:ii+M-1));
% end

% simulate "sequencially starting point"
for ii=1:fix(N/M/2)
    rslt2(ii)=mean(sig_pn((ii-1)*M+1:ii*M));
end


% for ii=1:2000
%     rslt3(ii)=mean(sig_pn(randi(5000,50,1)));
% end
%%
% figure
% hist(rslt1,20)
% title('sequencially pick')
% figure
% hist(rslt2,20)
% title('sequencially pick')
% figure
% hist(rslt3,20)
% title('randomly pick')

%%
% clear ii
% for ii=1:fix(N/M)
%     rslt4(ii)=mean(cal((ii-1)*M+1:ii*M).^3);
% end

% %% sequencially pick illustration
% figure
% plot(real(sig_bb(1:500,1)),'c','linewidth',2);
% hold
% plot(100:150,real(sig_bb(100:150,1)),'kx');
% grid
% 
% 
% %% randomly pick illustration
% index=randi(500,1,50);
% figure
% plot(real(sig_bb(1:500,1)),'c','linewidth',2);
% hold
% plot(sort(index),real(sig_bb(sort(index))),'kx');
% grid

[H,gaussfit(jj)]=chi2gof(rslt2);
end
plot(50,gaussfit,'x')