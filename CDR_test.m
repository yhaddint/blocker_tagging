% Ratio of constructively added components versus distructively added ones
clear;clc;clf;close all
warning off

checkPN=0;

% determined whether to apply estimation of power and detection
plot_stat=0;

runtimes=2e4;

Npool=4:20;
M=40;
blocker_num=40;

%  Define a reasonable length of signal. It means each simulation we can
%  get 60 statistics
stat_num=60;
L=M*20*stat_num;
%% blocker
for ii=1:blocker_num
upsam=10;
symbols=fix(L*2/upsam);

clearvars data temp_data
hmod = modem.pskmod('M', 4, 'InputType', 'integer');
hdesign  = fdesign.pulseshaping(upsam,'Square Root Raised Cosine');
hpulse = design(hdesign);
data = randi(4,symbols,1)-1;
data = modulate(hmod, data);
data = upsample(data,upsam);
temp_data = conv(data,hpulse.Numerator);
sig_bb = temp_data(end-L+1:end);

blocker(:,ii)=real(sig_bb)./sqrt(mean(real(sig_bb)).^2)*0.1;
end
%% cal
cal0_temp=PNgenerator_v1(110);
cal0=zeros(20,10);
for ii=1:10
    cal0(:,ii)=cal0_temp(ii:19+ii);
end

cal = LowerRate(cal0, 1, L);
rrpool = randi(L-30,1,runtimes);
cal1pool = randi(10,1,runtimes);
cal2pool = cal1pool+randi(9,1,runtimes);
cal2pool = cal2pool-(cal2pool>10)*10;
bpool = randi(blocker_num,1,runtimes);
%%
for nn=1:length(Npool)
    for runindex=1:runtimes
   
    N=Npool(nn);
        
    cal1index=cal1pool(runindex);
    cal2index=cal2pool(runindex);
    rr=rrpool(runindex);
    sig_noPN=blocker(rr:rr+N-1,bpool(runindex));
    sig_PN=0.5.*sig_noPN.^2+sig_noPN.*cal(1:N,cal1index);

    destructPower(runindex,nn)=(sum(sig_PN.*cal(1:N,cal2index)));
    constructPower(runindex,nn)=(sum(sig_noPN));
    
    end
end

for nn=1:length(Npool)
    destructPower(:,nn)=destructPower(:,nn)-2*destructPower(:,nn).*(destructPower(:,nn)<0);
    constructPower(:,nn)=constructPower(:,nn)-2*constructPower(:,nn).*(constructPower(:,nn)<0);
    temp=10*log10(constructPower(:,nn)./destructPower(:,nn));
    %temp=temp.*(temp<(40))+40*(temp>=(40));
    lb=fix(length(temp)*0.1);
    hb=fix(length(temp)*0.9);
    temp=sort(temp);
    CDR(nn)=mean(temp(lb:hb));
    
end
%%
figure
oddi=1:2:length(Npool);
eveni=2:2:length(Npool);

plot(Npool,CDR,'rx-','linewidth',2);hold on
%plot(Npool(eveni),CDR(eveni),'bx-','linewidth',2);hold on
grid on
xlabel('Parameter N')
ylabel('Constructive Distructive Ratio')
%%
figure(99)
nn=15;
plot(sort(10*log10(constructPower(1:1e4,nn)./destructPower(1:1e4,nn))),'r');hold on
plot(sort(10*log10(constructPower(1+1e4:end,nn)./destructPower(1+1e4:end,nn))),'b');
xlabel('Sample Number')
ylabel('Ratio')
