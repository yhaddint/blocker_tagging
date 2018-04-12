%
%  07/29/2014
%  We assume ADC works in 500MSamples/s and PN sequences have bandiwdth of
%  5MHz. We are dealing with blockers with BW less than 500KHz because
%  otherwise 1FFT 


clear;clc;clf;close all
warning off

checkPN=0;

% determined whether to apply estimation of power and detection
plot_stat=0;
runtimes=50;

upsam=20;

figindex=99;

N=20;
M=30;

%  tune power of each blockers
blocker_num=10;

% 3 blockers in -40 to -30dBm, 2 blockers in -50 to -40dBm
blocker_pow_pool(:,1)=1/100;


%  Define a reasonable length of signal. It means each simulation we can
%  get 60 statistics
stat_num=60;
L=M*N*stat_num;


symbols=L*2/upsam;


%  Blocker envelop
for ii=1:blocker_num*2
    clear data
    hmod = modem.pskmod('M', 4, 'InputType', 'integer');
    hdesign  = fdesign.pulseshaping(upsam,'Square Root Raised Cosine');
    hpulse = design(hdesign);
    data = randi(4,symbols,1)-1;
    data = modulate(hmod, data);
    data = upsample(data,upsam);
    sig_bb(:,ii) = conv(data,hpulse.Numerator);
end

%  Normalization of Power
blocker(:,1)=real(sig_bb(end-L+1:end,1))./sqrt(mean(real(sig_bb(end-L+1:end,1)).^2))*blocker_pow_pool(1);


cal0_temp=PNgenerator_v4(110);
cal0=zeros(N,blocker_num);
for ii=1:blocker_num
cal0(:,ii)=cal0_temp(ii:N-1+ii);
end
% Equivalently lower rate of PN signal
cal=LowerRate_v2(cal0,L);

%% Check PN Cross-Correlation
if checkPN
    corrMat=zeros(blocker_num,blocker_num);
    for ii=1:blocker_num
        for jj=ii:blocker_num
            corrMat(ii,jj)=sum(cal0(:,ii).*cal0(:,jj));
        end
    end
    corrMat
end

%% down sample
Ndown=N;
%% standard
%  We define Blocker with certain power, let's say 10, as standard.
%  It can be regarded as training signal to determine TH 

standard_raw=0.5*blocker(:,1).^2+blocker(:,1).*cal(:,1);
standard_MNsorted=tagging_v1(standard_raw,cal(:,1),Ndown,M);
TH=standard_MNsorted(0.1*length(standard_MNsorted))*2.5;
 
%% parallel PN correlation and averaging
achive=zeros(blocker_num,runtimes);


for runindex=1:runtimes
%  permutation of calibration sequences
cal=cal(:,randperm(blocker_num));

%  Target:2-10, Trivial:0.2-0.9, Critical FA:0.1, 0.05, 0.01
blocker_pow_pool(:,2:4)=(rand(1,3)*8+2)/100;
blocker_pow_pool(:,5:7)=(rand(1,3)*0.7+0.2)/100;
blocker_pow_pool(:,8)=0.1/100;
blocker_pow_pool(:,9)=0.05/100;
blocker_pow_pool(:,10)=0.01/100;

%  QPSK Blocker generator with power scaled
for ii=1:blocker_num
    tempindex=randi(20);
    blocker(:,ii)=real(sig_bb(end-L+1:end,tempindex))./sqrt(mean(real(sig_bb(end-L+1:end,tempindex)).^2))*blocker_pow_pool(ii);
end

%  calculate statistics (including intermediate)
sig_raw=zeros(M*Ndown*stat_num,1);
for ii=1:blocker_num
    sig_raw=sig_raw+0.5*blocker(:,ii).^2+blocker(:,ii).*cal(:,ii);
end

for ii=1:blocker_num
    sig_MNave(:,ii)=tagging_v2(sig_raw,cal(:,ii),Ndown,M);
end

%  initialization of tagged_table
tagged_table_1=zeros(blocker_num,stat_num);
tagged_table_2=zeros(blocker_num,stat_num);
for kk=1:stat_num
    %  Find those blockers who have been tagged
    for ii=1:blocker_num
        if sig_MNave(kk,ii)>TH
            tagged_table_1(ii,kk)=1;
            blocker((kk-1)*M*Ndown+1:kk*M*Ndown,ii)=zeros(M*Ndown,1);
        end
    end
end

    %  Stage-2 Tagging Starting from here
    
    %  Since the PN is 10dB less, equivalent blocker power will be 10 times
    %  as before
    blocker=blocker*10;
    
    sig_raw=zeros(M*Ndown*stat_num,1);
    for ii=1:blocker_num
        sig_raw=sig_raw+0.5*blocker(:,ii).^2+blocker(:,ii).*cal(:,ii);
    end

    for ii=1:blocker_num
        sig_MNave(:,ii)=tagging_v2(sig_raw,cal(:,ii),Ndown,M);
    end
    
    
    %  Find those blockers who have been tagged (tagging stage 2)
    for kk=1:stat_num
        for ii=1:blocker_num
            if tagged_table_1(ii,kk)==0 && sig_MNave(kk,ii)>TH
                tagged_table_2(ii,kk)=1;
            end
        end
    end               



% watch stat histogram and TH_standard
if plot_stat
    for ii=1:blocker_num
    figure(figindex)
    subplot(2,5,ii)
    stem(sig_MNave(:,ii));hold on
    title(['Power =, ',num2str(100*blocker_pow_pool(ii))])
    end
    figindex=figindex+1;
end

    for ii=1:blocker_num
        achive_total(ii,runindex)=sum(tagged_table_1(ii,:)+tagged_table_2(ii,:))/stat_num;
        achive_1(ii,runindex)=sum(tagged_table_1(ii,:))/stat_num;
        achive_2(ii,runindex)=sum(tagged_table_2(ii,:))/stat_num;
    end
    
%  Leakage will cause FA, we determined a certain TH and calculate per of
%  FA, as well as missing
%  Obviously, here the TH is finely tuned.

end



