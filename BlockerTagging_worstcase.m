 
%  08/06/2014
%  Worst case: all blockers existing with full power, except one, whose
%  power will be swept to check PD

%  3 stages to deal with dynamic range problem


clear;clc;clf;close all
warning off

checkPN=0;

% determined whether to apply estimation of power and detection
plot_stat=0;
runtimes=50;

usepowdetection=0;

figindex=99;

N=31;
M=60;
P=M*N;

%  tune power of each blockers
blocker_num=6;

% 3 blockers in -40 to -30dBm, 2 blockers in -50 to -40dBm
blocker_pow_pool(:,1)=1/100;

%  Define a reasonable length of signal. It means each simulation we can
%  get 60 statistics
stat_num=60;
L=M*N*stat_num;
blockercan_num=blocker_num*2;
%  Blocker envelop
for ii=1:blockercan_num
    
    upsam=randi([20 200]);
    symbols=fix(L*2/upsam);
    
    clearvars data temp_data
    hmod = modem.pskmod('M', 4, 'InputType', 'integer');
    hdesign  = fdesign.pulseshaping(upsam,'Square Root Raised Cosine');
    hpulse = design(hdesign);
    data = randi(4,symbols,1)-1;
    data = modulate(hmod, data);
    data = upsample(data,upsam);
    temp_data = conv(data,hpulse.Numerator);
    sig_bb(:,ii) = temp_data(end-L+1:end);
end

%  Normalization of Power
blocker(:,1)=real(sig_bb(:,1))./sqrt(mean(real(sig_bb(:,1)).^2))*blocker_pow_pool(1);


cal0_temp=PNgenerator_v1(N);
cal0=zeros(N,blocker_num);
for ii=1:blocker_num
    cal0(:,ii)=[cal0_temp(ii:end);cal0_temp(1:ii-1)];
end
% Equivalently lower rate of PN signal
cal=LowerRate(cal0, 1, L);

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

%% standard
%  We define Blocker with certain power, let's say 10, as standard.
%  It can be regarded as training signal to determine TH 

standard_raw=0.5*blocker(:,1).^2+blocker(:,1).*cal(:,1);

% 
if usepowdetection
    for ii=1:stat_num
    standard_MN(ii) = correlation_and_pow(standard_raw((ii-1)*P+1:ii*P).*cal((ii-1)*P+1:ii*P,1),N,M);
    end
else
    standard_MN = correlation_and_abs(standard_raw,cal(:,1),N,M);
end
standard_MNsorted=sort(standard_MN);
TH=standard_MNsorted(0.1*length(standard_MNsorted));
 
%% parallel PN correlation and averaging
THpool=[1 1.5 2 2.5 3 3.5]*TH;
for THindex=1:length(THpool)
for runindex=1:runtimes
%  permutation of calibration sequences
cal=cal(:,randperm(blocker_num));

%  Target:2-10, Trivial:0.2-0.9, Critical FA:0.1, 0.05, 0.01
blocker_pow_pool(:,1:blocker_num-1)=100/1000;
blocker_pow_pool(:,blocker_num)=0;

%  QPSK Blocker generator with power scaled
for ii=1:blocker_num
    tempindex=randi(blockercan_num);
    blocker(:,ii)=real(sig_bb(:,tempindex))./sqrt(mean(real(sig_bb(:,tempindex)).^2))*blocker_pow_pool(ii);
end

%  calculate statistics (including intermediate)
sig_raw=zeros(M*N*stat_num,1);
for ii=1:blocker_num
    sig_raw=sig_raw+0.5*blocker(:,ii).^2+blocker(:,ii).*cal(:,ii);
end

for ii=1:blocker_num
    if usepowdetection
        for jj=1:stat_num
        sig_MNave(jj,ii)=correlation_and_pow(sig_raw((jj-1)*P+1:jj*P).*cal((jj-1)*P+1:jj*P,ii),N,M);
        end
    else
        sig_MNave(:,ii)=correlation_and_abs(sig_raw,cal(:,ii),N,M);
    end
end

%  initialization of tagged_table
tagged_table_1=zeros(blocker_num,stat_num);
tagged_table_2=zeros(blocker_num,stat_num);
tagged_table_3=zeros(blocker_num,stat_num);
for kk=1:stat_num
    %  Find those blockers who have been tagged
    for ii=1:blocker_num
        if sig_MNave(kk,ii)>TH
            tagged_table_1(ii,kk) = 1;
            blocker((kk-1)*M*N+1:kk*M*N,ii) = zeros(M*N,1);
        end
    end
end

    %  Stage-2 Tagging Starting from here
    
    %  Since the PN is 10dB less, equivalent blocker power will be 10 times
    %  as before
    blocker=blocker*10;
    
    sig_raw=zeros(M*N*stat_num,1);
    for ii=1:blocker_num
        sig_raw=sig_raw+0.5*blocker(:,ii).^2+blocker(:,ii).*cal(:,ii);
    end

    for ii=1:blocker_num
        if usepowdetection
            for jj=1:stat_num
            sig_MNave(jj,ii) = correlation_and_pow(sig_raw((jj-1)*P+1:jj*P).*cal((jj-1)*P+1:jj*P,ii),N,M);
            end
        else
            sig_MNave(:,ii) = correlation_and_abs(sig_raw,cal(:,ii),N,M);
        end
    end
    
    %  Find those blockers who have been tagged (tagging stage 2)
    for kk=1:stat_num
        for ii=1:blocker_num
            if tagged_table_1(ii,kk)==0 && sig_MNave(kk,ii)>TH
                tagged_table_2(ii,kk)=1;
                blocker((kk-1)*M*N+1:kk*M*N,ii)=zeros(M*N,1);

            end
        end
    end               

%  Stage-3 Tagging Starting from here
     
%     %  Since the PN is 10dB less, equivalent blocker power will be 10 times
%     %  as before
%     blocker=blocker*4.64;
%     
%     sig_raw=zeros(M*N*stat_num,1);
%     for ii=1:blocker_num
%         sig_raw=sig_raw+0.5*blocker(:,ii).^2+blocker(:,ii).*cal(:,ii);
%     end
% 
%     for ii=1:blocker_num
%         if usev3
%             for jj=1:stat_num
%             sig_MNave(jj,ii)=tagging_v3(sig_raw((jj-1)*P+1:jj*P).*cal((jj-1)*P+1:jj*P,ii),N,M);
%             end
%         else
%             sig_MNave(:,ii)=tagging_v2(sig_raw,cal(:,ii),N,M);
%         end
%     end
%     
%     
% %     Find those blockers who have been tagged (tagging stage 2)
%     for kk=1:stat_num
%         for ii=1:blocker_num
%             if tagged_table_1(ii,kk)==0 && tagged_table_2(ii,kk)==0 &&sig_MNave(kk,ii)>TH
%                 tagged_table_3(ii,kk)=1;
%             end
%         end
%     end
    
    
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
        achive_total(ii,runindex)=sum(tagged_table_1(ii,:)+tagged_table_2(ii,:)+tagged_table_3(ii,:))/stat_num;
        achive_1(ii,runindex)=sum(tagged_table_1(ii,:))/stat_num;
        achive_2(ii,runindex)=sum(tagged_table_2(ii,:))/stat_num;
        achive_3(ii,runindex)=sum(tagged_table_3(ii,:))/stat_num;
    end
    
%  Leakage will cause FA, we determined a certain TH and calculate per of
%  FA, as well as missing
%  Obviously, here the TH is finely tuned.

end
mean(achive_total')'
end


