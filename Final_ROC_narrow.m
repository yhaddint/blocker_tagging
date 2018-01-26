%  10/02/2014
%  15 BLKs with Power uniformly from [-50,-30]dBm. BW from 100KHz to 1MHz
%  Pd=0.9, Pfa=0.1 get

clear;
clc;
warning off

worstcase=0;

% Independent repeating experiment.
runtimes=1e2;

%  parameter M and N setting. N is also length of CAL Sequences
N=15;
M=25;
downSampleTo=1;

SampleNumberAve=N*downSampleTo;

P=M*SampleNumberAve;

%  How many BLK we want to tag simutaneously
BLK_num=15;

%  How many active BLK existing. To get reasonable number of different
%  envelope, we generate 3x of them.

BLKcan_num=50;

% We assume simplified model hold (Strong CAL).
lowerThanPNdB=40;
lowerThanPN=10^(lowerThanPNdB/10);

%  Define length of BLK we want.
stat_num=20;
L=M*SampleNumberAve*stat_num;

%  Averaging Time = xT


%  Normalized BLK generation. 
%  CAL BW/ BLK BW = 100.
for ii=1:BLKcan_num
    upsam=randi([SampleNumberAve*2 SampleNumberAve*20]);
    symbols=max(fix(L*2/upsam),10);
    
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

%% calibration sequences
cal0=PNgenerator_v5(N,N,downSampleTo);
CAL=LowerRate_v2(cal0,P);


%%
BLK_pow_pool_db=rand(BLK_num,runtimes)*20-50;
BLK_amp_pool=sqrt(10.^(BLK_pow_pool_db./10));


%% Interference Generation and Tagging

for runindex=1:runtimes

    %  permutation of calibration sequences
    CAL=CAL(:,randperm(N));


    %  QPSK BLK generation with power scaling
    for ii=1:BLK_num
        tempindex=randi(BLKcan_num);
        BLK_r(:,ii)=real(sig_bb(:,tempindex))./sqrt(mean(real(sig_bb(:,tempindex)).^2));
        BLK_i(:,ii)=imag(sig_bb(:,tempindex))./sqrt(mean(imag(sig_bb(:,tempindex)).^2));
            for mm=1:M
                Starindex=randi(L-P-1);
                BLK((mm-1)*SampleNumberAve+1:mm*SampleNumberAve,ii)=BLK_r(Starindex:Starindex+SampleNumberAve-1,ii)+BLK_i(Starindex:Starindex+SampleNumberAve-1,ii);
            end
    end
    
    %  Power Scaling
    storedbb_h0=zeros(M*SampleNumberAve,1);
    storedbb_h1=zeros(M*SampleNumberAve,1);

    
    for bb=2:BLK_num            
        storedbb_h0=storedbb_h0+BLK(:,bb).*CAL(:,bb)*BLK_amp_pool(bb,runindex);
    end
    storedbb_h1=storedbb_h0+BLK(:,1).*CAL(:,1)*BLK_amp_pool(1,runindex);


    
    %  Tagging begin!
    %  Here we only care about BLK1   

%         BLK_power(runindex)=tagging_v3(BLK(:,1).*CAL(:,1)*BLK_amp_pool(1,runindex),CAL(:,1),SampleNumberAve,1);
%         h0_power(runindex)=tagging_v3(storedbb_h0,CAL(:,1),SampleNumberAve,1);
%         h1_power(runindex)=tagging_v3(storedbb_h1,CAL(:,1),SampleNumberAve,1);
% 
%     
        sig_MNave_BLK(runindex)=tagging_v3(BLK(:,1).*CAL(:,1)*BLK_amp_pool(1,runindex),CAL(:,1),SampleNumberAve,M);
        sig_MNave_h0(runindex)=tagging_v3(storedbb_h0,CAL(:,1),SampleNumberAve,M);
        sig_MNave_h1(runindex)=tagging_v3(storedbb_h1,CAL(:,1),SampleNumberAve,M);


end

%% TH comparison
TH_num=100;
TH_pool=linspace(0,15,TH_num)*1e-4;
tagged_table_h0=zeros(TH_num,runtimes);
tagged_table_h1=zeros(TH_num,runtimes);


for runindex=1:runtimes
    for TH_index=1:TH_num
        if sig_MNave_h0(runindex)>TH_pool(TH_index)
            tagged_table_h0(TH_index,runindex)=1;
        end
        if sig_MNave_h1(runindex)>TH_pool(TH_index)
            tagged_table_h1(TH_index,runindex)=1;
        end
    end
end


% ROC generation

    for TH_index=1:TH_num
        Pfa_MC(TH_index)=sum(tagged_table_h0(TH_index,:))/runtimes;
        Pd_MC(TH_index)=sum(tagged_table_h1(TH_index,:))/runtimes;
    end



% plot
color=['b';'r';'c';'g';'m'];
figure(99)

plot(Pfa_MC,Pd_MC,'linewidth',2);hold on

grid on


title('General Case of all BLKs, M=25')
xlabel('Pfa');
ylabel('Pd');
%
temp_h0=sort(sig_MNave_h0);

temp_h1=sort(sig_MNave_h1);

figure
subplot(211)
stem(temp_h0)
subplot(212)
stem(temp_h1)
% 
% %%

10*log10(mean(sig_MNave_BLK)...
    /mean(sig_MNave_h0))

% 
