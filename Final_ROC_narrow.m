%  09/24/2014
%  adapt to Narrow CAL model. Fix problem on Power/Amplitude
%  Now complex blockers are considered.

%clearvars -Pfa_eq -Pd_eq;
clc;
warning off

worstcase=0;

% Independent repeating experiment.
runtimes=2e2;

%  parameter M and N setting. N is also length of CAL Sequences
N=15;
M=25;
downSampleTo=1;

SampleNumberAve=N*downSampleTo;

P=M*SampleNumberAve;

%  How many BLK we want to tag simutaneously
tag_num=20;

%  How many active BLK existing. To get reasonable number of different
%  envelope, we generate 3x of them.
if worstcase
    BLK_num=6;
else
    BLK_num=14;
end

BLKcan_num=20;

% We assume simplified model hold (Strong CAL).
lowerThanPNdB=40;
lowerThanPN=10^(lowerThanPNdB/10);

%  Define length of BLK we want.
stat_num=40;
L=M*SampleNumberAve*stat_num;

%  Averaging Time = xT


%  Normalized BLK generation. 
%  CAL BW/ BLK BW = 100.
for ii=1:BLKcan_num
    upsam=SampleNumberAve*2;
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

%% calibration sequences
cal0 = PN_generator_downsample(N,N,downSampleTo);
CAL = LowerRate(cal0, 1, P);

%% initialization of tagging tables

TH_num=100;
tagged_table_h0=zeros(TH_num,stat_num,runtimes);
tagged_table_h1=zeros(TH_num,stat_num,runtimes);


%% Interference Generation and Tagging

if worstcase
    BLK_pow_pool_dB=ones(BLK_num,BLK_num-1)*(-40);
    BLK_pow_pool_dB(1,:)=ones(1,BLK_num-1)*(-10);
    for bb=1:BLK_num-1
        BLK_pow_pool_dB(2:bb+1,bb)=ones(bb,1)*(10);
    end
else
    BLK_pow_pool_dB=ones(BLK_num,(BLK_num-2)/4+1)*(-40);
    BLK_pow_pool_dB(1,:)=ones(1,(BLK_num-2)/4+1)*(-10);
    for bb=2:BLK_num
        BLK_pow_pool_dB(bb,fix((bb+5)/4):end)=3;
    end
end


BLK_pow_pool=(10.^(BLK_pow_pool_dB./10))./lowerThanPN;

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
    storedbb_h0=zeros(M*SampleNumberAve,BLK_num);
    storedbb_h1=zeros(M*SampleNumberAve,BLK_num);

    for pp=1:size(BLK_pow_pool,2)
        for bb=2:BLK_num            
            storedbb_h0(:,pp)=storedbb_h0(:,pp)+BLK(:,bb).*CAL(:,bb)*sqrt(BLK_pow_pool(bb,pp));
        end
        storedbb_h1(:,pp)=storedbb_h0(:,pp)+BLK(:,1).*CAL(:,1)*sqrt(BLK_pow_pool(1,pp));
    end
    
    
    %  Tagging begin!
    %  Here we only care about BLK1   
    for pp=1:size(BLK_pow_pool,2)
        sig_MNave_BLK(runindex,pp) = correlation_and_pow(BLK(:,1).*CAL(:,1)*sqrt(BLK_pow_pool(1,pp)),CAL(:,1),SampleNumberAve,M);
        sig_MNave_h0(runindex,pp) = correlation_and_pow(storedbb_h0(:,pp),CAL(:,1),SampleNumberAve,M);
        sig_MNave_h1(runindex,pp) = correlation_and_pow(storedbb_h1(:,pp),CAL(:,1),SampleNumberAve,M);
    end

end

%% TH comparison
TH_pool=linspace(0,2,TH_num)*1e-4;
tagged_table_h0=zeros(TH_num,runtimes,BLK_num-1);
tagged_table_h1=zeros(TH_num,runtimes,BLK_num-1);

for pp=1:size(BLK_pow_pool,2)
    for runindex=1:runtimes
        for TH_index=1:TH_num
            if sig_MNave_h0(runindex,pp)>TH_pool(TH_index)
                tagged_table_h0(TH_index,runindex,pp)=1;
            end
            if sig_MNave_h1(runindex,pp)>TH_pool(TH_index)
                tagged_table_h1(TH_index,runindex,pp)=1;
            end
        end
    end


% ROC generation

    for TH_index=1:TH_num
        Pfa_MC(TH_index,pp)=sum(tagged_table_h0(TH_index,:,pp))/runtimes;
        Pd_MC(TH_index,pp)=sum(tagged_table_h1(TH_index,:,pp))/runtimes;
    end
end


%% plot
color=['b';'r';'c';'g';'m'];
figure(99)
for pp=1:size(BLK_pow_pool,2)
plot(Pfa_MC(:,pp),Pd_MC(:,pp),color(pp,:),'linewidth',2);hold on
end
grid on

if worstcase
    legend('2 Active BLK','3 Active BLK','4 Active BLK','5 Active BLK','6 Active BLK');
else
    legend('2 Active BLK','6 Active BLK','10 Active BLK','14 Active BLK','18 Active BLK');   
end
title('-50dBm BLK with -37dBm Interference, M=25')
xlabel('Pfa');
ylabel('Pd');
%%
temp_h0=sort(sig_MNave_h0(:,4));

temp_h1=sort(sig_MNave_h1(:,4));

figure
subplot(211)
stem(temp_h0)
subplot(212)
stem(temp_h1)

%%
for ii=1:4
10*log10(mean(sig_MNave_BLK(:,ii))...
    /mean(sig_MNave_h0(:,ii)))
end

