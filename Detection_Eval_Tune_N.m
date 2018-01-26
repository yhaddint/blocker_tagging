%  02/25/2015
%  pfa is kept at 0.1 to see pd in different of M (averaging time for
%  energy detection)
% The figure used in ARR2015 poster
clear;
clc;
warning off


% Independent repeating experiment.
runtimes=1e3;

N_range=[15,31,63];
for nindex=1:3
    clearvars -except runtimes N_range nindex pd_withFixed_pfa
%  parameter M and N setting. N is also length of CAL Sequences
N=N_range(nindex);
M=50;
downSampleTo=1;

SampleNumberAve=N*downSampleTo;

P=M*SampleNumberAve;

%  How many BLK we want to tag simutaneously
tag_num=20;

%  How many active BLK existing. To get reasonable number of different
%  envelope, we generate 2x of them.
BLK_num=10;
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
cal0=PN_generator_downsample(N,N,downSampleTo);
CAL=LowerRate(cal0,1,P);

%% initialization of tagging tables

TH_num=100;

%% Interference Generation and Tagging
%powerdiff_range=[2,5,8];
powerdiff_range=[8];


%% tune M
mm_range=2:4:M;
mm_num=length(mm_range);

%%
for pindex=1:length(powerdiff_range)
    

    BLK_pow_pool_dB=ones(BLK_num,1)*(-40);
    BLK_pow_pool_dB(1,:)=ones(1,1)*(-10);

    BLK_pow_pool_dB(2:10)=ones(9,1)*powerdiff_range(pindex);

    %BLK_pow_pool_dB = fliplr(BLK_pow_pool_dB);



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
            storedbb_h0(:,pp)=storedbb_h0(:,pp)+BLK(:,bb).*CAL(:,bb)*sqrt(BLK_pow_pool(bb));
        end
        storedbb_h1(:,pp)=storedbb_h0(:,pp)+BLK(:,1).*CAL(:,1)*sqrt(BLK_pow_pool(1));
    end
    
    
    %  Tagging begin!
    %  Here we only care about BLK1   
    for pp=1:size(BLK_pow_pool,2)
        for mm=1:mm_num
            M_sweep=mm_range(mm);
            %sig_MNave_BLK(runindex,mm)=tagging_v3(BLK(1:M*SampleNumberAve,1).*CAL(1:M*SampleNumberAve,1)*sqrt(BLK_pow_pool(1,pp)),CAL(1:M*SampleNumberAve,1),SampleNumberAve,M);
            sig_MNave_h0(runindex,mm) = correlation_and_pow(storedbb_h0(1:M_sweep*SampleNumberAve,pp),CAL(1:M_sweep*SampleNumberAve,1),SampleNumberAve,M_sweep);
            sig_MNave_h1(runindex,mm) = correlation_and_pow(storedbb_h1(1:M_sweep*SampleNumberAve,pp),CAL(1:M_sweep*SampleNumberAve,1),SampleNumberAve,M_sweep);
        end
    end


end

%% TH comparison

Pfa_MC=zeros(TH_num,mm_num);
Pd_MC=zeros(TH_num,mm_num);

for mm=1:mm_num
    TH_pool=linspace(min(sig_MNave_h0(:,mm)),max(sig_MNave_h1(:,mm)),TH_num);
    for TH_index=1:TH_num
        Pfa_MC(TH_index,mm) = sum(sig_MNave_h0(:,mm)>TH_pool(TH_index))/runtimes;
        Pd_MC(TH_index,mm)=sum(sig_MNave_h1(:,mm)>TH_pool(TH_index))/runtimes;
    end
end


%%  fixed Pfa = 0.1
for mm=1:mm_num
    [value,index]=min(abs(Pfa_MC(:,mm)-0.1));
    pd_withFixed_pfa(mm,nindex)=Pd_MC(index,mm);
end

end

end
%%
% get Pd performance with Pfa = 0.1 in each cases
color=['b-o';'r-o';'c-o';'g-o';'m-o'];
xmm=mm_range(1:1:end);
figure
plot_setting()
for nindex=1:3
plot(xmm,pd_withFixed_pfa(1:1:end,nindex),color(nindex,:));hold on
end
grid on
legend('N = 15','N = 31','N = 63')
ylabel('Probability of Detection')
xlabel('M')
ylim([0,1.01])
% %%
% figure;
% subplot(211)
% plot(sig_MNave_h0(:,25))
% subplot(212)
% plot(sig_MNave_h1(:,25))