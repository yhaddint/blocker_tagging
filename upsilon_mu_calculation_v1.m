% 03/08/2015
% mean and variance of both S^2 and N^2 has been compared between anal. and
% sim. 
% rrc shape will be take into consideration.

clear;
clc;
warning off


% Independent repeating experiment.
runtimes=1e3;

%N_range=[4];
%N_num = length(N_range);
M=100;

    % BLK
R_range = 32;
N_range=[31,63,127,255];
N_num=length(N_range);

Rindex = 1;  
SampleperSymbol=R_range(Rindex);
SampleperFrame=M*SampleperSymbol;


hdesign  = fdesign.pulseshaping(SampleperSymbol,'Square Root Raised Cosine');
hpulse = design(hdesign);
rrc=hpulse.Numerator./sqrt(mean(hpulse.Numerator.^2));%./sqrt(length(hpulse.Numerator));
for ii=1:length(rrc)
    h(ii)=rrc(ii:end)*rrc(1:end-ii+1)'/(length(rrc)-ii+1);
end

%  How many active BLK existing. To get reasonable number of different
%  envelope, we generate 2x of them.
BLK_num=5;
BLKcan_num=60;
SampleNumberAve = R_range(1);
%  Define length of BLK we want.
stat_num=1e3;
L=SampleperFrame*stat_num;
clearvars sig_bb BLK_r BLK_i BlK_pool
for ii=1:BLKcan_num
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
    unnormalized = real(temp_data(end-L+1:end));
    sig_r(:,ii) = unnormalized./sqrt(sum(unnormalized.^2)/L);
    unnormalized = imag(temp_data(end-L+1:end));
    sig_i(:,ii) = unnormalized./sqrt(sum(unnormalized.^2)/L);
    BLK_pool(:,ii) = (sig_r(:,ii)+sig_i(:,ii))/sqrt(2);

end
%%
for Nindex=1:N_num
    
N=N_range(Nindex);

BLK=zeros(N,2);

cal0=PNgenerator_v5(N,N,1);
CAL=LowerRate_v2(cal0,N);

%
for runindex=1:runtimes
    %  permutation of calibration sequences
    CAL_index1=randi(N);
    CAL_index2=randi(N);
    while CAL_index1==CAL_index2
        CAL_index1=randi(N);
        CAL_index2=randi(N);
    end

    tempindex=randi(BLKcan_num);

    Starindex=randi(size(BLK_pool,1)-N-1);
    
    BLK(:,1)=BLK_pool(Starindex:Starindex+N-1,tempindex);%+sig_bb_i(Starindex:Starindex+N-1,tempindex);
    %BLK(:,2)=BLK(Starindex:Starindex+N-1,tempindex);

    Power_Target(runindex) = tagging_v3(BLK(:,1),ones(size(BLK(:,1))),N,1);
    Power_NonTarget(runindex) = tagging_v3(BLK(:,1).*CAL(:,CAL_index1),CAL(:,CAL_index2),N,1);
    %Power_mixed(runindex) = tagging_v3(BLK(:,1).*CAL(:,CAL_index1)+BLK(:,2).*CAL(:,CAL_index2),CAL(:,CAL_index1),N,1);
end


mus2(Rindex,Nindex) = mean(Power_Target);
mun2(Rindex,Nindex) = mean(Power_NonTarget);

mus2_anal(Rindex,Nindex)=anal_mean_fun_v2(SampleperSymbol,N,h,'t');
mun2_anal(Rindex,Nindex)=anal_mean_fun_v2(SampleperSymbol,N,h,'n');


end


%%


figure;plot_setting();
subplot(211)
plot(N_range,mus2_anal(:,:));hold on
plot(N_range,mus2(:,:),'--');hold on
grid on
title('Mean of S_i^2')
xlabel('N')
%legend('N = 15 (Anal.)','N = 15 (Sim.)','N = 31 (Anal.)','N = 31 (Sim.)','N = 63 (Anal.)','N = 63 (Sim.)')



subplot(212)
plot(N_range,mun2_anal);hold on
plot(N_range,mun2,'--');hold on
grid on
title('Mean of N_i^2')
xlabel('N')
%legend('N = 15 (Anal.)','N = 15 (Sim.)','N = 31 (Anal.)','N = 31 (Sim.)','N = 63 (Anal.)','N = 63 (Sim.)')



