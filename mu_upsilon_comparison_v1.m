% 05/21/2015
% theoretically evaluate \mu and \upsilon 
% rrc shape will be take into consideration.
% comparison: MC, Anal, Bound
clear;
clc;
warning off


% Independent repeating experiment.
runtimes=2e4;
rrcPulseShape=0;

%N_range=[4];
%N_num = length(N_range);
M=100;

% BLK
R_range = fliplr([2 4 8 16 32 64 128 256 512]);
N_range=[15,31,63];
%N_range=31;
N_num=length(N_range);
for Rindex=1:length(R_range)
 Rindex   
SampleperSymbol=R_range(Rindex);
SampleperFrame=M*SampleperSymbol;


hdesign  = fdesign.pulseshaping(SampleperSymbol,'Square Root Raised Cosine',30,1);
%hdesign  = fdesign.pulseshaping(SampleperSymbol,'Raised Cosine','Nsym,Beta',22,1);

hpulse = design(hdesign);
rrc=hpulse.Numerator./sqrt(mean(hpulse.Numerator.^2));%./sqrt(length(hpulse.Numerator));
for ii=1:length(rrc)
    h(ii)=rrc(ii:end)*rrc(1:end-ii+1)'/(length(rrc));
end

%  How many active BLK existing. To get reasonable number of different
%  envelope, we generate 2x of them.
BLK_num=50;
BLKcan_num=50;

%  Define length of BLK we want.
stat_num=40;
L=SampleperFrame*stat_num;
if rrcPulseShape==1
    clearvars BLK_r BLK_i BLK_sum sig_bb
    for ii=1:BLKcan_num
        symbols=500;
        hmod = modem.qammod('M', 4, 'InputType', 'integer');
        data = randi(4,symbols,1)-1;
        data = modulate(hmod, data);
        data = upsample(data,SampleperSymbol);
        temp_data = conv(data,hpulse.Numerator);
        sig_bb(:,ii) = temp_data(end-(symbols-5)*SampleperSymbol+1:end-5*SampleperSymbol+1);

        BLK_r(:,ii)=real(sig_bb(:,ii))./sqrt(mean(real(sig_bb(:,ii)).^2));
        BLK_i(:,ii)=imag(sig_bb(:,ii))./sqrt(mean(imag(sig_bb(:,ii)).^2));
        BLK_sum(:,ii) = (BLK_r(:,ii)+BLK_i(:,ii))/sqrt(2);
    end
else
    clearvars sig_bb_r sig_bb_i BLK_sum
    for ii=1:BLKcan_num
        sig_bb_r(:,ii) = kron(randi(2,M*stat_num,1)*2-3,ones(SampleperSymbol,1));
        sig_bb_i(:,ii) = kron(randi(2,M*stat_num,1)*2-3,ones(SampleperSymbol,1));
        BLK_sum(:,ii) = (sig_bb_r(:,ii)+sig_bb_i(:,ii))/sqrt(2);
    end
end

for Nindex=1:N_num
    
N=N_range(Nindex);

BLK=zeros(N,2);

cal0=PNgenerator_v5(N,N,1);
CAL=LowerRate_v2(cal0,N);

%%
for runindex=1:runtimes
    %  permutation of calibration sequences
    CAL_index1=randi(N);
    CAL_index2=randi(N);
    while CAL_index1==CAL_index2
        CAL_index1=randi(N);
        CAL_index2=randi(N);
    end
    Code = randi(2,N,2)*2-3;
    tempindex=randi(BLKcan_num);

    Starindex=randi(size(BLK_sum,1)-N-1);
    
    BLK(:,1)=BLK_sum(Starindex:Starindex+N-1,tempindex);%+sig_bb_i(Starindex:Starindex+N-1,tempindex);
    BLK(:,2)=BLK_sum(Starindex:Starindex+N-1,tempindex);

    Power_Target(runindex) = tagging_v3(BLK(:,1),ones(size(BLK(:,1))),N,1);
    Power_NonTarget(runindex) = tagging_v3(BLK(:,1).*CAL(:,CAL_index1),CAL(:,CAL_index2),N,1);
    Power_mixed(runindex) = tagging_v3(BLK(:,1).*CAL(:,CAL_index1)+BLK(:,2).*CAL(:,CAL_index2),CAL(:,CAL_index1),N,1);
    Power_NonTarget_randomcode(runindex) = tagging_v3(BLK(:,1).*Code(:,1),Code(:,2),N,1);
end


mus2(Rindex,Nindex) = mean(Power_Target);
mun2(Rindex,Nindex) = mean(Power_NonTarget);
mun2_randomcode(Rindex,Nindex) = mean(Power_NonTarget_randomcode);

if rrcPulseShape
    mus2_anal(Rindex,Nindex)=anal_mean_fun_v2(SampleperSymbol,N,h,'t');
    mun2_anal(Rindex,Nindex)=anal_mean_fun_v2(SampleperSymbol,N,h,'n');
else
    mus2_anal(Rindex,Nindex)=anal_mean_fun(SampleperSymbol,N,'t');
    mun2_anal(Rindex,Nindex)=anal_mean_fun(SampleperSymbol,N,'n');  
end
%sigman2(Rindex,Nindex) = var(Power_NonTarget);
%sigmas2(Rindex,Nindex) = var(Power_Target);

%sigmas2_anal(Rindex,Nindex) = sigma_s2_fun(N,SampleperSymbol)-mus2_anal(Rindex,Nindex)^2;
%sigman2_anal(Rindex,Nindex) = sigma_n2_fun(N,SampleperSymbol)-mun2_anal(Rindex,Nindex)^2;
end
end
    
%% plot
color_anal_s=['b--';'r--';'c--'];
color_sim_s=['bs';'r^';'co'];
color_anal_n=['b';'r';'c'];
color_sim_n=['b+';'rp';'cx'];
color_bound=['b-.';'r-.';'c-.';];

color_sim_rc=['bd';'r.';'c*'];
color_anal_rc=['b:';'r:';'c:'];

figure;plot_setting();
for Nindex=1:N_num
    plot(R_range,mus2_anal(:,Nindex),color_anal_s(Nindex,:));hold on
    plot(R_range,mus2(:,Nindex),color_sim_s(Nindex,:));hold on
end

for Nindex=1:N_num
    plot(R_range,mun2_anal(:,Nindex),color_anal_n(Nindex,:));hold on
    plot(R_range,mun2(:,Nindex),color_sim_n(Nindex,:));hold on
    plot(R_range,ones(1,length(R_range))*1/N_range(Nindex)^2,color_bound(Nindex,:));hold on
    plot(R_range,mun2_randomcode(:,Nindex),color_sim_rc(Nindex,:));hold on
    plot(R_range,(1/N_range(Nindex))*ones(1,length(R_range)),color_anal_rc(Nindex,:));hold on
end
grid on
xlabel('Interferer Bandwidth');
legend('\upsilon, N = 15 (Anal.)','\upsilon, N = 15 (Sim.)',...
    '\upsilon, N = 31 (Anal.)','\upsilon, N = 31 (Sim.)',...
    '\upsilon, N = 63 (Anal.)','\upsilon, N = 63 (Sim.)',...
    '\mu, N = 15 (Anal.)','\mu, N = 15 (Sim.)','\mu, N = 15 (Lim.)', '\mu, N = 15 (RC.Anal.)', '\mu, N = 15 (RC.Sim.)',...
    '\mu, N = 31 (Anal.)','\mu, N = 31 (Sim.)','\mu, N = 31 (Lim.)', '\mu, N = 31 (RC.Anal.)', '\mu, N = 31 (RC.Sim.)',...
    '\mu, N = 63 (Anal.)','\mu, N = 63 (Sim.)','\mu, N = 63 (Lim.)', '\mu, N = 63 (RC.Anal.)', '\mu, N = 63 (RC.Sim.)');


