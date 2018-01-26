% 03/08/2015
% mean and variance of both S^2 and N^2 has been compared between anal. and
% sim. 
% Square waveform is considered here.

clear;
clc;
warning off


% Independent repeating experiment.
runtimes=1e3;

%N_range=[4];
%N_num = length(N_range);
M=50;

    % BLK
R_range = 10:10:100;
N_range=[15,31,63];
N_num=length(N_range);
for Rindex=1:length(R_range)
 Rindex   
SampleperSymbol=R_range(Rindex);
SampleperFrame=M*SampleperSymbol;


%  How many BLK we want to tag simutaneously
tag_num=20;

%  How many active BLK existing. To get reasonable number of different
%  envelope, we generate 2x of them.
BLK_num=10;
BLKcan_num=20;

%  Define length of BLK we want.
stat_num=40;
L=SampleperFrame*stat_num;
clearvars sig_bb_r sig_bb_i
for ii=1:BLKcan_num
    sig_bb_r(:,ii) = kron(randi(2,M*stat_num,1)*2-3,ones(SampleperSymbol,1));
    sig_bb_i(:,ii) = kron(randi(2,M*stat_num,1)*2-3,ones(SampleperSymbol,1));
end

for Nindex=1:N_num
    
N=N_range(Nindex);

BLK=zeros(N,2);

cal0=PN_generator_downsample(N,N,1);
CAL=LowerRate(cal0,1,N);

%%
for runindex=1:runtimes
    %  permutation of calibration sequences
    CAL_index1=randi(N);
    CAL_index2=randi(N);
    while CAL_index1==CAL_index2
        CAL_index1=randi(N);
        CAL_index2=randi(N);
    end

    tempindex=randi(BLKcan_num);

    Starindex=randi(L-SampleperFrame-1);
    
    BLK(:,1)=sig_bb_r(Starindex:Starindex+N-1,tempindex);%+sig_bb_i(Starindex:Starindex+N-1,tempindex);
    BLK(:,2)=sig_bb_i(Starindex:Starindex+N-1,tempindex);

    Power_Target(runindex) = correlation_and_pow(BLK(:,1),ones(size(BLK(:,1))),N,1);
    Power_NonTarget(runindex) = correlation_and_pow(BLK(:,1).*CAL(:,CAL_index1),CAL(:,CAL_index2),N,1);
    Power_mixed(runindex) = correlation_and_pow(BLK(:,1).*CAL(:,CAL_index1)+BLK(:,2).*CAL(:,CAL_index2),CAL(:,CAL_index1),N,1);
end


mus2(Rindex,Nindex) = mean(Power_Target);
mun2(Rindex,Nindex) = mean(Power_NonTarget);

mus2_anal(Rindex,Nindex)=anal_mean_fun(SampleperSymbol,N,'t');
mun2_anal(Rindex,Nindex)=anal_mean_fun(SampleperSymbol,N,'n');

sigman2(Rindex,Nindex) = var(Power_NonTarget);
sigmas2(Rindex,Nindex) = var(Power_Target);

sigmas2_anal(Rindex,Nindex) = sigma_s2_fun(N,SampleperSymbol)-mus2_anal(Rindex,Nindex)^2;
sigman2_anal(Rindex,Nindex) = sigma_n2_fun(N,SampleperSymbol)-mun2_anal(Rindex,Nindex)^2;
end
end

%% mean analysical value
clearvars mus2_anal mu02_anal
R_range_fine=linspace(min(R_range),max(R_range),length(R_range));
R_num_fine = length(R_range_fine);
for Rindex=1:R_num_fine
    for Nindex=1:N_num
        R=R_range_fine(Rindex);
        N=N_range(Nindex);
        mus2_anal(Rindex,Nindex)=anal_mean_fun(R,N,'t');
        mu02_anal(Rindex,Nindex)=anal_mean_fun(R,N,'n');
    end
end
    
%% plot
color_anal=['b';'r';'c'];
color_sim=['bx';'rx';'cx'];

figure;plot_setting();
subplot(221)
for Nindex=1:N_num
    plot(R_range_fine,mus2_anal(:,Nindex),color_anal(Nindex,:));hold on
    plot(R_range,mus2(:,Nindex),color_sim(Nindex,:));hold on
end
grid on
title('Mean of S_i^2')
xlabel('R')
legend('N = 15 (Anal.)','N = 15 (Sim.)','N = 31 (Anal.)','N = 31 (Sim.)','N = 63 (Anal.)','N = 63 (Sim.)')

subplot(222)
for Nindex=1:N_num
    plot(R_range,sigmas2_anal(:,Nindex),color_anal(Nindex,:));hold on
    plot(R_range,sigmas2(:,Nindex),color_sim(Nindex,:));hold on
end
grid
title('Variance of S_i^2')
xlabel('R')
legend('N = 15 (Anal.)','N = 15 (Sim.)','N = 31 (Anal.)','N = 31 (Sim.)','N = 63 (Anal.)','N = 63 (Sim.)')


subplot(223)
for Nindex=1:N_num
    plot(R_range_fine,mun2_anal(:,Nindex),color_anal(Nindex,:));hold on
    plot(R_range,mun2(:,Nindex),color_sim(Nindex,:));hold on
end
grid on
title('Mean of N_i^2')
xlabel('R')
legend('N = 15 (Anal.)','N = 15 (Sim.)','N = 31 (Anal.)','N = 31 (Sim.)','N = 63 (Anal.)','N = 63 (Sim.)')

subplot(224)
for Nindex=1:N_num
    plot(R_range,sigman2_anal(:,Nindex),color_anal(Nindex,:));hold on
    plot(R_range,sigman2(:,Nindex),color_sim(Nindex,:));hold on
end
grid
title('Variance of N_i^2')
xlabel('R')
legend('N = 15 (Anal.)','N = 15 (Sim.)','N = 31 (Anal.)','N = 31 (Sim.)','N = 63 (Anal.)','N = 63 (Sim.)')
%legend('N = 15 (Sim.)','N = 31 (Sim.)','N = 63 (Sim.)')


