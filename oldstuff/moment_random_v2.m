% 05/14/2015
% mean and variance of both S^2 and N^2 has been compared between anal. and
% sim. 
% Square waveform is considered here.
% 4-QAM, 16QAM, 64QAM are compared here.
% here random codes are used instead of PN codes.

clear;
clc;clf;close all
warning off


% Independent repeating experiment.
runtimes=1e4;

%N_range=[4];
%N_num = length(N_range);
M=50;

    % BLK
R_range = 10:10:200;
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

%%
for runindex=1:runtimes
    %  permutation of calibration sequences
    CAL=randi(2,N,2)*2-3;
    
    tempindex=randi(BLKcan_num);

    Starindex=randi(L-SampleperFrame-1);
    
    BLK(:,1)=(sig_bb_r(Starindex:Starindex+N-1,tempindex)+sig_bb_i(Starindex:Starindex+N-1,tempindex))*sqrt(2)/2;
    BLK(:,2)=sig_bb_i(Starindex:Starindex+N-1,tempindex);

    Power_Target(runindex) = tagging_v3(BLK(:,1),ones(size(BLK(:,1))),N,1);
    Power_NonTarget(runindex) = tagging_v3(BLK(:,1).*CAL(:,1),CAL(:,2),N,1);
    %Power_mixed(runindex) = tagging_v3(BLK(:,1).*CAL(:,1)+BLK(:,2).*CAL(:,1),CAL(:,2),N,1);
end


mus2(Rindex,Nindex) = mean(Power_Target);
mun2(Rindex,Nindex) = mean(Power_NonTarget);

% mus2_anal(Rindex,Nindex)=anal_mean_fun(SampleperSymbol,N,'t');
% mun2_anal(Rindex,Nindex)=anal_mean_fun(SampleperSymbol,N,'n');

sigman2(Rindex,Nindex) = var(Power_NonTarget);
sigmas2(Rindex,Nindex) = var(Power_Target);

% sigmas2_anal(Rindex,Nindex) = sigma_s2_fun(N,SampleperSymbol)-mus2_anal(Rindex,Nindex)^2;
% sigman2_anal(Rindex,Nindex) = sigma_n2_fun(N,SampleperSymbol)-mun2_anal(Rindex,Nindex)^2;
end
end

  
%% plot
color_anal=['b';'r';'c'];
color_sim=['bx';'rx';'cx'];

figure(1);
% plot_setting
subplot(221)
for Nindex=1:N_num
    %plot(R_range_fine,mus2_anal(:,Nindex),color_anal(Nindex,:));hold on
    plot(R_range,mus2(:,Nindex),color_sim(Nindex,:));hold on
end
grid on
title('Mean of S_i^2')
xlabel('R')
legend('N = 15 (Anal.)','N = 15 (Sim.)','N = 31 (Anal.)','N = 31 (Sim.)','N = 63 (Anal.)','N = 63 (Sim.)')

subplot(222)
for Nindex=1:N_num
    %plot(R_range,sigmas2_anal(:,Nindex),color_anal(Nindex,:));hold on
    plot(R_range,sigmas2(:,Nindex),color_sim(Nindex,:));hold on
end
grid
title('Variance of S_i^2')
xlabel('R')
legend('N = 15 (Anal.)','N = 15 (Sim.)','N = 31 (Anal.)','N = 31 (Sim.)','N = 63 (Anal.)','N = 63 (Sim.)')


subplot(223)
for Nindex=1:N_num
    %plot(R_range_fine,mun2_anal(:,Nindex),color_anal(Nindex,:));hold on
    plot(R_range,mun2(:,Nindex),color_sim(Nindex,:));hold on
end
grid on
title('Mean of N_i^2')
xlabel('R')
legend('N = 15 (Anal.)','N = 15 (Sim.)','N = 31 (Anal.)','N = 31 (Sim.)','N = 63 (Anal.)','N = 63 (Sim.)')

subplot(224)
for Nindex=1:N_num
    %plot(R_range,sigman2_anal(:,Nindex),color_anal(Nindex,:));hold on
    plot(R_range,sigman2(:,Nindex),color_sim(Nindex,:));hold on
end
grid
title('Variance of N_i^2')
xlabel('R')
legend('N = 15 (Anal.)','N = 15 (Sim.)','N = 31 (Anal.)','N = 31 (Sim.)','N = 63 (Anal.)','N = 63 (Sim.)')
%legend('N = 15 (Sim.)','N = 31 (Sim.)','N = 63 (Sim.)')

%%
for Nindex=1:N_num
    rs(:,Nindex) = sigmas2(:,Nindex)-2*mus2(:,Nindex).^2;
    r0(:,Nindex) = sigman2(:,Nindex)-2*mun2(:,Nindex).^2;
end

figure;
for Nindex=1:N_num
    plot(R_range,rs(:,Nindex),'b');hold on
    plot(R_range,r0(:,Nindex),'r');hold on
    
end