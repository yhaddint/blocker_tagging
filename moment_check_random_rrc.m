% 05/17/2015
% mean and variance of both S^2 and N^2 has been compared between anal. and
% sim. 
% rrc shape will be take into consideration.
% use random code rather than PN codes

clear;
clc;
warning off


% Independent repeating experiment.
runtimes=3e4;

%N_range=[4];
%N_num = length(N_range);
M=100;

    % BLK
R_range = 10:20:210;
%R_range = 30;
N_range=[15,31,63];
N_num=length(N_range);
for Rindex=1:length(R_range)
 Rindex
 clearvars sig
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
BLK_num=20;
BLKcan_num=20;

%  Define length of BLK we want.
stat_num=40;
L=SampleperFrame*stat_num;
clearvars sig_bb BLK_r BLK_i
for ii=1:BLKcan_num
    symbols=200;
    hmod = modem.pskmod('M', 4, 'InputType', 'integer');
    data = randi(4,symbols,1)-1;
    data = modulate(hmod, data);
    data = upsample(data,SampleperSymbol);
    temp_data = conv(data,hpulse.Numerator);
    sig_bb(:,ii) = temp_data(end-195*SampleperSymbol+1:end-5*SampleperSymbol+1);
     
    BLK_r(:,ii)=real(sig_bb(:,ii))./sqrt(mean(real(sig_bb(:,ii)).^2));
    BLK_i(:,ii)=imag(sig_bb(:,ii))./sqrt(mean(imag(sig_bb(:,ii)).^2));
    sig(:,ii)=(BLK_r(:,ii)+BLK_i(:,ii))/sqrt(2);
end

for Nindex=1:N_num
    
N=N_range(Nindex);

BLK=zeros(N,2);

%%
for runindex=1:runtimes
    CAL=randi(2,N,2)*2-3;
    Starindex=randi(size(BLK_r,1)-N-1);
    tempindex=randi(BLKcan_num);
    
    BLK(:,1)=sig(Starindex:Starindex+N-1,tempindex);%+sig_bb_i(Starindex:Starindex+N-1,tempindex);
    BLK(:,2)=sig(Starindex:Starindex+N-1,tempindex);

    Power_Target(runindex) = tagging_v3(BLK(:,1),ones(size(BLK(:,1))),N,1);
    Power_NonTarget(runindex) = tagging_v3(BLK(:,1).*CAL(:,1),CAL(:,2),N,1);
%     Power_mixed(runindex) = tagging_v3(BLK(:,1).*CAL(:,1)+BLK(:,2).*CAL(:,CAL_index2),CAL(:,CAL_index1),N,1);
end


mus2(Rindex,Nindex) = mean(Power_Target);
mun2(Rindex,Nindex) = mean(Power_NonTarget);

% mus2_anal(Rindex,Nindex)=anal_mean_fun_v2(SampleperSymbol,N,h,'t');
% mun2_anal(Rindex,Nindex)=anal_mean_fun_v2(SampleperSymbol,N,h,'n');

sigman2(Rindex,Nindex) = var(Power_NonTarget);
sigmas2(Rindex,Nindex) = var(Power_Target);

%sigmas2_anal(Rindex,Nindex) = sigma_s2_fun(N,SampleperSymbol)-mus2_anal(Rindex,Nindex)^2;
%sigman2_anal(Rindex,Nindex) = sigma_n2_fun(N,SampleperSymbol)-mun2_anal(Rindex,Nindex)^2;
end
end

%% mean analysical value
% clearvars mus2_anal mu02_anal
% R_range_fine=linspace(min(R_range),max(R_range),length(R_range));
% R_num_fine = length(R_range_fine);
% for Rindex=1:R_num_fine
%     for Nindex=1:N_num
%         R=R_range_fine(Rindex);
%         N=N_range(Nindex);
%         mus2_anal(Rindex,Nindex)=anal_mean_fun(R,N,'t');
%         mu02_anal(Rindex,Nindex)=anal_mean_fun(R,N,'n');
%     end
% end
%     
%% plot
color_anal=['b';'r';'c'];
color_sim=['bx';'rx';'cx'];

figure;plot_setting();
subplot(221)
for Nindex=1:N_num
%     plot(R_range,mus2_anal(:,Nindex),color_anal(Nindex,:));hold on
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
%     plot(R_range,mun2_anal(:,Nindex),color_anal(Nindex,:));hold on
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

