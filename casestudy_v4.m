%  06/19/2015
% case study. its simumation setting details is listted in the slides
% here power of interferers are uniformly random without scaling
% and P_TH is not fixed either
clear;clc;%clf;close all


runtimes_sim = 1e4;
Pc=1000;


blk_num = 31;
M = 30;
N = 31;
SampleNumberAve = 32;
stat_num = 5e2;
P = M*N;
L = M*SampleNumberAve*stat_num;
blk_order=1:blk_num;

for ii = 1:blk_num
%     sig_r(:,ii) = kron(randi(2,M*stat_num,1)*2-3,ones(SampleNumberAve,1));
%     sig_i(:,ii) = kron(randi(2,M*stat_num,1)*2-3,ones(SampleNumberAve,1));
%     sig(:,ii) = (sig_r(:,ii)+sig_i(:,ii))/sqrt(2);
    upsam = SampleNumberAve;
    symbols = fix(L*2/upsam);
    
    clearvars data temp_data
    hmod = modem.pskmod('M', 4, 'InputType', 'integer');
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
    sig_i(:,ii) = zeros(size(sig_r(:,ii)));
    sig(:,ii) = (sig_r(:,ii)+sig_i(:,ii))/sqrt(2);

end

%% calibration sequences
% cal0=PNgenerator_v5(N,N,1);
% CAL_noperm=LowerRate_v2(cal0,P);
% CAL_perm=zeros(P,N);
% for mm=1:M
%     CAL_perm((mm-1)*N+1:mm*N,:)=CAL_noperm((mm-1)*N+1:mm*N,randperm(N));
% end

%% tagging: first stage
pow_req_dB=-20;
stage_total=0;
result_all = zeros(1,runtimes_sim*20);
tag_table=zeros(runtimes_sim,blk_num);
counter=1;
for runindex=1:runtimes_sim
    CAL_perm=randi(2,N*M,N)*2-3;
    
    result = zeros(blk_num,3);
    
    sig_pow_dB(:,runindex) = sort(rand(31,1)*90-90,'descend');
    sig_pow_temp = 10.^(sig_pow_dB(:,runindex)/10);
    sig_pow = sig_pow_temp./sum(sig_pow_temp);
%     sig_pow = sig_pow_temp;
    
    start_point = randi(L-P-1,blk_num);
    for ii=1:blk_num
%         for mm=1:M
%             sig_cache((mm-1)*N+1:mm*N,ii) = sig(start_point(ii,mm):start_point(ii,mm)+N-1,ii).*sqrt(sig_pow(ii));
%         end
            sig_cache(:,ii) = sig(start_point(ii):start_point(ii)+P-1,ii).*sqrt(sig_pow(ii));

    end
    
    stage(runindex)=1;
    endflag=0;
    while endflag==0
        % measure power summation
        pow_sum=0;
        for ii=1:blk_num
            if tag_table(runindex,ii)~=1
                pow_sum=pow_sum+sig_pow(ii);
            end
        end
        
        % power estimation
        sig_mix_perm = zeros(P,1);
        for ii=1:blk_num
            if tag_table(runindex,ii)~=1
                sig_mix_perm = sig_mix_perm + sig_cache(:,ii).*CAL_perm(:,ii)+0.5/sqrt(Pc)*sig_cache(:,ii).^2;
            end
        end
        
        pow_sum_hat = sum(sig_mix_perm.^2)/P*2;
        % update P_TH 
        if 10*log10(pow_sum)-9<=pow_req_dB
            P_TH=10^((pow_req_dB-4)/10);
            endflag=1;
        else
            P_TH=pow_sum*0.0265*2;
        end
% %         
        % update P_TH

%         P_TH=pow_sum*0.0398;
        % r0= [8,16,32,64], N=31
%         S0 = [0.2386,0.4378,0.7165,0.9054];
        % r0=32, N = 63
%         S0 = 0.4319;
%         P_TH=(pow_sum/N+0.005*S0)*(qfuncinv(0.05)/sqrt(M/2)+1)/2;
%         P_TH=(pow_sum_hat/N+0.005)*(qfuncinv(0.05)/sqrt(M/2)+1)/2;

        

        for ii=1:blk_num
            result(ii,stage(runindex)) = tagging_v3(sig_mix_perm,CAL_perm(1:N*M,ii),N,M);
        end
        % detection
        tag_table(runindex,:) = tag_table(runindex,:)+(result(:,stage(runindex))>=P_TH)';
        % control while
        stage(runindex)=stage(runindex)+1;
        stage_total=stage_total+1;
        if stage(runindex)==5
            endflag=1;
        end
    end
    
    % detection summary
    tagged_blk = find(blk_order.*tag_table(runindex,:));
    
    result_all(counter:counter+length(tagged_blk)-1)=sig_pow_dB(tagged_blk,runindex);
    counter = counter+length(tagged_blk);
end
%%
result_all_nz = result_all(find(result_all));
for ii=1:41
    TH=-(ii-1);
    count_total(ii) = sum(sum(sig_pow_dB(:,:)<=TH));
    count_tagged(ii) = sum(result_all_nz<=TH);
end
for ii=1:40
    ydata(ii) = (count_tagged(ii)-count_tagged(ii+1))/(count_total(ii)-count_total(ii+1));
end
%%
ydata_pad = [ydata(1),ydata,ydata(end)];
for jj=1:length(ydata) 
    ydata_smooth(jj) = mean(ydata_pad(jj:jj+2));
end
stage_ave = stage_total/runtimes_sim;
figure
plot_setting
plot(-(0.5:1:39.5),ydata);hold on
% plot(-(0.5:1:39.5),ydata);hold on
plot(-(0.5:1:16.5),ones(1,17)*0.9,'k--');hold on
plot(-(23.5:1:39.5),ones(1,17)*0.1,'k--');hold on
xlabel('Interferers Power (dBm)')
ylabel('Probability of being tagged')
legend('Sim.','Designed Bound')
grid on
%% detection
% pdpfa_s1 = mean(tag_table_s1(:,:));
% pdpfa_s2 = mean(tag_table_s2(:,:));
% pdpfa_all = mean(tag_table_s1(:,:)+tag_table_s2(:,:));
%% figure
% figure
% plot_setting;
% plot(pdpfa_s1,'b');hold all
% plot(pdpfa_s2,'r');hold all
% plot(pdpfa_all,'c');hold all
% xlabel('index')
% ylabel('Probability of Detection')

%%
% for ii=1:blk_num
%     sig_pow_hat(ii) = sum(sig_cache(:,ii).^2)/P;
% end
% figure
% plot(10*log10(sig_pow_hat)+3);hold on
% plot(10*log10(sig_pow))
%% theoretical curve for R = 32 , N = 31
clf
M=15:120;
pur_anal1 = 10*log10(((1/31+0.005*0.7165)*(sqrt(M/2)+qfuncinv(0.1))./(sqrt(M/2)-qfuncinv(0.1))-1/31)/0.7165)+3;
figure(88)
plot(M*30,pur_anal1,'b');hold on
plot([15 30 60 90 120]*30,[-8,-9.9,-11.75,-12.5,-13],'bo--')
% theoretical curve for R = 32 , N = 63
M=15:120;
pur_anal2 = 10*log10(((1/63+0.005*0.4319)*(sqrt(M/2)+qfuncinv(0.1))./(sqrt(M/2)-qfuncinv(0.1))-1/63)/0.4319)+3;
figure(88)
plot(M*60,pur_anal2,'r');hold on
plot([15 30 45 60 90 120]*60,[-8.75,-10.25,-11.3,-12.3,-13.2,-13.7],'ro--')

figure(88)
plot([30 60 90 120]*30,[-9.9,-13.02,-15.7,-16],'go--')
legend('Single-Stage, N = 31 (Anly.)','Single-Stage, N = 31 (Sim.)','Single-Stage, N = 63 (Anly.)','Single-Stage, N = 63 (Sim.)','Multi-Stage, N = 31 (Sim.)')
xlim([30,120]*30)
xlabel('Number of Frames (M)')
ylabel('P_u_r (dBm)')