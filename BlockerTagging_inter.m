%
%  06/21/2014
%  watch statistics of each stage
%  06/22/2014
%  actually, this method doesn't help to surpress interference a lot.

clear;clc;clf;close all
warning off

symbols=5e4;
upsam=500;


N=200;
M=20;

%  tune power of each blockers
blocker_num=3;
int_power_pool=[5 10 20 50 100 200]/100';
%int_power_pool=0.1';

blocker_pow_pool=zeros(length(int_power_pool),blocker_num);
blocker_pow_pool(:,1)=1;
blocker_pow_pool(:,2)=int_power_pool;
blocker_pow_pool(:,3)=int_power_pool;

% determined whether to apply estimation of power and detection
do_est=0;
do_det=1;

% number of data accumulated to watch histogram.
height=200;

% initialization of statistics store vector for int_power_pool
rslt_matrix_h0=zeros(length(int_power_pool),height);
rslt_matrix_h1=zeros(length(int_power_pool),height);

mark=zeros(1,length(int_power_pool));   

L=4e5;

for ii=1:3
    clear data
    hmod = modem.pskmod('M', 4, 'InputType', 'integer');
    hdesign  = fdesign.pulseshaping(upsam,'Square Root Raised Cosine');
    hpulse = design(hdesign);
    data = randi(4,symbols,1)-1;
    data = modulate(hmod, data);
    data = upsample(data,upsam);
    sig_bb(:,ii) = conv(data,hpulse.Numerator);
    sig_bb(:,ii)=sig_bb(:,ii)/(sqrt(mean(abs(sig_bb(:,ii)).^2)));
end

for int_power_index=1:length(int_power_pool)
    
    %  QPSK Blocker generator with power scaled    
    for ii=1:3
        blocker(ii,:)=real(sig_bb(1:L,1))./mean(real(sig_bb(:,1)).^2)*blocker_pow_pool(int_power_index,ii);
    end


    cal0=PNgenerator_v1(M*N*4);
    cal(1,:)=kron(ones(1e3*M/N,1),cal0(M*N*3+1:end));
    cal(2,:)=kron(ones(1e3*M/N,1),cal0(M*N+1:2*M*N));
    cal(3,:)=kron(ones(1e3*M/N,1),cal0(M*N*2+1:3*M*N));
    

    %  calculate statistics
    stat_h0=0.5*(cal(1,:)*2+(blocker(2,:).^2+blocker(3,:).^2).*cal(1,:))+...
        blocker(3,:).*cal(1,:).*cal(3,:)+blocker(2,:).*cal(1,:).*cal(2,:);
    stat_h1=0.5*(cal(1,:)+blocker(1,:).^2.*cal(1,:))+blocker(1,:)+stat_h0;
    temp1=0.5*(cal(1,:)+blocker(1,:).^2.*cal(1,:));
    temp2=0.5*(cal(1,:)*2+(blocker(2,:).^2+blocker(3,:).^2).*cal(1,:));
    temp3=blocker(3,:).*cal(1,:).*cal(3,:)+blocker(2,:).*cal(1,:).*cal(2,:);
    temp4=blocker(1,:);
    
    for ii=1:fix(L/N/2)
        temp1_sum(ii)=sum(temp1((ii-1)*N+1:ii*N)/N);
        temp2_sum(ii)=sum(temp2((ii-1)*N+1:ii*N)/N);
        temp3_sum(ii)=sum(temp3((ii-1)*N+1:ii*N)/N);
        temp4_sum(ii)=sum(temp4((ii-1)*N+1:ii*N)/N);
        rslt_h0(ii)=abs(sum(stat_h0((ii-1)*N+1:ii*N)/N));
        rslt_h1(ii)=abs(sum(stat_h1((ii-1)*N+1:ii*N)/N));
    end
    
%     %%
%     figure
%     subplot(2,2,1)
%     stem(temp1(100:200));
%     subplot(2,2,2)
%     stem(temp2(100:200));
%     subplot(2,2,3)
%     stem(temp3(100:200));
%     subplot(2,2,4)
%     stem(temp4(100:200));
    %%
    
    S_sta_h0(:,int_power_index)=rslt_h0;
    S_sta_h1(:,int_power_index)=rslt_h1;
    S_mean_h0(int_power_index)=mean(S_sta_h0(:,int_power_index));
    S_mean_h1(int_power_index)=mean(S_sta_h1(:,int_power_index));
    S_var_h0(int_power_index)=var(S_sta_h0(:,int_power_index));
    S_var_h1(int_power_index)=var(S_sta_h1(:,int_power_index));
    
    l_need=fix(length(rslt_h0)/M);
    if M==1
        temp_h0=rslt_h0;
        temp_h1=rslt_h1;
    else
        temp_h0=mean(reshape(rslt_h0(1:M*l_need),M,l_need));
        temp_h1=mean(reshape(rslt_h1(1:M*l_need),M,l_need));
    end
    
    
    if mark(int_power_index)<height
    % it means we should update new data into it
        if mark(int_power_index)+length(temp_h0)<=height
            rslt_matrix_h0(int_power_index,mark(int_power_index)+1:mark(int_power_index)+length(temp_h0))=temp_h0;
            rslt_matrix_h1(int_power_index,mark(int_power_index)+1:mark(int_power_index)+length(temp_h1))=temp_h1;
            mark(int_power_index)=mark(int_power_index)+length(temp_h0);
        else
            rslt_matrix_h0(int_power_index,mark(int_power_index)+1:height)=temp_h0(1:height-mark(int_power_index));
            rslt_matrix_h1(int_power_index,mark(int_power_index)+1:height)=temp_h1(1:height-mark(int_power_index));
            mark(int_power_index)=height;
        end
    end
    
    
end


%% detection and plot
% use do_det to control this section
if do_det
color=['b  ';'r  ';'c  ';'k  ';'r--';'c--'];
    
for int_power_index=1:length(int_power_pool)
    
    
    temp_h0=squeeze(rslt_matrix_h0(int_power_index,:));
    temp_h1=squeeze(rslt_matrix_h1(int_power_index,:));
    th_pool=linspace(0,max(temp_h0),4000);

    for th_index=1:length(th_pool)
        th=th_pool(th_index);
        prob_fa(th_index)=sum(temp_h0>th)/mark(int_power_index);
        prob_d(th_index)=sum(temp_h1>th)/mark(int_power_index);
    end

    figure(1)
    plot(prob_fa,prob_d,color(int_power_index,:),'linewidth',2);hold on
    grid on
    
    clearvars temp_h0 temp_h1
    
    figure(2)
    subplot(2,3,int_power_index)
    temp_h0=squeeze(rslt_matrix_h0(int_power_index,:))';
    temp_h1=squeeze(rslt_matrix_h1(int_power_index,:))';
    hist_axis=linspace(0,fix(max(temp_h1)),10);
    hist([temp_h0 temp_h1],hist_axis);hold on
    title(['Int = ' num2str(int_power_pool(int_power_index))]);

end

figure(1)
title('ROC curve');
legend('int = 0.0625','int = 0.125','int = 0.25','int = 0.5','int = 1','int = 2')

end




%% hist of stat S
figure(10)
for ii=1:length(int_power_pool)
subplot(2,3,ii)
hist(S_sta_h0(:,ii));hold on
end

figure(11)
for ii=1:length(int_power_pool)
subplot(2,3,ii)
hist(S_sta_h1(:,ii));hold on
end
