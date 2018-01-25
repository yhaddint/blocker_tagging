
%
%  06/20/2014
%  We keep M = 16 and N =10, blocker_pow_pool is not a matrix 

clear;clc;clf;close all
warning off

trueCal=0;

symbols=5e4;
upsam=8;


N=5;
M=40;

%  tune power of each blockers
blocker_num=3;
int_power_pool=[1 5 10 20 50 100]/10';

blocker_pow_pool=zeros(length(int_power_pool),blocker_num);
blocker_pow_pool(:,1)=10;
blocker_pow_pool(:,2)=int_power_pool;
blocker_pow_pool(:,3)=int_power_pool;

% determined whether to apply estimation of power and detection
do_est=0;
do_det=1;

% number of data accumulated to watch histogram.
height=625;

% initialization of statistics store vector for int_power_pool
rslt_matrix_h0=zeros(length(int_power_pool),height);
rslt_matrix_h1=zeros(length(int_power_pool),height);

mark=zeros(1,length(int_power_pool));   

L=4e5;


for int_power_index=1:length(int_power_pool)
    
    %  QPSK Blocker generator with power scaled    
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
    blocker(ii,:)=real(sig_bb(1:L,1))./mean(real(sig_bb(:,1)).^2)*blocker_pow_pool(int_power_index,ii);
    end


    %  calibration seq, we use shifted version to approximate
    if trueCal
    cal0=PNgenerator_v1(L+2);
    cal(1,:)=cal0(1:end-2);
    cal(2,:)=cal0(2:end-1);
    cal(3,:)=cal0(3:end);
    end
    
    %  Another version of cal
    if trueCal==0
    cal0=PNgenerator_v1(M*N*4);
    cal(1,:)=kron(ones(2000,1),cal0(M*N*3+1:end));
    cal(2,:)=kron(ones(2000,1),cal0(M*N+1:2*M*N));
    cal(3,:)=kron(ones(2000,1),cal0(M*N*2+1:3*M*N));
    
    end

    %  calculate statistics
    stat_h0=0.5*(cal(1,:)*2+(blocker(2,:).^2+blocker(3,:).^2).*cal(1,:))+blocker(3,:).*cal(1,:).*cal(3,:)+blocker(2,:).*cal(1,:).*cal(2,:);
    stat_h1=0.5*(cal(1,:)+blocker(1,:).^2.*cal(1,:))+blocker(1,:)+stat_h0;
    

    for ii=1:fix(L/N/2)
        rslt_h0(ii)=abs(sum(stat_h0((ii-1)*N+1:ii*N))/N);
        rslt_h1(ii)=abs(sum(stat_h1((ii-1)*N+1:ii*N))/N);
    end

    
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
legend('int = 1','int = 5','int = 10','int = 20','int = 50','int = 100')

end




%% estimation of blocker power
if do_est

end