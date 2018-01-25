
% 04/25/2014
% Keep M=1 and tune N from small to infinite to see whether statistic will
% approach blocker power. Use different blocker power to watch clearly.


clear;clc;clf;close all
warning off

symbols=5e4;
upsam=8;

% tune parameter M, N and relative power of blockers v.s. PN sequences
Npool=[5 10 20 40 80 160 320 640 1280];
Mpool=[1 2 4 8];
sig_pow_pool=[1 2 4 8 16];
noisepower=0;

% determined whether to apply estimation of power and detection
do_est=1;
do_det=0;

% number of data accumulated to watch histogram.
height=5e3;

% initialization of data store matrix for M, N pair
rslt_matrix_h0=zeros(length(Npool),length(sig_pow_pool),height);
rslt_matrix_h1=zeros(length(Npool),length(sig_pow_pool),height);

mark=zeros(length(Npool),length(sig_pow_pool));   


for sig_pow_index=1:length(sig_pow_pool)

for runtime=1:1
    for i=1:1
    clear data
    hmod = modem.pskmod('M', 4, 'InputType', 'integer');
    hdesign  = fdesign.pulseshaping(upsam,'Square Root Raised Cosine');
    hpulse = design(hdesign);
    data = randi(4,symbols,1)-1;
    data = modulate(hmod, data);
    data = upsample(data,upsam);
    data = conv(data,hpulse.Numerator);

    sig_bb(:,i)=data;
    sig_bb(:,i)=sig_bb(:,i)/(sqrt(mean(abs(sig_bb(:,i)).^2)));
    end


    %% signal processing
    L=4e5;
    sig_normal=real(sig_bb(1:L,1))./mean(real(sig_bb(:,1)).^2);
    cal=PNgenerator_v1(L);
    noise_0=randn(1,L)*noisepower;
    noise_1=randn(1,L)*noisepower;

    sig_pow=sig_pow_pool(sig_pow_index);
    sig_use=sig_normal*sig_pow;
    
    sig_h0=noise_0';
    sig_h1=sig_use+noise_1';
    
    stat_h0=0.5*(cal+sig_h0.^2.*cal)+sig_h0;
    stat_h1=0.5*(cal+sig_h1.^2.*cal)+sig_h1;

    
    for nn=1:length(Npool)
    for mm=1%:length(Mpool)
        clc
        display(['Run number',num2str(runtime)]);
        display(['Blockers vs PN',num2str(sig_pow)]);
        display(['N = ',num2str(Npool(nn)) ', M = ',num2str(Mpool(mm))]);
    clear rslt_cancel rslt_noncancel rslt_PN ii i
    
    N=Npool(nn);
    M=Mpool(mm);

    for ii=1:fix(L/N/2)
        rslt_h0(ii)=abs(mean(stat_h0((ii-1)*N+1:ii*N)));
        rslt_h1(ii)=abs(mean(stat_h1((ii-1)*N+1:ii*N)));
    end

    
    l_need=fix(length(rslt_h0)/M);
    if M==1
        temp_h0=rslt_h0;
        temp_h1=rslt_h1;
    else
        temp_h0=mean(reshape(rslt_h0(1:M*l_need),M,l_need));
        temp_h1=mean(reshape(rslt_h1(1:M*l_need),M,l_need));
    end
    
    
    if mark(nn,sig_pow_index)<height
    % it means we should update new data into it
        if mark(nn,sig_pow_index)+length(temp_h0)<=height
            rslt_matrix_h0(nn,sig_pow_index,mark(nn,sig_pow_index)+1:mark(nn,sig_pow_index)+length(temp_h0))=temp_h0;
            rslt_matrix_h1(nn,sig_pow_index,mark(nn,sig_pow_index)+1:mark(nn,sig_pow_index)+length(temp_h1))=temp_h1;
            mark(nn,sig_pow_index)=mark(nn,sig_pow_index)+length(temp_h0);
        else
            rslt_matrix_h0(nn,sig_pow_index,mark(nn,sig_pow_index)+1:height)=temp_h0(1:height-mark(nn,sig_pow_index));
            rslt_matrix_h1(nn,sig_pow_index,mark(nn,sig_pow_index)+1:height)=temp_h1(1:height-mark(nn,sig_pow_index));
            mark(nn,sig_pow_index)=height;
        end
    end
    % figure
    % hist(rslt2,20)
    % title('sequencially pick')


    
    % clear ii
    % for ii=1:fix(N/M)
    %     rslt4(ii)=mean(cal((ii-1)*M+1:ii*M).^3);
    % end

    % %% sequencially pick illustration
    % figure
    % plot(real(sig_bb(1:500,1)),'c','linewidth',2);
    % hold
    % plot(100:150,real(sig_bb(100:150,1)),'kx');
    % grid
    % 
    % 
    % %% randomly pick illustration
    % index=randi(500,1,50);
    % figure
    % plot(real(sig_bb(1:500,1)),'c','linewidth',2);
    % hold
    % plot(sort(index),real(sig_bb(sort(index))),'kx');
    % grid


    %[H,gaussfit(mm,nn,:)]=chi2gof(rslt_matrix(mm,nn,:));
    end
    end
end
    % Before next signal power, accumulate needed data
    for nn=1:length(Npool)
        h1_stat_mean(sig_pow_index,nn)=mean(squeeze(rslt_matrix_h1(nn,sig_pow_index,:)));
    end
end

%% detection and plot
% use do_det to control this section
if do_det
color=['b  ';'r  ';'c  ';'k  ';'r--';'c--'];
th_pool=0:0.01:5;
for nn=1:length(Npool)
    subplot(2,2,nn)
    for mm=1:length(Mpool)
        temp_h0=squeeze(rslt_matrix_h0(nn,mm,:));
        temp_h1=squeeze(rslt_matrix_h1(nn,mm,:));
        for th_index=1:length(th_pool)
            th=th_pool(th_index);
            prob_fa(th_index)=sum(temp_h0>th)/mark(nn,mm);
            prob_d(th_index)=sum(temp_h1>th)/mark(nn,mm);
        end
        figure(1)
        plot(prob_fa,prob_d,color(mm,:),'linewidth',2);hold on
        grid on
        
        if Npool(nn)*Mpool(mm)==40
            figure(99)
            plot(prob_fa,prob_d,color(mm,:),'linewidth',2);hold on
        end
    end
    figure(1)
    title(['N = ' num2str(Npool(nn))])
    legend('M = 1','M = 2','M = 4','M = 8')

end

figure(99)
title('With Same Sample Number')
legend('N=5, M=8','N=10, M=4','N=20, M=2','N=40, M=1');
grid on

% temp=fix(10*log10(sig_pow));
% saveas(figure(1),['ratio' num2str(temp) 'sweep.eps'],'epsc');
% saveas(figure(99),['ratio' num2str(temp) 'sample.eps'],'epsc');

% %% hist of noncancelling term
% figure
% for ii=1:4
% subplot(2,2,ii)
% hist(squeeze(rslt_matrix(ii,1,:)))
% title(['N=',num2str(Npool(ii))]);
% end

% plot hist of h0
figure;
for mm=1:4
subplot(2,2,mm);
hist(squeeze(rslt_matrix_h0(mm,1,:)),10);
title(['N = 5, M=' num2str(Mpool(mm))]);
end

end

%% estimation of blocker power
if do_est
    AP=1;
    color=['b-o';'r-o';'c-o';'k-o';'m-o'];
    color_AP=['b--';'r--';'c--';'k--';'m--'];
    figure(99)
    for sig_pow_index=1:length(sig_pow_pool)
        semilogy(Npool,h1_stat_mean(sig_pow_index,:),color(sig_pow_index,:),'linewidth',2);
        hold on
    end
    grid on
    legend('S/PN = 0dB','S/PN = 3dB','S/PN = 6dB','S/PN = 9dB','S/PN = 12dB')
    xlabel('N')
    ylabel('statistic value')
    if AP
        figure(99)
        for sig_pow_index=1:length(sig_pow_pool)
            semilogy(Npool,ones(1,length(Npool))*sig_pow_pool(sig_pow_index),color_AP(sig_pow_index,:),'linewidth',2);
            hold on
        end
    end
end
