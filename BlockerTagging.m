
% 04/18/2014
% Different m and n are swept to watch detection performance only single blocker is considered, and thus results could be trivial noise is neglected here since compared with 0 dBm blocker and several dB weaker PN signal, noise floor is too weak.

% 04/25/2014 
% Since larger Blocker/PN ratio leads to better performance, why not lower PN power to noise floor? 

clear;clc;clf;close all
warning off
symbols=5e4;
upsam=8;
Npool=[5 10 20 40];
Mpool=[1 2 4 8];
% Npool=[100];
% Mpool=[4];

noisepower=1;
sig_pow=100;

needmoredata=1;
height=5e3;

rslt_matrix_h0=zeros(length(Npool),length(Mpool),height);
rslt_matrix_h1=zeros(length(Npool),length(Mpool),height);

mark=zeros(length(Npool),length(Mpool));

for runtime=1
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
    sig_use=sig_normal*sig_pow;
    cal=PNgenerator_v1(L);
    
    noise_0=randn(1,L)*noisepower;
    noise_1=randn(1,L)*noisepower;
    
    sig_h0=noise_0';
    sig_h1=sig_use+noise_1';
    
    stat_h0=0.5*(cal+sig_h0.^2.*cal)+sig_h0;
    stat_h1=0.5*(cal+sig_h1.^2.*cal)+sig_h1;
    
    for nn=1:length(Npool)
    for mm=1:length(Mpool)
        clc
        display(['Run number',num2str(runtime)]);
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
    
    
    if mark(nn,mm)<height
    % it means we should update new data into it
        if mark(nn,mm)+length(temp_h0)<=height
            rslt_matrix_h0(nn,mm,mark(nn,mm)+1:mark(nn,mm)+length(temp_h0))=temp_h0;
            rslt_matrix_h1(nn,mm,mark(nn,mm)+1:mark(nn,mm)+length(temp_h1))=temp_h1;
            mark(nn,mm)=mark(nn,mm)+length(temp_h0);
        else
            rslt_matrix_h0(nn,mm,mark(nn,mm)+1:height)=temp_h0(1:height-mark(nn,mm));
            rslt_matrix_h1(nn,mm,mark(nn,mm)+1:height)=temp_h1(1:height-mark(nn,mm));
            mark(nn,mm)=height;
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
%% plot
color=['b  ';'r  ';'c  ';'k  ';'r--';'c--'];
th_pool=0:0.1:50;
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
%% plot hist of h0
figure;
for mm=1:4
subplot(2,2,mm);
hist(squeeze(rslt_matrix_h0(mm,1,:)),10);
title(['N = 5, M=' num2str(Mpool(mm))]);
end


