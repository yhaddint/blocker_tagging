%  03/31/2016
%  interference level in the conventional method
clear;clc;%clf;close all
rand('seed',3)
blk_num = 4;
M = 1e2;
N = 63;
SampleNumberAve = 32;
stat_num = 2e2;
P = M*N;
L = M*SampleNumberAve*stat_num;


other_pow = 1e-4;
sig_pow0 = ones(blk_num,1)*other_pow;
sig_pow1 = ones(blk_num,1)*other_pow;

sig_pow1([1,4]) = [0.31,1];
sig_pow0([1,4]) = [0,1];

for ii = 1:blk_num
%     sig_r(:,ii) = kron(randi(2,M*stat_num,1)*2-3,ones(SampleNumberAve,1));
%     sig_i(:,ii) = kron(randi(2,M*stat_num,1)*2-3,ones(SampleNumberAve,1));
%     sig(:,ii) = (sig_r(:,ii)+sig_i(:,ii))/sqrt(2);
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
    %sig_i(:,ii) = zeros(size(sig_r(:,ii)));
    sig(:,ii) = (sig_r(:,ii)+sig_i(:,ii))/sqrt(2);

end
%% calibration sequences
cal0=PNgenerator_v5(N,N,1);
CAL_noperm=LowerRate_v2(cal0,P);
CAL=zeros(P,N);
for mm=1:M
    CAL((mm-1)*N+1:mm*N,:)=CAL_noperm((mm-1)*N+1:mm*N,randperm(N));
end
%% tagging
runtimes_sim = 1;
result_noproc_sim_perm=zeros(runtimes_sim,M,blk_num);

for runindex=1:runtimes_sim
    if mod(runindex,100)==0
        runindex/100
        CAL=zeros(P,N);
        for mm=1:M
            CAL((mm-1)*N+1:mm*N,:)=CAL_noperm((mm-1)*N+1:mm*N,randperm(N));
        end
    end
    %CAL=CAL(:,randperm(N));
    start_point = randi(L-P-1,blk_num,M);
    for ii=1:blk_num
        sig_cache(:,ii) = sig(start_point(ii,mm):start_point(ii,mm)+P-1,randi(blk_num));
    end
    sig_mix0 = zeros(P,1);
    sig_mix1 = zeros(P,1);

    for ii=1:blk_num
        for jj=1:blk_num
            if abs(ii-jj)<=10
            sig_mix0 = sig_mix0 + sig_cache(:,ii).*CAL(:,jj).*sin(2*pi*rand+(ii-jj)/32*2*pi*(1:P)')*sqrt(sig_pow0(ii));
            sig_mix1 = sig_mix1 + sig_cache(:,ii).*CAL(:,jj).*sin(2*pi*rand+(ii-jj)/32*2*pi*(1:P)')*sqrt(sig_pow1(ii));

            end
        end
    end
    
    for ii=1:blk_num 
        
        ydata0 = abs(fft(sig_mix0(1:N*M).*CAL(1:N*M,ii)));
        ydata1 = abs(fft(sig_mix1(1:N*M).*CAL(1:N*M,ii)));
        l_ydata = length(ydata1)
        xdata = linspace(-1,1,l_ydata);
        figure;
        plot(xdata,[ydata1(l_ydata/2:end);ydata1(1:l_ydata/2-1)],'k');hold on
        plot(xdata,[ydata0(l_ydata/2:end);ydata0(1:l_ydata/2-1)],'g');
        DLPF = abs(xdata)<1/16;
        plot(xdata,800*DLPF,'k--');
        ylim([0,1000])
        xlabel('Frequency (\pi)')
        ylabel('relative amplitude')
    end
end


%% figure

figure
plot(abs(fft(sig_cache(:,1))));hold on
xlabel('Number of Frames (M)')
ylabel('Identification Probability')
%%

ydata = abs(fft(sig_mix0(1:N*M).*CAL(1:N*M,ii)));
l_ydata = length(ydata)
figure;
plot(ydata)
figure;
plot([ydata(l_ydata/2:end);ydata(1:l_ydata/2)]);  
