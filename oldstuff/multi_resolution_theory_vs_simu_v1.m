%  12/17/2015
%  compare simulation results using proposed method with simulation using
%  approximate noise
%  

clear;clc;%clf;close all

blk_num = 3;
M = 1e2;
N = 32;
SampleNumberAve = 128;

stat_num = 2e2;
P = M*N;
L = N*stat_num;

sig_pow_dB = [-9,-6,4];
sig_pow = 10.^(sig_pow_dB/10);
rand('seed',4);
for ii = 1:blk_num
%     sig_r(:,ii) = kron(randi(2,M*stat_num,1)*2-3,ones(N,1));
%     sig_i(:,ii) = kron(randi(2,M*stat_num,1)*2-3,ones(N,1));
%     sig(:,ii) = (sig_r(:,ii)+1j*sig_i(:,ii))/sqrt(2);
    upsam = SampleNumberAve;
    symbols = fix(L*2/upsam)*5;
    
    clearvars data temp_data
    hmod = modem.qammod('M', 4, 'InputType', 'integer');
    hdesign  = fdesign.pulseshaping(upsam,'Square Root Raised Cosine');
    hpulse = design(hdesign);
    data = randi(4,symbols,1)-1;
    data = modulate(hmod, data);
    data = upsample(data,upsam);
    temp_data = conv(data,hpulse.Numerator);
    sig0(:,ii) = temp_data./sqrt(sum(temp_data.*conj(temp_data))/length(temp_data));
    sig(:,ii) = sig0(2000:end-2000,ii);
end
freq_mask0 = transpose(exp(1j*2*2*pi*(1:length(sig(:,1)))/N));
sig(:,1) = sig(:,1).*freq_mask0;

freq_mask1 = transpose(exp(1j*6*2*pi*(1:length(sig(:,1)))/N));
sig(:,2) = sig(:,2).*freq_mask1;

freq_mask2 = transpose(exp(1j*2*2*pi*(1:length(sig(:,1)))/N));
sig(:,3) = sig(:,3).*freq_mask2;

best_bin_ref = [19,7];
%% a simple test with 3 interferers.
%  First two are in the same band but different subband, the third is in
%  another band
fft_level = 1;
Nn = fft_level*N;
runtimes = 32*24*20;

bins_bands = zeros(runtimes,Nn/2,2);
for runindex = 1:runtimes
    phi = exp(1j*2*pi*rand(1,3));
    %phi = ones(1,3);
    sig_sum = zeros(P,1);
    CAL = randi(2,P,3)*2-3;
    sig_cache = zeros(P,3);
    start_point = randi(length(sig(:,1))-P,1,3);
    for ii=1:3
        sig_cache(:,ii) = sig(start_point(ii):start_point(ii)+P-1,ii);
        if ii ~= 3
            y(:,ii) = 2*real(sig_cache(:,ii).*phi(ii)).*CAL(:,1)*sqrt(sig_pow(ii));
        else
            y(:,ii) = 2*real(sig_cache(:,ii).*phi(ii)).*CAL(:,2)*sqrt(sig_pow(ii));
        end
        sig_sum = sum(y,2);
    end

    for kk=1:3
        sig_sum_corr = sig_sum.*CAL(:,kk);
        periodo_gram = fft(sig_sum_corr(1:Nn))/Nn;
        periodo_square = (abs(periodo_gram).^2)';
        bins_bands(runindex,:,kk) = periodo_square(1:Nn/2);
    end
end
figure(1)
for kk=1:3
    temp = squeeze(bins_bands(:,:,kk));
    ydata = mean(abs(temp),1);
    subplot(3,1,kk)
    plot(10*log10(ydata));hold on
    grid on
end
%\
% figure(99)
% for kk=1:2
%     for ii=1:16
%         bins_bands0(:,ii,kk) = bins_bands(:,(ii-1)*fft_level+1,kk);
%     end
%     temp = squeeze(bins_bands0(:,:,kk));
%     ydata = mean(abs(temp),1);
%     subplot(2,1,kk)
%     plot(10*log10(ydata));hold on
%     grid on
% end
%
ave_num_pool = [8,16,24,32,40,48,56,64]/fft_level;
pd = zeros(1,length(ave_num_pool));
for ave_index = 1:length(ave_num_pool)
    ave_num = ave_num_pool(ave_index);
    pd_counter = 0;
    for ii = 1:runtimes/ave_num
        for kk=1:3
            A(:,kk) = mean(bins_bands((ii-1)*ave_num+1:ii*ave_num,:,kk));
        end
        [sortval,sortindex] = sort((A(:)),'descend');
        bestbin(ii,:)=sortindex(1:2);
        if sort(sortindex(1:2))'==sort(best_bin_ref(1:2))
            pd_counter = pd_counter+1;
        end
    end
    pd(ave_index) = pd_counter/(runtimes/ave_num);
end
figure(2)
plot(ave_num_pool*Nn,pd,'x-');
hold on
grid on
xlabel('Number of Samples')
ylabel('Ident. Prob.')
ylim([0,1])


%% use noise with equivalent power to approximate
runtimes = 32*24*50;
leakage_table=ones(1,16)*(10.^(-33.5/10));
leakage_table(1:4) = 10.^([0.035,-22.73,-28.42,-31.51]/10);

bins_bands_aprx = zeros(runtimes,16,2);
bins_bands_aprx0 = (randn(runtimes,16,2)+1j*randn(runtimes,16,2))/sqrt(2);
% initialization of power table
s_power_table = zeros(16,2);
intra_power_table = zeros(16,2);
inter_power_table = zeros(16,2);
% self power
s_power_table(3,1) = sig_pow(1);
s_power_table(7,1) = sig_pow(2);
s_power_table(3,2) = sig_pow(3);
% interference to other bands (n_inter)
inter_power_table(:,1) = 2*sum(s_power_table(:,2),1)/N;
inter_power_table(:,2) = 2*sum(s_power_table(:,1),1)/N;
% interference within the band (n_intra)
for kk=1:2
    for ii=1:16
        for jj=1:16
            intra_power_table(ii,kk) = intra_power_table(ii,kk)+leakage_table(1+abs(ii-jj))*s_power_table(jj,kk);
        end
    end
end
power_table = inter_power_table+intra_power_table;
%
figure(1)
subplot(211)
plot(10*log10(power_table(:,1)),'r');hold on
subplot(212)
plot(10*log10(power_table(:,2)),'r');hold on
% scale the noise
for kk=1:2
    for ii=1:16
        if kk==1
            if ii<=16
                bins_bands_aprx(:,ii,kk) = abs(bins_bands_aprx0(:,ii,kk)).^2.*power_table(ii,kk);
            end
        else
            if ii<=16
                bins_bands_aprx(:,ii,kk) = abs(bins_bands_aprx0(:,ii,kk)).^2.*power_table(ii,kk);
            end
        end
    end
end
%
ave_num_pool = [8,16,24,32];
pd = zeros(1,length(ave_num_pool));
    
for ave_index = 1:length(ave_num_pool)
    ave_num = ave_num_pool(ave_index);
    pd_counter = 0;
    for ii = 1:runtimes/ave_num
        A_aprx(:,1) = mean(bins_bands_aprx((ii-1)*ave_num+1:ii*ave_num,:,1));
        A_aprx(:,2) = mean(bins_bands_aprx((ii-1)*ave_num+1:ii*ave_num,:,2));
        [sortval,sortindex] = sort((A_aprx(:)),'descend');
        bestbin_aprx(ii,:)=sortindex(1:2);
        if sort(sortindex(1:2))'==sort(best_bin_ref(1:2))
            pd_counter = pd_counter+1;
        end
    end
    pd(ave_index) = pd_counter/(runtimes/ave_num);
end
figure(2)
plot(ave_num_pool*N,pd,'o--');
hold on
grid on
xlabel('Number of Samples')
ylabel('Ident. Prob.')
ylim([0,1])
%% theoretical approximation
n=31;
N=32;
% data0 = zeros(1e3,n);
% ave_num_pool = 1;
ave_num_pool = 8:64;
% for ave_index = 1:length(ave_num_pool)
%     ave_num = ave_num_pool(ave_index);
%     for ii=1:1e3
%         data0(ii,:) = 10^(-1.16)*chi2rnd(ave_num*2,1,n)/ave_num/2;
%         chi2(ii) = max(data0(ii,:));
%     end
%     mu(ave_index) = mean(chi2);
%     sigma(ave_index) = var(chi2);
% end
% mu_aprx = 10^(-11.6/10).*(1+2.05.*sqrt(1./[8,16,24,32]));
mu_aprx = 10^(-11.2494/10).*(1+2.23.*sqrt(1./ave_num_pool));

% figure(98)
% xdata = 0.05:0.005:0.15;
% plot(xdata,normcdf(xdata,mu(end),sqrt(sigma(end))),'b');hold on
% [a,b] = ecdf(chi2);
% figure(98)
% plot(b,a,'r-')
% figure(99)
% hist(chi2,50)
%
mu_s = (1+0.2)*2/32+0.2;
% mu_s = (1+0.1)*2/32+0.1;
% mu_s = 0.1419;
% mu_s = (1+0.05)*2/32+0.05;
figure(2)
% plot(ave_num_pool*Nn,(1-qfunc((0.1296-0.0631)./sqrt((0.1296^2+0.0631^2)./ave_num_pool))))
for ave_index = 1:length(ave_num_pool)
    ave_num = ave_num_pool(ave_index);
    temp0(ave_index) = qfunc((mu_aprx(ave_index)-mu_s)/(sqrt(mu_s^2/ave_num)));
end
plot(ave_num_pool*N,temp0)
%%
% for ii = 1:runtimes/8
%     temp(ii,:) = mean(bins_bands_aprx((ii-1)*8+1:ii*8,:,2));
% end


%%
 %clear;clc;clf

%%
ave_num=32;
chi = zeros(runtimes/ave_num,16);
for runindex=1:runtimes/ave_num
    chi(runindex,:) = mean(squeeze(bins_bands_aprx((runindex-1)*ave_num+1:runindex*ave_num,:,1)));
end
for runindex=1:runtimes/ave_num
    record(runindex) = max(chi(runindex,[1:6,8:16]));
end
mean(record)
var(record)
