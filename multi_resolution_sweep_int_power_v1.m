%  12/23/2015
%  test how many samples are required for STPE method in finding strongest
%  blocker.
%  



%%
% fist we use theoretical formula to plot two curve of sample v.s. power of int. #3
%
% curve in case 1
clear;clc
Ps_dB = fliplr(-18:0.1:-7);
Ps = 10.^(Ps_dB./10);
P_tot = Ps+0.25+1;
% P_tot = Ps+0.25+1;

M1 = qfuncinv(0.1)^2*((Ps+P_tot/16).^2/2+(0.25+P_tot/16).^2/2)./(Ps-0.25).^2;
figure(88)
plot(Ps_dB,M1*32,'k-.');hold on
xlim([-18,-8])
%% curve in case 2

% mu_s = P_tot./16+0.25;
% 
% % M2 = ((2.23*(P_tot./16)-mu_s.*qfuncinv(0.9))./0.25).^2;
% M2 = ((2.78*(P_tot./16)-mu_s.*qfuncinv(0.9))./0.25).^2;
% figure(88)
% plot(Ps_dB,M2*32,'r--');hold on
%% curve in case 2 (better one)
mu_s = P_tot./16+0.25;
M=linspace(0,10,1e5);
for ii=1:length(mu_s)
    mu_t = P_tot(ii)/16.*(1+1./sqrt(M).*qfuncinv(5/(8*(8*16-2)+2)));
    sigma_t = (P_tot(ii)/16)^2*pi^2/12./M/log(16*8-2);
    curve = 1-qfunc(sqrt((M.*(mu_t-mu_s(ii)).^2)./(sigma_t+1./M*(mu_s(ii))^2)));
    [bestscore(ii),best(ii)] = min(abs(curve-0.9));
end
figure(88)
plot(Ps_dB,M(best)*32,'k--');hold on
% mu_aprx = 10^(-11.2494/10).*(1+2.23.*sqrt(1./ave_num_pool));
% mu_s = (1+0.2)*2/32+0.2;
% mu_s = (1+0.1)*2/32+0.1;
% mu_s = 0.1419;
% mu_s = (1+0.05)*2/32+0.05;
% figure(88)
% % plot(ave_num_pool*Nn,(1-qfunc((0.1296-0.0631)./sqrt((0.1296^2+0.0631^2)./ave_num_pool))))
% for ave_index = 1:length(ave_num_pool)
%     ave_num = ave_num_pool(ave_index);
%     temp0(ave_index) = 1-qfunc((mu_s-mu_aprx(ave_index))/(sqrt(mu_s^2/ave_num)));
% end
% plot(ave_num_pool*N,temp0)
%%
figure(2)
% plot(ave_num_pool*Nn,(1-qfunc(sqrt(ave_num_pool).*(0.0927)./sqrt((0.2465^2/2+0.3380^2/2)))))
% plot(ave_num_pool*Nn,(1-qfunc(sqrt(ave_num_pool).*(0.125)./sqrt((0.3359^2/2+0.2109^2/2)))))
% plot(ave_num_pool*Nn,(1-qfunc(sqrt(ave_num_pool).*(0.125)./sqrt((0.25^2+0.125^2)))))
plot(ave_num_pool*Nn,(1-qfunc(sqrt(ave_num_pool).*(0.25-0.1585)./sqrt((0.25^2+0.1585^2)))))
%% simulated value on required sample in order to reach 0.9 ident. prob.
xdata = [-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8];
% ydata = [126.1,127.9,128.4,132.6,134.3,140.9,148.1,168.5,200,238.7,556.2];
ydata = [91.3,91.3,92,93.6,94.7,98,105.3,115.1,127.3,167,476.7];
figure(88)
plot(xdata,ydata,'bo-');hold on


%%

blk_num = 3;
M = 1e2;
N = 32;
SampleNumberAve = 128;

stat_num = 2e2;
P = M*N;
L = N*stat_num;

sig_pow_dB = [-7,-6,0];
sig_pow = 10.^(sig_pow_dB/10);
rand('seed',3);
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
% a simple test with 3 interferers.
%  First two are in the same band but different subband, the third is in
%  another band
fft_level = 1;
Nn = fft_level*N;
runtimes = 32*24*20;

bins_bands = zeros(runtimes,Nn/2,4);
for runindex = 1:runtimes
%     phi = exp(1j*2*pi*rand(1,3));
    %phi = ones(1,3);
    sig_sum = zeros(P,1);
    CAL = randi(2,P,4)*2-3;
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

    for kk=1:4
        sig_sum_corr = sig_sum.*CAL(:,kk);
        periodo_gram = fft(sig_sum_corr(1:Nn))/Nn;
        periodo_square = (abs(periodo_gram).^2)';
        bins_bands(runindex,:,kk) = periodo_square(1:Nn/2);
    end
end
figure(1)
for kk=1:4
    temp = squeeze(bins_bands(:,:,kk));
    ydata = mean(abs(temp),1);
    subplot(4,1,kk)
    plot(10*log10(ydata));hold on
    grid on
end

%
ave_num_pool = 2:1:16;
pd = zeros(1,length(ave_num_pool));
for ave_index = 1:length(ave_num_pool)
    ave_num = ave_num_pool(ave_index);
    pd_counter = 0;
    for ii = 1:fix(runtimes/ave_num)
        for kk=1:4
            A(:,kk) = mean(bins_bands((ii-1)*ave_num+1:ii*ave_num,:,kk));
        end
        [sortval,sortindex] = sort((A(:)),'descend');
        bestbin(ii,:)=sortindex(1:2);
        if sort(sortindex(1:2))'==sort(best_bin_ref(1:2))
            pd_counter = pd_counter+1;
        end
    end
    pd(ave_index) = pd_counter/(fix(runtimes/ave_num));
end
%
figure(2)
plot(ave_num_pool*Nn,pd,'x-');
hold on
grid on
xlabel('Number of Samples')
ylabel('Ident. Prob.')
ylim([0,1])
