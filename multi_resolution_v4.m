%  12/01/2015
%  test multi-resolution idea
%  Pseudo-spectrum is used. 

clear;clc;clf;close all

blk_num = 31;
M = 1e2;
N = 32;

stat_num = 2e2;
P = M*N;
L = N*stat_num;

for ii = 1:blk_num
    sig_r(:,ii) = kron(randi(2,M*stat_num,1)*2-3,ones(N,1));
    sig_i(:,ii) = kron(randi(2,M*stat_num,1)*2-3,ones(N,1));
    sig(:,ii) = (sig_r(:,ii)+1j*sig_i(:,ii))/sqrt(2);
%     upsam = SampleNumberAve;
%     symbols = fix(L*2/upsam);
%     
%     clearvars data temp_data
%     hmod = modem.qammod('M', 4, 'InputType', 'integer');
%     hdesign  = fdesign.pulseshaping(upsam,'Square Root Raised Cosine');
%     hpulse = design(hdesign);
%     data = randi(4,symbols,1)-1;
%     data = modulate(hmod, data);
%     data = upsample(data,upsam);
%     temp_data = conv(data,hpulse.Numerator);
%     sig(:,ii) = temp_data./sqrt(sum(temp_data.*conj(temp_data))/length(temp_data));
    
end
sig_table = zeros(3,3);
sig_table(:,1) = [1;1;2];
sig_table(:,2) = [1;6;3];
sig_pow_dB = [-3,-3,-10];
sig_table(:,3) = 10.^(sig_pow_dB/10);

for ii=1:3
    freq_mask = transpose(exp(1j*(sig_table(ii,2)-1)*2*pi*(1:length(sig(:,ii)))/N));
    sig(:,ii) = sig(:,ii).*freq_mask*sqrt(sig_table(ii,3));
end
%% a simple test with 3 interferers.
%  First two are in the same band but different subband, the third is in
%  other bands

Nn = N;
L = N;
runtimes = 5e3;
bins_bands = zeros(runtimes,L/2,2);

fake_noise = (randn(L,runtimes)+1j*randn(L,runtimes)).*sqrt(1/N);

for runindex = 1:runtimes
    sig_sum = zeros(Nn,1);
    CAL = randi(2,Nn,3)*2-3;
    sig_cache = zeros(Nn,3);
    start_point = randi(P*M,1,3);
    for ii=1:3
        sig_cache(:,ii) = sig(start_point(ii):start_point(ii)+Nn-1,ii);
        sig_sum = sig_sum+2*real(sig_cache(:,ii)).*CAL(:,sig_table(ii,1));
    end
    
    % use fft-first method
    for ii=1:2
        sig_sum_corr = sig_sum.*CAL(:,ii);
        periodo_gram = fft(sig_sum_corr(1:Nn))/Nn;
%         periodo_gram = fft(sig_sum_corr(1:Nn))/Nn+fake_noise(:,runindex);
        periodo_square = (abs(periodo_gram).^2)';
        bins_bands_fft(runindex,:,ii) = [periodo_square(1),periodo_square(2:L/2)];
    end
    
    % use corr first method (slightly larger data in use)

end

% for runindex = 1:runtimes
%     sig_sum = zeros(Nn+L,1);
%     CAL = randi(2,Nn+L,3)*2-3;
%     sig_cache = zeros(Nn+L,3);
%     start_point = randi(P*M,1,3);
%     for ii=1:3
%         sig_cache(:,ii) = sig(start_point(ii):start_point(ii)+Nn+L-1,ii);
%         sig_sum = sig_sum+2*real(sig_cache(:,ii)).*CAL(:,sig_table(ii,1));
%     end
%     
%     for ii=1:2
%         sig_sum_corr = sig_sum.*CAL(:,ii);
%         sig_conv(1)=0;
%         for ll=2:L
%             sig_conv(ll) = sum(sig_sum_corr(1:Nn).*(sig_sum_corr((1:Nn)+ll-1)))/Nn;
%         end
%         bins_bands_uncom = fft(sig_conv)/L;
%         bins_bands_corr(runindex,:,ii) = ([bins_bands_uncom(1),bins_bands_uncom(2:L/2)+fliplr(bins_bands_uncom(L/2+2:L))]);
%     end
% end
%% watch the noise is gaussian
% figure
% [a,b] = ecdf(real(record_noise));
% plot(b,a)
% 
% figure
% [a,b] = ecdf(imag(record_noise));
% plot(b,a)

%%
ave_num_pool = [1,2,4,8,16,32];
pd = zeros(1,length(ave_num_pool));
for ave_index = 1:length(ave_num_pool)
    clearvars bestscore bestbin
    ave_num = ave_num_pool(ave_index);
    if ave_num == 1
        for ii=1:runtimes
            [bestscore(ii),bestbin(ii)] = max(bins_bands_fft(ii,:,2));
        end
    else
        for ii = 1:runtimes/ave_num
            [bestscore(ii),bestbin(ii)] = max(mean(bins_bands_fft((ii-1)*ave_num+1:ii*ave_num,:,2)));
        end
    end

    pd(ave_index) = sum(bestbin ==3)/(runtimes/ave_num);
end
figure(2)
plot(ave_num_pool*Nn,pd,'x-');
hold on
grid on
xlabel('Number of Samples')
ylabel('Ident. Prob.')
ylim([0,1])

% 
% pd = zeros(1,length(ave_num_pool));
% for ave_index = 1:length(ave_num_pool)
%     clearvars bestscore bestbin
%     ave_num = ave_num_pool(ave_index);
%     if ave_num == 1
%         for ii=1:runtimes
%             [bestscore(ii),bestbin(ii)] = max(bins_bands_corr(ii,:,2));
%         end
%     else
%         for ii = 1:runtimes/ave_num
%             [bestscore(ii),bestbin(ii)] = max(mean(bins_bands_corr((ii-1)*ave_num+1:ii*ave_num,:,2)));
%         end
%     end
% 
%     pd(ave_index) = sum(bestbin ==3)/(runtimes/ave_num);
% end
% figure(2)
% plot(ave_num_pool*Nn+L,pd,'x-');
% hold on
% grid on
% xlabel('Number of Samples')
% ylabel('Ident. Prob.')
% ylim([0,1])
%%
% figure
figure(1)
for ii=1:2
    temp = squeeze(bins_bands_fft(:,:,ii));
    ydata = mean(abs(temp),1);
    subplot(2,1,ii)
    plot(10*log10(ydata));hold on
    grid on
end