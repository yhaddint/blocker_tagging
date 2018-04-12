%  12/18/2015
%  when signals have (block) fading, proposed method is supposed to be more robust
%  beacuse it takes more symbol duration to estimate power
%  
% take look at CDF of residual power

clear;clc;%clf;close all

band_num = 8;
blk_num = 4;
M = 1e2;
N = 32;
SampleNumberAve = 32;

stat_num = 2e2;
P = M*N;
L = N*stat_num*M-P;

fading_flag = 1;
noise_power = 1e-3;
P_c = 10^(40/10);

rand('seed',4);
for ii = 1:blk_num
%     sig_r(:,ii) = kron(randi(2,M*stat_num,1)*2-3,ones(N,1));
%     sig_i(:,ii) = kron(randi(2,M*stat_num,1)*2-3,ones(N,1));
%     sig0(:,ii) = (sig_r(:,ii)+1j*sig_i(:,ii))/sqrt(2);
    upsam = SampleNumberAve;
    symbols = fix(L*2/upsam);
    
    clearvars data temp_data
    qam_level = 4;
    hmod = modem.qammod('M', qam_level, 'InputType', 'integer');
    hdesign  = fdesign.pulseshaping(upsam,'Square Root Raised Cosine');
    hpulse = design(hdesign);
    data = randi(qam_level,symbols,1)-1;
    data = modulate(hmod, data);
    data = upsample(data,upsam);
    temp_data = conv(data,hpulse.Numerator);
    sig0(:,ii) = temp_data./sqrt(sum(temp_data.*conj(temp_data))/length(temp_data));
end

%% Formulating multi-band signals
band_realization_num = 1000;
runtimes = 20;
F_pool = [1,2,3];
F_num = length(F_pool);
ave_num_pool = [8,16,48];
pd = zeros(length(ave_num_pool),F_num,band_realization_num);

b_num=zeros(1,10);
pow_res = zeros(runtimes,length(ave_num_pool),F_num,band_realization_num);
for band_realization_index=1:band_realization_num
    disp(band_realization_index/band_realization_num*100)
    % signal table (band | subband | power)
    sig_table = zeros(blk_num,3);
    % sig in which band
    sig_table(:,1) = randi(band_num,blk_num,1);
    % sig in which subband
    sig_table(:,2) = randi(N/2-1,blk_num,1);
    % power of sig
    sig_pow_dB = sort(rand(poissrnd(2),1)*30-30,'descend');
%     sig_pow_dB = sort(rand(randi(4)-1,1)*15-15,'descend');
%     sig_pow_dB = [-3,-6,-9,-20,-20,-20];
    sig_pow0 = 10.^(sig_pow_dB/10);
    sig_pow = sig_pow0/sum(sig_pow0);
    pow_tot(band_realization_index) = sum(sig_pow);
    sig_table(1:length(sig_pow),3) = sig_pow';
    % we need to know how many blockers have been generated
    for bb=1:size(b_num,2)
        if sum(sig_table(:,3)>0)>=bb
            b_num(bb) = b_num(bb)+1;
        end
    end

    for ff=1:max(F_pool)
        best_bin_ref(ff) = (sig_table(ff,1)-1)*N/2+sig_table(ff,2);
    end
    
    pow_ref = ones(1,128)*1e-10;
    for bb=1:blk_num
        pow_ref((sig_table(bb,1)-1)*N/2+sig_table(bb,2))=sig_table(bb,3);
    end

    Nn = N;
    

    % fading_coeff = raylrnd(1,blk_num,runtimes);
    pd_counter = zeros(runtimes,length(F_pool),length(ave_num_pool));
%     pow_res = zeros(runtimes,length(F_pool),length(ave_num_pool));
    
    for runindex = 1:runtimes
        LL = N*max(ave_num_pool);
        y = zeros(LL,blk_num);
        CAL = randi(2,LL,band_num)*2-3;
        for ii=1:blk_num
            start_P = randi(L);
    %         sig_cache0 = sig(start_P:start_P+N-1,ii)*sqrt(10^(-randn));
            sig_cache0 = sig0(start_P:start_P+LL-1,ii);
            sig_cache1 = sig_cache0.*transpose(exp(1j*(rand*2*pi+(sig_table(ii,2)-1)*2*pi*(1:LL)/N)))*sqrt(sig_table(ii,3));
            if fading_flag
%                 chan = rayleighchan(1/6e6,100);
%                 sig_cache = filter(chan,sig_cache1);
                for zz=1:LL/N
                    sig_cache((zz-1)*N+1:zz*N,1) = randn*real(sig_cache1((zz-1)*N+1:zz*N))+1j*randn*imag(sig_cache1((zz-1)*N+1:zz*N));
                end
            else
                sig_cache = sig_cache1;
            end
            y(:,ii) = 2*real(sig_cache).*CAL(:,sig_table(ii,1))+abs(sig_cache).^2/sqrt(P_c);
        end
        sig_sum = sum(y,2)+sqrt(noise_power)*(rand(LL,1)*2-1);
        
        % fft-based periodogram estimation
        bins_bands_temp = zeros(max(ave_num_pool),N/2,band_num);
        for kk=1:band_num
            sig_sum_corr = sig_sum.*CAL(:,kk);
            for aa=1:max(ave_num_pool)
                periodo_gram = (abs(fft(sig_sum_corr((aa-1)*N+1:aa*N))/N).^2)';
                bins_bands_temp(aa,:,kk) = periodo_gram(1:N/2);
            end
        end
        
        for ave_index = 1:length(ave_num_pool)
            ave_num = ave_num_pool(ave_index);
            for kk=1:band_num
                bins_bands(:,kk) = sum(squeeze(bins_bands_temp(1:ave_num,:,kk)),1)'/ave_num;
            end

            % detection
            
            [sortval,sortindex] = sort((bins_bands(:)),'descend');
            for ff=1:length(F_pool)
                F=F_pool(ff);
%                 bestbin(runindex,:)=sortindex(1:F);
                if sort(sortindex(1:F))'==sort(best_bin_ref(1:F))
                    pd_counter(runindex,ff,ave_index) = pd_counter(runindex,ff,ave_index)+1;
                end
                pow_res(runindex,ave_index,ff,band_realization_index) = pow_tot(band_realization_index)-sum(pow_ref(sortindex(1:F)));
            end
        end
    end

%     for ave_index = 1:length(ave_num_pool)
%         for ff=1:length(F_pool)
%             pd(ave_index,ff,band_realization_index) = sum(pd_counter(:,ff,ave_index))/runtimes;
% %             pow_res1(ave_index,ff,band_realization_index) = sum(pow_res(:,ff,ave_index))/runtimes;
%         end
%     end
end
%% average over different realization of frequency settings
figure
plot(10*log10(1e-6+sort(reshape(squeeze(pow_res(:,1,3,:)),runtimes*band_realization_num,1))),linspace(0,1,runtimes*band_realization_num));hold on
plot(10*log10(1e-6+sort(reshape(squeeze(pow_res(:,end,3,:)),runtimes*band_realization_num,1))),linspace(0,1,runtimes*band_realization_num));hold on
grid on
%%
% figure
% plot(10*log10(1e-6+sort(pow_tot)),linspace(0,1,1000));hold on
% grid on