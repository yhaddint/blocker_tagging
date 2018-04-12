%  12/16/2015
%  test how many samples are required for STPE method in finding strongest
%  blocker.
%  We add rayleigh fading to interferers
%  % instead of looking for average residual power, Now we watch CDF
% the residual power




clear;clc;%clf;close all

band_num = 8;
blk_num = 10;
M = 1e2;
N = 32;
SampleNumberAve = 32;

shadowing=0;

stat_num = 3e2;
P = M*N;
L = N*stat_num;

noise_power = 1e-3;

rand('seed',4);
for ii = 1:blk_num
%     sig_r(:,ii) = kron(randi(2,M*stat_num,1)*2-3,ones(N,1));
%     sig_i(:,ii) = kron(randi(2,M*stat_num,1)*2-3,ones(N,1));
%     sig0(:,ii) = (sig_r(:,ii)+1j*sig_i(:,ii))/sqrt(2);
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
    sig(:,ii) = temp_data./sqrt(sum(temp_data.*conj(temp_data))/length(temp_data));

end

%% Formulating multi-band signals
band_realization_num = 1000;
runtimes = 20;
F_pool = [1 2 3];
F_num = length(F_pool);
ave_pool = 1:8;
max_ave = max(ave_pool);
ave_num = length(ave_pool);
pd = zeros(ave_num,F_num,band_realization_num);
b_num=zeros(1,10);
pow_res = zeros(runtimes,ave_num,F_num,band_realization_num);

for band_realization_index=1:band_realization_num
    disp(band_realization_index/band_realization_num*100)
    % Formulating multi-band signals
    % signal table (band | subband | power)
    sig_table = zeros(blk_num,3);
    % sig in which band
    sig_table(:,1) = randi(band_num,blk_num,1);
    % sig in which subband
    sig_table(:,2) = randi(N/2-1,blk_num,1);
    % power of sig
    
    P_lr_dB = -10;
    P_lr = 10^(P_lr_dB/10);
    sig_pow = zeros(blk_num-1,1);
    while sum(sig_pow>=P_lr)<3
        sig_pow_dB = sort(rand(3,1)*(-P_lr_dB-2)+P_lr_dB,'descend');
        sig_pow0 = 10.^(sig_pow_dB/10);
        sig_pow = sig_pow0./sum(sig_pow0);
    end
    
%     sig_pow_dB = sort(rand(poissrnd(2),1)*30-30,'descend');
%     sig_pow_dB = sort(rand(3,1)*10-10,'descend');
%     sig_pow0 = 10.^(sig_pow_dB/10);
%     sig_pow = sig_pow0/sum(sig_pow0);
    pow_tot(band_realization_index) = sum(sig_pow);
    sig_table(1:length(sig_pow),3) = sig_pow';
    pow_ref = ones(1,128)*1e-10;
    for bb=1:blk_num
        pow_ref((sig_table(bb,1)-1)*N/2+sig_table(bb,2))=sig_table(bb,3);
    end
    
    % we need to know how many blockers have been generated
    for bb=1:size(b_num,2)
        if sum(sig_table(:,3)>0)>=bb
            b_num(bb) = b_num(bb)+1;
        end
    end
    
    
    %sig_conv = zeros(runtimes,L);
    bins_bands = zeros(runtimes,N/2,band_num,ave_num);
    for runindex = 1:runtimes
        %start_point = randi(P*M,1,blk_num);
        for kk=1:band_num
            sig_sum0 = zeros(N*max_ave,1);
            for ii=1:blk_num
                if sig_table(ii,1)==kk
                    for aa=1:max_ave
                        start_P = randi([10*SampleNumberAve,size(sig,1)-10*SampleNumberAve-1]);
                        y_temp0 = sig(start_P:start_P+N-1,ii);
                        if shadowing
                            y_temp1 = randn*real(y_temp0)+1j*randn*imag(y_temp0);
                        else
                            y_temp1 = y_temp0;
                        end
                        y1((aa-1)*N+1:aa*N) = y_temp1.*transpose(exp(1j*(rand*2*pi+(sig_table(ii,2)-1)*2*pi*(1:N)/N)))*...
                            sqrt(sig_table(ii,3));
                    end
%                     y1 = y;
                    sig_sum0 = sig_sum0 + transpose(y1);
                end
            end
            sig_sum = sig_sum0+sqrt(noise_power)*(rand(max_ave*N,1)*2-1);
%             sig_sum = sig_sum0;
            for ave_index = 1:max_ave
                periodo_gram(:,ave_index) = (abs(fft(sig_sum((ave_index-1)*N+1:ave_index*N)/N)).^2)';
            end
            for ave_index = 1:ave_num
                ave = ave_pool(ave_index);
                bins_bands(runindex,:,kk,ave_index) = sum(periodo_gram(1:N/2,1:ave),2)/ave;
            end
        end
    end

    % ident. prob. calculation
    best_bin_ref = zeros(1,max(F_pool));
    for ff=1:max(F_pool)
        best_bin_ref(ff) = (sig_table(ff,1)-1)*N/2+sig_table(ff,2);
    end
    for F_index = 1:F_num
        F = F_pool(F_index);

        bestbin = zeros(runtimes,F);
        A = zeros(N/2,band_num);
        for ave_index = 1:ave_num
            clearvars bestscore bestbin
            pd_counter = 0;
            for ii=1:runtimes
               A = abs(squeeze(bins_bands(ii,:,:,ave_index)));
               [sortval,sortindex] = sort((A(:)),'descend');
               bestbin(ii,:)=sortindex(1:F);
               if sort(sortindex(1:F))'==sort(best_bin_ref(1:F))
                    pd_counter = pd_counter+1;
               end
               pow_res(ii,ave_index,F_index,band_realization_index) = pow_tot(band_realization_index)-sum(pow_ref(sortindex(1:F)));
            end
            pd(ave_index,F_index,band_realization_index) = pd_counter/runtimes;
        end
    end
    


end

%% average over different realization of frequency settings
temp0 = 10*log10(1e-6+sort(reshape(squeeze(pow_res(:,1,3,:)),band_realization_num*runtimes,1)));
temp1 = [ones(5,1)*temp0(1);temp0;temp0(end);ones(5,1)*temp0(end)];
for ii=1:length(temp0)
    temp2(ii) = mean(temp1(ii:ii+10));
end
figure
plot(temp2,linspace(0,1,band_realization_num*runtimes));
hold on

temp0 = 10*log10(1e-6+sort(reshape(squeeze(pow_res(:,end,3,:)),band_realization_num*runtimes,1)));
temp1 = [ones(5,1)*temp0(1);temp0;temp0(end);ones(5,1)*temp0(end)];
for ii=1:length(temp0)
    temp2(ii) = mean(temp1(ii:ii+10));
end

plot(temp2,linspace(0,1,band_realization_num*runtimes));
hold on
grid on
%%
figure
plot([ave_pool]*N*band_num,mean(pd,3))
