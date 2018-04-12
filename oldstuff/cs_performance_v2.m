clear;clc;%clf;close all

f_start = 500e6;
f_end = 5000e6;
N = 22500;
ave_num=1;
tone_num = linspace(f_start,f_end,N);
% t_window = 1e-3;
fs = 5000e6;

rand('seed',4);
N2 = 22500/4;

for ii = 1:7
%     sig_r(:,ii) = kron(randi(2,M*stat_num,1)*2-3,ones(N,1));
%     sig_i(:,ii) = kron(randi(2,M*stat_num,1)*2-3,ones(N,1));
%     sig(:,ii) = (sig_r(:,ii)+1j*sig_i(:,ii))/sqrt(2);
    upsam = 100;
    symbols = 3e2;
    
    clearvars data temp_data
    hmod = modem.qammod('M', 4, 'InputType', 'integer');
    hdesign  = fdesign.pulseshaping(upsam,'Square Root Raised Cosine');
    hpulse = design(hdesign);
    data = randi(4,symbols,1)-1;
    data = modulate(hmod, data);
    data = upsample(data,upsam);
    temp_data = conv(data,hpulse.Numerator);
    sig0(:,ii) = temp_data./sqrt(sum(temp_data.*conj(temp_data))/length(temp_data));
    sig_env(:,ii) = real(sig0(2000:end-2000,ii));
end
L = length(sig_env);
runtimes = 50;
band_realization_num = 500;
M_num = 12;

pow_res = zeros(3,M_num,runtimes,band_realization_num);
pd = zeros(3,3,M_num,band_realization_num);

for ii=1:N
    carrier(:,ii) = sin(2*pi*tone_num(ii)/fs*(1:N2));
%     carrier(:,ii) = sin(rand*2*pi+2*pi*tone_num(ii)/fs*(1:N));
end

subband_num = 20;

candidate = [];
for ii=1:8
    candidate = [candidate,ii*1000+1:ii*1000+subband_num];
end

candidate_noexcess = [];
for ii=1:8
    candidate_noexcess = [candidate_noexcess,(ii-1)*subband_num+3:ii*subband_num-2];
end

for b_index=1:band_realization_num
    disp(b_index/band_realization_num*100)
    blk_num = min(poissrnd(2),7);
    blk_num_rec(b_index) = blk_num;
    sig_table = zeros(blk_num,2);
    sig_pow_dB = sort(rand(blk_num,1)*30-30,'descend');
%     sig_pow_dB = sort(rand(4,1)*30-30,'descend');
    sig_pow0 = 10.^(sig_pow_dB/10);
    sig_pow = sig_pow0/sum(sig_pow0);
    sig_table(1:blk_num,3) = sig_pow;
    band_ref = randi(8,blk_num,1);
    subband_ref = randi(16,blk_num,1)+2;
    sig_table(1:blk_num,1) = band_ref*1000+subband_ref;

    pow_ref = zeros(200,1);
    sig_table(:,2) = subband_ref+[(band_ref-1)*subband_num];
    pow_ref(sig_table(:,2)) = sig_pow;
    pow_tot(b_index) = sum(sig_pow);
    
    

    for M_index = 1:M_num
        M = 25*M_index;
        pd_counter = zeros(3,3);
        y = zeros(M,3);
        for runindex=1:runtimes
            if mod(runindex,50)==1
                PN_code = randi(2,M,N2)*2-3;
                measure_mat = PN_code*carrier(:,candidate);
            end

            for ii=1:blk_num
                start_P = randi(L-300);
%                 env = sig_env(start_P:start_P+2*ave_num-1,ii)';
%                 sig(:,ii) = interp1(1:2*ave_num,env,linspace(1,2*ave_num,N2*ave_num)).*...
%                     sin(rand*2*pi+2*pi*tone_num(sig_table(ii,1))/fs*(1:N2*ave_num))*sqrt(sig_table(ii,3));
                sig(:,ii) = sig_env(start_P,ii)*sin(rand*2*pi+2*pi*tone_num(sig_table(ii,1))/fs*(1:N2*ave_num))*sqrt(sig_table(ii,3));
%                 sig(:,ii) = sin(rand*2*pi+2*pi*tone_num(sig_table(ii,1))/fs*(1:N2*ave_num))*sqrt(sig_table(ii,3));
            end

            y(:,1) = PN_code*(sum(reshape(sum(sig,2),N2,ave_num),2)/ave_num);
    %         for mm=1:M
    %             sig_spread(mm,:) = PN_code(mm,:).*sig;
    %             y(mm,1) = sum(sig_spread(mm,:));

    %         end
            %% reconstruction
            for ii=1:3
                phi(:,ii) = measure_mat'*y(:,ii);
                [bestscore(ii),best(ii,runindex)] = max(abs(phi(:,ii)));
    % %             regular MP updata
    %             coeff(ii) = pinv(measure_mat(:,best(ii)))*y(:,ii);
    %             y(:,ii+1) = y(:,ii) - coeff(ii)*measure_mat(:,best(ii));
                y(:,ii+1) = (eye(M)-measure_mat(:,best(ii,runindex))*measure_mat(:,best(ii,runindex))'/norm(measure_mat(:,best(ii,runindex)))^2)*y(:,ii);
            end
            index = best(:,runindex);
            
            for cc=1:3
                switch blk_num
                    case 0
                    case 1
                        refindex = sig_table(1,2);
                        flag(1) = length(intersect(index(1),[refindex(1)-(cc-1):refindex(1)+(cc-1)]))>0;
                        pd_counter(1,cc) = pd_counter(1,cc)+flag(1); 
                    case 2
                        refindex = sig_table(1:2,2);
                        flag(1) = length(intersect(index(1),[refindex(1)-(cc-1):refindex(1)+(cc-1)]))>0;
                        flag(2) = length(intersect(index(1:2),[refindex(2)-(cc-1):refindex(2)+(cc-1)]))>0;
                        pd_counter(1,cc) = pd_counter(1,cc)+flag(1);
                        pd_counter(2,cc) = pd_counter(2,cc)+flag(1)*flag(2); 
                    otherwise
                        refindex = sig_table(1:3,2);
                        flag(1) = length(intersect(index(1),[refindex(1)-(cc-1):refindex(1)+(cc-1)]))>0;
                        flag(2) = length(intersect(index(1:2),[refindex(2)-(cc-1):refindex(2)+(cc-1)]))>0;
                        flag(3) = length(intersect(index(1:3),[refindex(3)-(cc-1):refindex(3)+(cc-1)]))>0;
                        pd_counter(1,cc) = pd_counter(1,cc)+flag(1);
                        pd_counter(2,cc) = pd_counter(2,cc)+flag(1)*flag(2);
                        pd_counter(3,cc) = pd_counter(3,cc)+flag(1)*flag(2)*flag(3);  
                        
                end
                pow_res(cc,M_index,runindex,b_index) = sum(sig_pow);
                pow_ref0 = pow_ref;
                for ii=1:3
                    if sum(index(ii)==candidate_noexcess)
                        pow_res(cc,M_index,runindex,b_index) = pow_res(cc,M_index,runindex,b_index)-sum(pow_ref0(index(ii)-(cc-1):index(ii)+(cc-1)));
                        pow_ref0(index(ii)-(cc-1):index(ii)+(cc-1)) = 0;
                    end
                end
            end
        end
        for cc=1:3
            for ff=1:3
                pd(ff,cc,M_index,b_index) = pd_counter(ff,cc)/runtimes;
            end
        end
    end
end
%%
for ff=1:3
    for cc=1:3
        for M_index=1:M_num
            pd_summary(ff,cc,M_index) = mean(squeeze(pd(ff,cc,M_index,:)));
        end
    end
end
figure
plot((1:M_num)*25*ave_num,squeeze(pd_summary(1,1,:))*band_realization_num/sum(blk_num_rec>=1));hold on
plot((1:M_num)*25*ave_num,squeeze(pd_summary(2,1,:))*band_realization_num/sum(blk_num_rec>=2));hold on
plot((1:M_num)*25*ave_num,squeeze(pd_summary(3,1,:))*band_realization_num/sum(blk_num_rec>=3));hold on
% plot((1:M_num)*25*ave_num,squeeze(pd_summary(1,3,:))*band_realization_num/sum(blk_num_rec>=1));hold on
% plot((1:M_num)*25*ave_num,squeeze(pd_summary(2,3,:))*band_realization_num/sum(blk_num_rec>=2));hold on
% plot((1:M_num)*25*ave_num,squeeze(pd_summary(3,3,:))*band_realization_num/sum(blk_num_rec>=3));hold on
grid
% xlabel('Number of Samples')
xlabel('Number of Samples')
ylabel('Ident. Prob.')
ylim([0,1])
% %%
% figure
% plot(abs(temp(:,1)));hold on
% figure
% plot(abs(temp(:,2)))
%%
for mm=1:M_num
    for cc=1:3
        pow_res1(cc,mm) = (mean(mean(squeeze(pow_res(cc,mm,:,:)))));
    end
end
figure
plot((1:M_num)*25*ave_num,10*log10(pow_res1));hold on
plot((1:M_num)*25*ave_num,ones(1,M_num)*10*log10(mean(pow_tot)));
grid
%%
figure
plot(32*(1:10),[-6.704,-9.283,-10.67,-11.52,-12.04,-12.48,-12.76,-13.01,-13.2,-13.4])