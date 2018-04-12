clear;clc;%clf;close all

f_start = 500e6;
f_end = 5000e6;
N = 22500;
ave_num=1;
tone_num = linspace(f_start,f_end,N);
% t_window = 1e-3;
fs = 5000e6;
sig_table(:,1) = [2000;6000;12000;2100;6100];
candidate = [1951:2050,5951:6050,11951:12050];

N2 = 22500/4;

for ii=1:N
    carrier(:,ii) = sin(2*pi*tone_num(ii)/fs*(1:N2));
%     carrier(:,ii) = sin(rand*2*pi+2*pi*tone_num(ii)/fs*(1:N));
end

for ii = 1
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

sig_pow_dB = [-3,-3,-3,-20,-20]';
sig_table(1:3,3) = 10.^(sig_pow_dB/10);

M_num = 6;
for M_index = 1:M_num
    M = 25*M_index;
    
    runtimes = 5e3;
    pd_counter = zeros(1,3);
    y = zeros(M,3);
    for runindex=1:runtimes
        if mod(runindex,50)==1
            PN_code = randi(2,M,N2)*2-3;
            measure_mat = PN_code*carrier(:,candidate);
        end
        
        for ii=1:5
            start_P = randi(L-300);
            env = sig_env(start_P:start_P+2*ave_num-1,1)';
            sig(:,ii) = interp1(1:2*ave_num,env,linspace(1,2*ave_num,N2*ave_num)).*...
                sin(rand*2*pi+2*pi*tone_num(sig_table(ii,1))/fs*(1:N2*ave_num))*sqrt(sig_table(ii,3));
%             sig(:,ii) = sin(rand*2*pi+2*pi*tone_num(sig_table(ii,1))/fs*(1:N2*ave_num))*sqrt(sig_table(ii,3));
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
        index = sort(best(:,runindex));
        for cc=1:3
            if (abs(index(1)-50)<=(cc-1))*(abs(index(2)-150)<=(cc-1))*(abs(index(3)-250)<=(cc-1))
                pd_counter(cc) = pd_counter(cc)+1;
            end
        end
    end
    for cc=1:3
        pd(cc,M_index) = pd_counter(cc)/runtimes;
    end
end
%%
figure
plot((1:M_num)*25*ave_num,pd);
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
% figure
% plot(PN_code(1,:))
% figure
% plot(sig_spread(1,:))