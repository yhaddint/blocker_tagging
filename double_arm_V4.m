%  01/17/2017
%  Scatter plot of tagging performance

clear;clc;%clf;close all
rng('default')

blk_num = 2;

% oversampling ratio, which is PN rate v.s. blocker rate
OverSampling = 8;

total_time = 1e-3; % maximum running time is 1ms
Lmax = fix(total_time/(1/6e6/OverSampling));

% Length of PN per period
M = 1600;
N = fix(Lmax/M);
P = M*N;
L = 30*P;

Nrange = fix(P*(0.1:0.1:0.9));

Pc = 10; % power of calibration signal (mW)

% generating complex baseband signal, which is the sig_bb used later
for ii = 1:blk_num
    upsam = OverSampling;
    bhi = fir1(400,1/upsam*0.7,'low');
    symbols = fix(L*2/upsam);
    clearvars data temp_data
    hmod = modem.qammod('M', 4, 'InputType', 'integer');
    hdesign  = fdesign.pulseshaping(upsam,'Square Root Raised Cosine');
    hpulse = design(hdesign);
    data = randi(4,symbols,1)-1;
    data = modulate(hmod, data);
    data = upsample(data,upsam);
    temp_data = conv(data,hpulse.Numerator);
    temp1_data = filter(bhi,1,temp_data);
    sig(:,ii) = temp1_data(end-L+1:end)./sqrt(temp1_data(end-L+1:end)'*temp1_data(end-L+1:end)/L);
end


%%
ReusePN = 1;
ofNum = blk_num+1;
MCtimes = 1e3;
results = zeros(ofNum,blk_num,MCtimes,length(Nrange));
% results_H1 = zeros(ofNum,blk_num,runtimes,length(Nrange));
% blocker power settings


sig_pow_all = zeros(blk_num,MCtimes);

for MCindex=1:MCtimes
    sig_on = rand(blk_num,1)<0.7;
    while sum(sig_on)==0
    sig_on = rand(blk_num,1)<0.7;
    end
    sig_pow_dB = rand(blk_num,1)*30-30;
    sig_pow = sig_on.*10.^(sig_pow_dB/10)+(1-sig_on).*10.^(ones(blk_num,1)*(-90/10));
    sig_pow_all(:,MCindex) = sig_pow;
    % generating PN sequence
    if ReusePN
        CALseed = randi(2,M+1,blk_num)*2-3;
        CAL_u0 = kron(ones(N,1),CALseed);
        CAL_u = CAL_u0(1:P,:);
        CALseed = randi(2,M+2,blk_num)*2-3;
        CAL_l0 = kron(ones(N,1),CALseed);
        CAL_l = CAL_l0(1:P,:);
    else
        CAL_u = randi(2,P,blk_num)*2-3;
        CAL_l = randi(2,P,blk_num)*2-3;
    end
    
    % pilots contains two co-channel PN sequences passing through diode
    d = randi(5e3);
    pilot_u = 0.5*Pc*[zeros(d,1);CALseed(:,1)].*conj([zeros(d,1);CALseed(:,1)])...
        +0.5*Pc*[zeros(d,1);CALseed(:,2)].*conj([zeros(d,1);CALseed(:,2)])...
        +Pc*real([zeros(d,1);CALseed(:,1)].*conj([zeros(d,1);CALseed(:,2)]));
    [peak,dhat] = max(abs(conv(pilot_u,flipud(CALseed(:,1).*CALseed(:,2)))));
    % estimated delay offset that will be used in tagging
    dOS = dhat-size(CALseed(:,1),1);

    % Simulating Diode after blockers are allowed to bt input
    sig_mix_u = zeros(2*P,1);
    sig_mix_l = zeros(2*P,1);
    start_point = randi(L-P-1,blk_num,1);
    for ii=1:blk_num
        sig_bb(:,ii) = sig(start_point(ii,1):start_point(ii,1)+P-1,randi(blk_num));
    end
    Phi = exp(1j*rand(1,blk_num));
    for kk=1:blk_num
        for ll=1:blk_num
            if kk==ll
                sig_mix_u(d+1:d+P) = sig_mix_u(d+1:d+P)...
                    +sqrt(Pc*sig_pow(kk))*real(sig_bb(1:P,kk)*Phi(kk)).*CAL_u(1:P,kk)...
                    +0.5*sig_pow(kk)*(sig_bb(1:P,kk).*conj(sig_bb(1:P,kk)));
                sig_mix_l(d+1:d+P) = sig_mix_l(d+1:d+P)...
                    +sqrt(Pc*sig_pow(kk))*real(sig_bb(1:P,kk)*Phi(kk)).*CAL_l(1:P,kk)...
                    +0.5*sig_pow(kk)*(sig_bb(1:P,kk).*conj(sig_bb(1:P,kk)));
            elseif abs(kk-ll)<OverSampling
                CFO = exp(1j*((kk-ll)/OverSampling*2*pi*(1:P)'));
                % Assuming band spacing DeltaF, such that fs/DeltaF =
                % OverSampling
                sig_mix_u(d+1:d+P) = sig_mix_u(d+1:d+P)...
                    + sqrt(Pc*sig_pow(kk))*CAL_u(1:P,ll).*real(sig_bb(1:P,kk).*CFO * Phi(kk))...
                    + sqrt(sig_pow(kk)*sig_pow(ll))*real(sig_bb(1:P,kk).*conj(sig_bb(1:P,ll)).*CFO * Phi(kk) * conj(Phi(ll)));
                sig_mix_l(d+1:d+P) = sig_mix_l(d+1:d+P)...
                    + sqrt(Pc*sig_pow(kk))*CAL_l(1:P,ll).*real(sig_bb(1:P,kk).*CFO * Phi(kk))...
                    + sqrt(sig_pow(kk)*sig_pow(ll))*real(sig_bb(1:P,kk).*conj(sig_bb(1:P,ll)).*CFO * Phi(kk) * conj(Phi(ll)));
            end
        end
    end      

    % Digital multiplication part, delay offset is adjusted
    bhi = fir1(200,1.1/OverSampling,'low');
    for bb=1:blk_num
        yout_u = filter(bhi,1,sqrt(2/Pc)*sig_mix_u(dOS+1:dOS+P).*CAL_u(1:P,bb));
        yout_l = filter(bhi,1,sqrt(2/Pc)*sig_mix_l(dOS+1:dOS+P).*CAL_l(1:P,bb));
        Sam_Num = Nrange(2);
        results(bb,MCindex) = yout_u(201:(200+Sam_Num))'*...
            yout_l(201:(200+Sam_Num))/Sam_Num;

    end

    % Double-arm multiplication for tagging
    for bb=1:blk_num
        Sam_Num = Nrange(2);
        sumpow_hat = sum(squeeze(results(:,MCindex)));
        sumpow_hat_res(bb,MCindex) = sumpow_hat;
        sigma_hat = sumpow_hat*(blk_num)/OverSampling...
                        /sqrt(Sam_Num/OverSampling);
        TH = qfuncinv(0.05)*sigma_hat;
        ptag(bb,MCindex) = results(bb,MCindex)>TH;
    end
end


%% Plotting tagging rate curve
figure
for MCindex=1:MCtimes
    for bb=1:blk_num
        xdata = 10*log10(sum(sig_pow_all(:,MCindex)));
        ydata = 10*log10(sig_pow_all(bb,MCindex));
        if ptag(bb,MCindex)
            plot(xdata,ydata,'ro');hold on
        else
            plot(xdata,ydata,'bx');hold on
        end
    end
end
xlabel('Total Power (dBm)')
xlim([-20,5])
ylabel('Interferer Power (dBm)');
grid on
%% false alarm counting
H0_count = 0;
fa_count = 0;
for MCindex=1:MCtimes
    for bb=1:blk_num
        if 10*log10(sig_pow_all(bb,MCindex))==-90
                H0_count = H0_count+1;
            if ptag(bb,MCindex)
                fa_count = fa_count+1;
            end
        end
    end
end
fa_count/H0_count