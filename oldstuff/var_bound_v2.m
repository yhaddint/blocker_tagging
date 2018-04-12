%  05/22/2015
%  this script test variance of estimation by Monte Carlo Simulation
%  toghther with my theoretically upper bound
clear all;clc;

runtimes_sim = 1e4;
rrcPulseShape=1;

blk_num = 31;
M = 20;
N = 31;
SampleNumberAve = 30;
stat_num = 2e3;
P = M*N;
L = M*SampleNumberAve*stat_num;

sig_pow_dB = linspace(0,60,blk_num)'-90;
sig_pow_dB = sort(sig_pow_dB,'descend');
sig_pow = 10.^(sig_pow_dB/10);

if rrcPulseShape==0
    for ii=1:blk_num
        sig_r(:,ii) = kron(randi(2,M*stat_num,1)*2-3,ones(SampleNumberAve,1));
        sig_i(:,ii) = kron(randi(2,M*stat_num,1)*2-3,ones(SampleNumberAve,1));
        sig_i(:,ii) = zeros(size(sig_r(:,ii)));
        sig(:,ii) = (sig_r(:,ii)+sig_i(:,ii))*sqrt(sig_pow(ii));
    end
else
    upsam = SampleNumberAve;
    symbols = fix(L*2/upsam);
    hmod = modem.pskmod('M', 4, 'InputType', 'integer');
    hdesign  = fdesign.pulseshaping(upsam,'Square Root Raised Cosine');
    hpulse = design(hdesign);
    for ii=1:blk_num
        ii
        clearvars data temp_data
        data = randi(4,symbols,1)-1;
        data = modulate(hmod, data);
        data = upsample(data,upsam);
        temp_data = conv(data,hpulse.Numerator);
        unnormalized = real(temp_data(end-L+1:end));
        sig_bb_r(:,ii) = unnormalized./sqrt(sum(unnormalized.^2)/L);
        unnormalized = imag(temp_data(end-L+1:end));
        sig_bb_i(:,ii) = unnormalized./sqrt(sum(unnormalized.^2)/L);
        sig(:,ii) = (sig_bb_r(:,ii)+sig_bb_i(:,ii))/sqrt(2)*sqrt(sig_pow(ii));
    end
end
%% sweeping one particular blk's bandwidth
R_range = 10:20:210;

for Rindex=1:length(R_range)

    SampleperSymbol=R_range(Rindex);
    SampleperFrame=M*SampleperSymbol;

    hdesign  = fdesign.pulseshaping(SampleperSymbol,'Square Root Raised Cosine');
    hpulse = design(hdesign);
    %  Define length of BLK we want.
    stat_num=40;
    L=SampleperFrame*stat_num;
    clearvars sig_bb BLK_r BLK_i BLK_sweep
    for ii=1:5
        symbols=1500;
        hmod = modem.pskmod('M', 4, 'InputType', 'integer');
        data = randi(4,symbols,1)-1;
        data = modulate(hmod, data);
        data = upsample(data,SampleperSymbol);
        temp_data = conv(data,hpulse.Numerator);
        sig_bb(:,ii) = temp_data(end-(symbols-5)*SampleperSymbol+1:end-5*SampleperSymbol+1);

        BLK_r(:,ii)=real(sig_bb(:,ii))./sqrt(mean(real(sig_bb(:,ii)).^2));
        BLK_i(:,ii)=imag(sig_bb(:,ii))./sqrt(mean(imag(sig_bb(:,ii)).^2));
        BLK_sweep(:,ii)=(BLK_r(:,ii)+BLK_i(:,ii))/sqrt(2)*sig_pow(1);
    end

% MC estimation variance
    for runindex=1:runtimes_sim
        if mod(runindex,fix(runtimes_sim/20))==0
            runindex/fix(runtimes_sim/20)*5
        end
        CAL = randi(2,M*N,blk_num)*2-3;

        start_point = randi(L-P-1,blk_num);
        for ii=1:blk_num
            sig_cache(:,ii) = sig(start_point(ii):start_point(ii)+N*M-1,ii);
        end
        blk_cache = BLK_sweep(start_point(ii):start_point(ii)+N*M-1,randi(5));
        sig_mix = blk_cache.*CAL(:,1);
        for ii=2:blk_num
            sig_mix = sig_mix + sig_cache(:,ii).*CAL(:,ii);
        end

        for ii=1:blk_num
            result_sim(runindex) = tagging_v3(sig_mix,CAL(:,1),N,M);
        end
    end

    var_sim(Rindex) = var(result_sim);

end
%% plot
figure
plot(R_range,sqrt(var_sim));hold on
grid on
xlabel('R');
ylabel('Var');

