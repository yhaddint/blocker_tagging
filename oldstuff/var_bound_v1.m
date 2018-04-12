%  05/14/2015
%  this script test variance of estimation by Monte Carlo Simulation
%  toghther with my theoretically upper bound
clear all;clc;

runtimes_sim = 1e4;
rrcPulseShape=1;

blk_num = 31;
M = 1e2;
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
        sig(:,ii) = (sig_bb_r(:,ii)+sig_bb_i(:,ii))*sqrt(sig_pow(ii));
    end
end
%% MC estimation variance
for runindex=1:runtimes_sim
    if mod(runindex,fix(runtimes_sim/20))==0
        runindex/fix(runtimes_sim/20)*5
    end
    CAL = randi(2,M*N,blk_num)*2-3;

    start_point = randi(L-P-1,blk_num);
    for ii=1:blk_num
        sig_cache(:,ii) = sig(start_point(ii):start_point(ii)+N*M-1,ii);
    end
    sig_mix = zeros(P,1);
    for ii=1:blk_num
        sig_mix = sig_mix + sig_cache(:,ii).*CAL(:,ii);
    end
    
    for ii=1:blk_num
        result_sim(runindex,ii) = tagging_v3(sig_mix,CAL(:,ii),N,M);
    end
end
%
for ii=1:blk_num
    mean_sim(ii) = mean(result_sim(:,ii));
    var_sim(ii) = var(result_sim(:,ii));
end
%% analytical upper bound for estiamtion variance
% mus = 0.6503;
% mu0 = 0.0323;
% sigmas = 0.1591;
% sigma0 = 0.002025;
mus = 1.54;
mu0 = 2/31;
sigmas = 1.215;
sigma0 = 0.01031;
for ii=1:blk_num
    vech=ones(blk_num,1)*mu0;
    vech(ii)=mus;
    vSigma=ones(blk_num,1)*(sigma0-2*mu0^2);
    vSigma(ii)=(sigmas-2*mus^2);
    pow_square_sum=sum(sig_pow.^2);
    var_anal1(ii)=(2*(sig_pow'*vech)^2+sig_pow(ii)^2*(-2*mus^2))/M;%+sig_pow'*diag(vSigma)*sig_pow)/M;
    var_anal2(ii)=(2*(sig_pow'*vech)^2+sig_pow(ii)^2*2)/M;%+sig_pow'*diag(vSigma)*sig_pow)/M;
end
var_bound = 2*mean_sim.^2/M;
% plot
figure
semilogy(sqrt(var_sim),'b');hold on
semilogy(sqrt(var_anal1),'r');hold on
semilogy(sqrt(var_anal2),'k');hold on
semilogy(sqrt(var_bound),'g');hold on
xlabel('BLK index');
ylabel('Var');
legend('sim.','anal.','bound')
%%
figure
semilogy(mean(result_sim),'-x')