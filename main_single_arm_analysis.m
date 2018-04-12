%  08/15/2016
%  estimate interference floor and use it to set threshold for identifying
%  dominant blockers

clear;clc;%clf;close all


% rand('seed',3)
blk_num = 6;
M = 1e2;
N = 63;
SampleNumberAve = 32;
stat_num = 2e2;
P = M*N;
L = M*SampleNumberAve*stat_num;


for ii = 1:blk_num
%     sig_r(:,ii) = kron(randi(2,M*stat_num,1)*2-3,ones(SampleNumberAve,1));
%     sig_i(:,ii) = kron(randi(2,M*stat_num,1)*2-3,ones(SampleNumberAve,1));
%     sig(:,ii) = (sig_r(:,ii)+sig_i(:,ii))/sqrt(2);
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
%     unnormalized = real(temp_data(end-L+1:end));
%     sig_r(:,ii) = unnormalized./sqrt(sum(unnormalized.^2)/L);
%     unnormalized = imag(temp_data(end-L+1:end));
%     sig_i(:,ii) = unnormalized./sqrt(sum(unnormalized.^2)/L);
%     sig_i(:,ii) = zeros(size(sig_r(:,ii)));
%     sig(:,ii) = (sig_r(:,ii)+sig_i(:,ii))/sqrt(2);
    sig(:,ii) = temp_data(end-L+1:end)./sqrt(temp_data(end-L+1:end)'*temp_data(end-L+1:end)/L);
end

    ii = 1;
%     sig_r(:,ii) = kron(randi(2,M*stat_num,1)*2-3,ones(SampleNumberAve,1));
%     sig_i(:,ii) = kron(randi(2,M*stat_num,1)*2-3,ones(SampleNumberAve,1));
%     sig(:,ii) = (sig_r(:,ii)+sig_i(:,ii))/sqrt(2);
    upsam = 128;
    L2 = M*upsam*stat_num;
    symbols = fix(L2*2/upsam);
    
    clearvars data temp_data
    hmod = modem.qammod('M', 4, 'InputType', 'integer');
    hdesign  = fdesign.pulseshaping(upsam,'Square Root Raised Cosine');
    hpulse = design(hdesign);
    data = randi(4,symbols,1)-1;
    data = modulate(hmod, data);
    data = upsample(data,upsam);
    temp_data = conv(data,hpulse.Numerator);
%     unnormalized = real(temp_data(end-L+1:end));
%     sig_r(:,ii) = unnormalized./sqrt(sum(unnormalized.^2)/L);
%     unnormalized = imag(temp_data(end-L+1:end));
%     sig_i(:,ii) = unnormalized./sqrt(sum(unnormalized.^2)/L);
%     sig_i(:,ii) = zeros(size(sig_r(:,ii)));
%     sig(:,ii) = (sig_r(:,ii)+sig_i(:,ii))/sqrt(2);
    sig_other = temp_data(end-L+1:end)./sqrt(temp_data(end-L+1:end)'*temp_data(end-L+1:end)/L);


%%
other_pow = 1e-4;
sig_pow0 = ones(blk_num,1)*other_pow;
sig_pow1 = ones(blk_num,1)*other_pow;

pow_pool = -10-linspace(0,18,10);
% pow_pool = -20
pow_num = 10;
TH = zeros(pow_num,1);
pd = zeros(pow_num,1);

for pow_index = 1:pow_num
j1_dB = pow_pool(pow_index);
sig_pow1(4) = 1;
sig_pow0(4) = 1;
% sig_pow1(2:6) = ones(5,1)*0.1;
% sig_pow0(2:6) = ones(5,1)*0.1;
sig_pow1(1) = 10^(j1_dB/10);
sig_pow0(1) = 0;

%% calibration sequences
% cal0=PNgenerator_v5(N,N,1);
% CAL_noperm=LowerRate_v2(cal0,P);
% CAL=zeros(P,N);
% for mm=1:M
%     CAL((mm-1)*N+1:mm*N,:)=CAL_noperm((mm-1)*N+1:mm*N,randperm(N));
% end

%% tagging
runtimes = 5e2;
for runindex=1:runtimes
%     CAL=zeros(P,N);
%     for mm=1:M
%         CAL((mm-1)*N+1:mm*N,:)=CAL_noperm((mm-1)*N+1:mm*N,randperm(N));
%     end
    CAL = randi(2,P,blk_num)*2-3;

    %CAL=CAL(:,randperm(N));
    start_point = randi(L-P-1,blk_num,1);
    for ii=1:blk_num
        sig_cache(:,ii) = sig(start_point(ii,1):start_point(ii,1)+P-1,randi(blk_num));
    end
    sig_cache(:,4) = sig_other(start_point(ii,1):start_point(ii,1)+P-1,1);
    
    sig_mix0_part1 = zeros(P,1);
    sig_mix1_part1 = zeros(P,1);

    sig_mix0_part2 = zeros(P,1);
    sig_mix1_part2 = zeros(P,1);

    sig_mix0 = zeros(P,1);
    sig_mix1 = zeros(P,1);
    PhaseShift = exp(1j*rand(1,blk_num));
    
    
    for ii=1:blk_num
        sig_mix0_part1 = sig_mix0_part1+sig_cache(:,ii).*CAL(:,ii)*PhaseShift(ii)*sqrt(sig_pow0(ii));
        sig_mix1_part1 = sig_mix1_part1+sig_cache(:,ii).*CAL(:,ii)*PhaseShift(ii)*sqrt(sig_pow1(ii));
    end

    for ii=1:blk_num
        for jj=1:blk_num
            if ii ~= jj
            carrier = exp(1j*((ii-jj)/32*2*pi*(1:P)'))*PhaseShift(ii);
            sig_mix0_part2 = sig_mix0_part2 + sqrt(sig_pow0(ii))*CAL(:,jj).*...
                real(sig_cache(:,ii).*carrier);
            sig_mix1_part2 = sig_mix1_part2 + sqrt(sig_pow1(ii))*CAL(:,jj).*...
                real(sig_cache(:,ii).*carrier);
            end
        end
    end
    sig_mix0 = sig_mix0_part1+sig_mix0_part2;
    sig_mix1 = sig_mix1_part1+sig_mix1_part2;
    
    %% original algorithm
        Sam_Num = 5e3;
        bhi = fir1(200,1/32,'low');
        ED_results1(pow_index,runindex) = norm(filter(bhi,1,sig_mix1(1:Sam_Num).*CAL(1:Sam_Num,1)));
        ED_results0(pow_index,runindex) = norm(filter(bhi,1,sig_mix0(1:Sam_Num).*CAL(1:Sam_Num,1))); 
    
end
%
temp = sort(ED_results0(pow_index,:),'ascend');
TH(pow_index) = temp(runtimes*0.9);
pd(pow_index) = sum(ED_results1(pow_index,:)>TH(pow_index))/runtimes;

% temp = sort(ED_results2(1,:),'ascend');
% TH2(pow_index) = temp(runtimes*0.9);
% pd2(pow_index) = sum(ED_results2(2,:)>TH2(pow_index))/runtimes;
% if pow_index==1
%    1 
% end
end
%%

figure
plot(pow_pool,pd);hold on
% plot(pow_pool,pd2);hold on
grid on
% legend('algorithm 1','algorithm 2')
%%
pow_index = 1;
figure
plot(sort(ED_results0(pow_index,:)));hold on
plot(sort(ED_results1(pow_index,:)));
%% Plot Detection Statistics in both hypothesis
pow_index = 1;
[f,xi] = ksdensity(ED_results0(pow_index,:));
mu = mean(ED_results0(pow_index,:));
sigma = sqrt(var(ED_results0(pow_index,:)));
pdfxdata = linspace(mu-4*sigma,mu+4*sigma,100);
figure
plot(pdfxdata, normpdf(pdfxdata,mu,sigma),'linewidth',2);hold on
plot(xi,f,'--','linewidth',2)

[f,xi] = ksdensity(ED_results1(pow_index,:));
mu = mean(ED_results1(pow_index,:));
sigma = sqrt(var(ED_results1(pow_index,:)));
pdfxdata = linspace(mu-4*sigma,mu+4*sigma,100);
plot(pdfxdata, normpdf(pdfxdata,mu,sigma),'linewidth',2);hold on
grid on
xlabel('Energy Detection Statistic')
ylabel('PDF')
plot(xi,f,'--','linewidth',2)
legend('H0, Theo.','H0, Sim.', 'H1, Theo.','H1, Sim.')
%%
% ydata = abs(fft(sig_mix1(1:Sam_Num).*CAL(1:Sam_Num,1)));
ydata = abs(fft(sig_cache(:,1)));
xdata = linspace(-pi,pi,Sam_Num);
figure
plot(xdata,[ydata(2501:5000);ydata(1:2500)])
% figure
% subplot(211)
% plot(abs(temp(51:1050)));hold on
% plot(abs(sig_cache(1:1000,4)))
% subplot(212)
% plot(abs(temp(51:1050)-sig_cache(1:1000,4)*PhaseShift(4)));
%%

% figure
% subplot(211)
% plot(angle(temp(51:1050)));hold on
% plot(angle(sig_cache(1:1000,4)*PhaseShift(4)))
%%
% figure
% subplot(311)
% plot(abs(sig_mix0(1:5000)));hold on
% subplot(312)
% plot(abs(y_j4_hat));
% subplot(313)
% plot(abs(y_j4_hat-sig_mix0(1:5000)));
