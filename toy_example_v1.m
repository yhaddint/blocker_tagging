clc;clear
runtimes = 2e3;
sig_pow_pool = -20:1:-11;
H1_twoarm = zeros(runtimes,10);
H0_twoarm = zeros(runtimes,10);
H1_ED = zeros(runtimes,10);
H0_ED = zeros(runtimes,10);

for pow_index = 1:10
    sig_pow_dB = sig_pow_pool(pow_index);
    sig_pow = 10^(sig_pow_dB/10);
    
    for runindex=1:runtimes
        P = fix(5e3/32);
        signal = (randn(P,1)+1j*randn(P,1))*sqrt(sig_pow/2);
        noise1 = (randn(P,1)+1j*randn(P,1))*sqrt(1/2);
        noise2 = (randn(P,1)+1j*randn(P,1))*sqrt(1/2);
        H1_twoarm(runindex,pow_index) = real((signal+noise1)'*(signal+noise2));
        H0_twoarm(runindex,pow_index) = real((noise1)'*(noise2));

        H1_ED(runindex,pow_index) = norm(signal+signal+noise1+noise2);
        H0_ED(runindex,pow_index) = norm(noise1+noise2);
    end

temp = sort(H0_twoarm,'ascend');
TH_twoarm = temp(runtimes*0.9);
pd_twoarm(pow_index) = sum(H1_twoarm(:,pow_index)>TH_twoarm)/runtimes;

temp = sort(H0_ED,'ascend');
TH_ED = temp(runtimes*0.9);
pd_ED(pow_index) = sum(H1_ED(:,pow_index)>TH_ED)/runtimes;
end
%%
figure
plot(sig_pow_pool,pd_twoarm);hold on
plot(sig_pow_pool,pd_ED)
xlabel('noise power (dB)')
legend('two arm','ED')
% plot(sort(H0_twoarm));hold on
% plot(sort(H1_twoarm))
%%
clear;clc
N = 5000;
pow_s = 0.05;
pow_n = 2;
theo_var = 2*pow_s*pow_n+2*pow_s^2+pow_n^2;
runtimes = 2e3;
for runindex = 1:runtimes
s = randn(N,1)*sqrt(pow_s);
n1 = randn(N,1)*sqrt(pow_n);
n2 = randn(N,1)*sqrt(pow_n);
ED_H0(runindex) = n1.'*n2/N;
ED_H1(runindex) = (s+n1).'*(s+n2)/N;
end
temp = sort(ED_H0,'ascend');
TH_ED = temp(runtimes*0.9)
TH_theo = qfuncinv(0.1)*pow_n/sqrt(N)
pd_theo = qfunc((-pow_s+TH_theo)/(sqrt((theo_var)/N)))
pd = sum(ED_H1>TH_ED)/runtimes
% mean((s+n1).*(s+n2))
% var((s+n1).*(s+n2))
%%
clear;clc
N = 6e3;
bhi = fir1(100,1/32,'low');
runtimes = 5e2;
for runindex = 1:runtimes
n0 = randn(N,1);

n0_ave = filter(bhi,1,n0);

mu(runindex) = norm(n0_ave(101:5100))^2/5e3;
end
%%
var(mu)
2*(mean(mu))^2/(5e3/32)



