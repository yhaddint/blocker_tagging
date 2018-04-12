clear;clc;
B0_pool=[1.5,3,6];
eta_pool=[23,20,17];
%for bb=1:length(B0_pool)
for bb=1:length(eta_pool)
    eta=10^(eta_pool(bb)/10);
    %B0=B0_pool(bb);
    B0=2.5;
    mu0=10.^(linspace(-3,-2.3,100));
    M(bb,:) = (1.68*sqrt(B0-1)./(2./(1+mu0.*(eta-1))-1)).^2;
end
%%
figure
plot_setting
loglog(mu0,M)
grid on
legend('d3 = 23dB','d3 = 20dB','d3 = 17dB')
xlabel('mu0')
ylabel('Required Number of Frams (M)')
ylim([0,500])
