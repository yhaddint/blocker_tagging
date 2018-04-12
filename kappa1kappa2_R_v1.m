clear;clc
f = openfig('figure/data_achive.fig');
h = findobj(f,'Type','line');
%xdata = get(h,'xdata');
ydata = get(h,'ydata');
close all
B0 = 3;
M = 30;
Qratio = (1+qfuncinv(0.05)*sqrt((B0-1)/M))/(1-qfuncinv(0.05)*sqrt((B0-1)/M));
line_type = ['--';'-.';'  '];
figure(1)
plot_setting

for Nindex=1:3
upsilon_pool = ydata{Nindex+3};
mu_pool = ydata{Nindex};
r_step_num = 41;
eta_dB = linspace(0,40,r_step_num);
eta = 10.^(eta_dB/10);
%%
r_num = length(mu_pool);

kappa_product = zeros(r_num,r_step_num);

for rindex=1:r_num
    upsilon0 = upsilon_pool(rindex);
    mu0 = mu_pool(rindex);

    kappa_product(rindex,:) = Qratio*(1+mu0.*(eta-1))./upsilon0;
end
%%
figure(1)
plot(eta_dB,10*log10(kappa_product(end-3:end,:)),line_type(Nindex,:));hold on
xlabel('Power Range \eta (dB)')
ylabel('Transition Margin \kappa (dB)')
end
figure(1)
grid on
legend('N=63, B_0=125KHz','N=63, B_0=250KHz','N=63, B_0=500KHz','N=63, B_0=1MHz',...
    'N=31, B_0=125KHz','N=31, B_0=250KHz','N=31, B_0=500KHz','N=31, B_0=1MHz',...
    'N=15, B_0=125KHz','N=15, B_0=250KHz','N=15, B_0=500KHz','N=15, B_0=1MHz')
