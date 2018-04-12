clear;clc
M=10:1:300;
kappa=(1+1.64*sqrt(2./M))./(1-1.64*sqrt(2./M));
figure
plot_setting
semilogx(M,10*log10(kappa))
grid on