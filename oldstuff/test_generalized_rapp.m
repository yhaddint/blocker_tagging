clear;clc;
sig_in = linspace(0,4,1e3);
Vsat = 1;
P_range = [1,2,10,100];
% sig_out = (1+(sig_in./Vsat).^(2*P)).^(1/P);
figure
plot(sig_in,sig_in.^2,'linewidth',2);hold on 
for pp = 1:length(P_range)
    sig_out = get_rapp_square( sig_in,Vsat,P_range(pp));
    plot(sig_in,sig_out,'linewidth',2);hold on 
end
ylim([0,1.1])
grid on
legend('True Square', 'P = 1', 'P = 2', 'P = 10', 'P = 100')