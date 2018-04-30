clear;clc;
PNvsBLK_ratio_range = 2.^linspace(2,6,2e2);
BLK_powerratio_range_dB = -linspace(10,30,2e2);
BLK_powerratio_range = 10.^(BLK_powerratio_range_dB/10);
pfa_goal = 0.05;
pm_goal = 0.05;

for rr = 1:length(PNvsBLK_ratio_range)
    for pp = 1:length(BLK_powerratio_range)
        PNvsBLK_ratio = PNvsBLK_ratio_range(rr);
        BLK_powerratio = BLK_powerratio_range(pp);
        sinr = PNvsBLK_ratio*BLK_powerratio;
        samp_require(pp,rr) = (PNvsBLK_ratio)*(qfuncinv(pfa_goal) - sqrt(1+sinr+2*sinr^2)*qfuncinv(1-pm_goal))^2/sinr^2;
    end
end
%% plot
[X,Y] = meshgrid(BLK_powerratio_range_dB,PNvsBLK_ratio_range);
figure
contour(X,Y,log10(samp_require.'),'ShowText','on')
xlabel('Power Ratio (dB)')
ylabel('PN/BLK Ratio')
xlim([-30,-12])
title('Required Number of Samples (log10)')
