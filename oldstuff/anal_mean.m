% 03/07/2015
% Basically write matrix form into a summation form
% It gives us more insight how parameter N and R change the mean of
% statistics.

clear;clc;
R_range = 10:2:200;
N_range = [15,31,63];
N_num = length(N_range);
R_num = length(R_range);
for Rindex = 1:R_num
    for Nindex = 1:N_num
        R=R_range(Rindex);
        N=N_range(Nindex);
        bound=min(N,R);
        critical=0;
        for ii=1:bound
            critical = critical+(R-ii)*(N-ii);
        end
        mean_same(Rindex,Nindex) = (N+critical*2/R)/N^2;
        mean_diff(Rindex,Nindex) = (N-1/N*critical*2/R)/N^2;
    end
end

color=['b';'r';'c'];
figure;plot_setting();
for Nindex=1:N_num
    subplot(211)
    plot(R_range,10*log10(mean_same(:,Nindex)),color(Nindex,:));hold on
    subplot(212)
    plot(R_range,10*log10(mean_diff(:,Nindex)),color(Nindex,:));hold on
end

subplot(211)
grid on
legend('N = 15','N = 31','N = 63')
title('Mean of S_i^2 Anal.')
xlabel('R')
ylabel('Suppression Figure (dB)')

subplot(212)
grid on
legend('N = 15','N = 31','N = 63')
title('Mean of N_i^2 Anal.')
xlabel('R')
ylabel('Suppression Figure (dB)')
