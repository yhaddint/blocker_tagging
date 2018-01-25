% 04/24/2014
% This is to observe autocorrelation property of PN sequences generated.
% Basically, PN sequence is correlated with shifted version of itself. It
% seems good.

clear;clc;clf;close all
L=4e5;
cal=PNgenerator_v1(L);
auto(1)=sum(cal.^2)/L;
for ii=2:1e3
auto(ii)=sum(cal(1:end-(ii-1)).*cal(ii:end))/(L-(ii-1));
end
%% plot
plot(auto,'linewidth',2);
title('autocorr')
xlabel('shifted sample');
grid on
xlim([-50,1050])