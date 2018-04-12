clear;clc;clf;close all
N = 1e3;
runtimes=1e3;
a=randn(1,N);
result = zeros(runtimes,2);
for ii=1:runtimes
    result(ii,1) = a*(randi(2,N,1)*2-3);
    result(ii,2) = a*(randi(2,N,1)*2-3);
end
figure
plot(result(:,1),result(:,2),'.')
figure
plot(randn(1,N),randn(1,N),'.r')