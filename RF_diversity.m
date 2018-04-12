clear;clc;
N = 1024;
SNR = 6;
noise_ampl = sqrt(10^(-SNR/10));
MCtimes = 5e3;
for MCindex = 1:MCtimes
    n1 = randn(N,1)*noise_ampl;
    n2 = randn(N,1)*noise_ampl;
    sig = randn(N,1);
    for ii= 1:10
        len = 1*2^ii;
        result_H0(ii,MCindex) = n1(1:len)'*n2(1:len)/len;
        result_H1(ii,MCindex) = (sig(1:len)+n1(1:len))'*(sig(1:len)+n2(1:len))/len;
    end
end
for ii=1:10
    temp = sort(result_H0(ii,:),'ascend');
    TH = temp(MCtimes*0.9);
    PD(ii) = sum(result_H1(ii,:)>TH)/MCtimes;
end
%
figure
semilogx(1*2.^(1:10),PD)
grid on
ylim([0,1])
%%
N = 5120;
SNR = -4;
noise_ampl = sqrt(10^(-SNR/10));
MCtimes = 5e3;
PD = zeros(100,1);
for MCindex = 1:MCtimes
    n1 = randn(N,1)*noise_ampl;
%     n2 = randn(N,1)*noise_ampl;
    sig = randn(N,1);
    for ii= 1:100
        len = ii*50;
        result_H0(ii,MCindex) = n1(1:len)'*n1(1:len)/len;
        result_H1(ii,MCindex) = (sig(1:len)+n1(1:len))'*(sig(1:len)+n1(1:len))/len;
    end
end
for ii=1:100
    temp = sort(result_H0(ii,:),'ascend');
    TH = temp(MCtimes*0.9);
    PD(ii) = sum(result_H1(ii,:)>TH)/MCtimes;
end
%
figure
semilogx(50*(1:100),PD)
grid on
ylim([0,1])