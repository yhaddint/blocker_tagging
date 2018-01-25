function [ sig_MNsorted ] = tagging_v1( sig,cal,N,M )
%TAGGING_V1 Summary of this function goes here
%   Detailed explanation goes here

L=length(sig);

sig_corr=sig.*cal;

for ii=1:fix(L/N/2)
    sig_Nave(ii)=abs(sum(sig_corr((ii-1)*N+1:ii*N)/N));
end

l_need=fix(length(sig_Nave)/M);

sig_MNave=mean(reshape(sig_Nave(1:M*l_need),M,l_need));

sig_MNsorted=sort(sig_MNave);
% TH_standard=sorted_standard(fix(length(sorted_standard)*0.2));  

end

