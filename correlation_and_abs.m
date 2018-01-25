function [ sig_MNave ] = correlation_and_abs( sig,cal,N,M )

%  version 3 is energy detection while version 2 is absolute value
%  detection

L=length(sig);

sig_corr=sig.*cal;

%  N average
ll_need=L/N;
sig_Nave=abs(mean(reshape(sig_corr,N,ll_need)));

%  M average
l_need=fix(length(sig_Nave)/M);
sig_MNave=mean(reshape(sig_Nave(1:M*l_need),M,l_need));


% TH_standard=sorted_standard(fix(length(sorted_standard)*0.2));  

end

