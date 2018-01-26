function [ sig_MNave ] = correlation_and_pow( sig,cal,N,M )
%  version 3 is energy detection while version 2 is absolute value
%  detection

L=length(sig);
sig_corr=sig.*cal;
sig_Nave=(sum(reshape(sig_corr,N,(L/N)))/N).^2;
if M==1
    sig_MNave=sig_Nave;
else
    sig_MNave=sum(sig_Nave)/M;
end


end
