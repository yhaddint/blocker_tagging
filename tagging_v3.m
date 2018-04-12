function [ sig_MNave ] = tagging_v3( sig,cal,N,M )

L=length(sig);
sig_corr=sig.*cal;
%% Energy
sig_Nave=(sum(reshape(sig_corr,N,(L/N)))/N).^2;
%% absolute value
%temp=sum(reshape(sig_corr,N,(L/N)))/N;
%sig_Nave=(temp>0).*temp+(temp<0).*(-temp);

if M==1
    sig_MNave=sig_Nave;
else
    sig_MNave=sum(sig_Nave)/M;
end


end
