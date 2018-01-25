function [ cal ] = LowerRate_v1( cal0,Dsampling,L )
%LOWERRATE_V1 Summary of this function goes here
%   Detailed explanation goes here
blocker_num=length(cal0(1,:));
N=length(cal0(:,1));

for ii=1:blocker_num
    % reuse PN sequences of length M*N
    cal_temp(:,ii)=kron(ones(L/N/Dsampling,1),cal0(:,ii));
    % Rate of PN is lower than Fs, which is equivalent to upsamling PN sequence
    cal(:,ii)=kron(cal_temp(:,ii),ones(Dsampling,1));
end



end

