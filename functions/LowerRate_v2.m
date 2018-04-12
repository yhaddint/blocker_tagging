function [ cal ] = LowerRate_v2( cal0,L )
%LOWERRATE_V1 Summary of this function goes here
%   Detailed explanation goes here
blocker_num=length(cal0(1,:));
N=length(cal0(:,1));

for ii=1:blocker_num
    % reuse PN sequences of length M*N
    cal(:,ii)=kron(ones(L/N,1),cal0(:,ii));
end



end
