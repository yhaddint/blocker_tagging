function [ output ] = anal_mean_fun( R,N,option)
% function [ output ] = anal_mean_fun( R,N,option)
% option: 't' calculate mean of S^2 (target)
%         'n' calculate mean of N^2 (nontarget)


bound=min(N,R);
critical=0;
for ii=1:bound
    critical = critical+(R-ii)*(N-ii);
end
mean_same = (N+critical*2/R)/N^2;
mean_diff = (N-1/N*critical*2/R)/N^2;

if option == 't'
    output=mean_same;
elseif option == 'n'
    output=mean_diff;
end


end

