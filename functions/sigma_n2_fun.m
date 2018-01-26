function [ out ] = sigma_n2_fun( N,R )
%SIGMA_S2_FUN Summary of this function goes here
%   Detailed explanation goes here
sigma=0;
load('keycomb.mat');
if N==15
    keycomb=keycomb15;
    num=28;
elseif N==31
    keycomb=keycomb31;
    num=140;
elseif N==63
    keycomb=keycomb63;
    num=620;
end


if N>=R
    bound =R-1;
else 
    bound =N-1;
end

for a=0:bound
    part0=R-a;
    for c=0:bound
        if a+c<N
        for b=0:N-a-c-1
            kb=fix((a+b)/R)+1;
            kc=fix((a+b+c)/R)+1;
            if ([a,b,c]==0)==[1 1 1]
                atten=1;
            elseif ([a,b,c]==0)==[1 0 1]
                atten=1;
            elseif max(sum((kron(ones(num,1),[a,b,c])==keycomb),2))==3
                atten=1;
            else
                atten=-1/N;
            end
            if (kc==kb)
                part1=kc*R-(a+b+c);
                part2=c;
                if part0<=part1
                    result=part0/R;
                elseif part0<=part1+part2
                    result=part1/R;
                else
                    result=(part1+(part0-part1-part2))/R;
                end
            else
                part1=kb*R-(a+b);
                part2=R-c;
                if part0<=part1
                    result=0;
                elseif part0<=part1+part2
                    result=(part0-part1)/R;
                else
                    result=part1/R;
                end
            end
        sigma=sigma+result*atten*count_fun(N,a,b,c);    
        end
        end
    end
end

out=sigma/(N^4);
    
            

end

