function [ out ] = sigma_s2_fun( N,R )
%SIGMA_S2_FUN Summary of this function goes here
%   Detailed explanation goes here
sigma=0;
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
        sigma=sigma+result*count_fun(N,a,b,c);    
        end
        end
    end
end

out=sigma/(N^4);
    
            

end

