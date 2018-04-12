clear;clc
%SIGMA_S2_FUN Summary of this function goes here
%   Detailed explanation goes here
N=7;R=10;
sigma=0;
load('keycomb.mat');
if N==7
    keycomb=keycomb7;
    num=4;
elseif N==15
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

for a=0:N-1
    for c=0:N-1
        if a+c<N
        for b=0:N-a-c-1
            if ([a,b,c]==0)==[1 1 1]
                atten=1;
            elseif ([a,b,c]==0)==[1 0 1]
                atten=1;
            elseif max(sum((kron(ones(num,1),[a,b,c])==keycomb),2))==3
                atten=1;
            else
                atten=-1/N;
            end
        result=1;
        sigma=sigma+result*atten*count_fun(N,a,b,c);    
        end
        end
    end
end

out=sigma/(N^4)
   
%%
%cal0=ones(N,N);
cal0=PN_generator_downsample(N,N,1);
count=0;
for CAL_index1=1:N;
    for CAL_index2=1:N
        if CAL_index1 ~= CAL_index2
            count=count+1;
            result(count)=sum(cal0(:,CAL_index1).*cal0(:,CAL_index2))/N;
        end
    end
end

fourth_sim = var(result)+mean(result)^2
       
%%
sigma=0;
for x1=1:N
    for x2=1:N
        for x3=1:N
            for x4=1:N
                temp=sort([x1 x2 x3 x4]);
                a=temp(2)-temp(1);
                b=temp(3)-temp(2);
                c=temp(4)-temp(3);
            if ([a,b,c]==0)==[1 1 1]
                atten=1;
            elseif ([a,b,c]==0)==[1 0 1]
                atten=1;
            elseif max(sum((kron(ones(num,1),[a,b,c])==keycomb),2))==3
                atten=1;
            else
                atten=-1/N;
            end
        result=1;
        sigma=sigma+result*atten;    
        end
        end
    end
end
sigma/(N^4)

