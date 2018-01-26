function [ output ] = count_fun(N,a,b,c)
%COUNT_FUN Summary of this function goes here
%   Detailed explanation goes here
choose=bi2de([a,b,c]>0,'left-msb');
count(1)=nchoosek(N,1);
count(2)=(N-c)*nchoosek(4,1);
count(3)=(N-b)*nchoosek(4,2);
count(4)=(N-b-c)*3*4;
count(5)=(N-a)*nchoosek(4,1);
count(6)=(N-a-c)*3*4;
count(7)=(N-a-b)*3*4;
count(8)=(N-a-b-c)*factorial(4);
output=count(choose+1);

end

