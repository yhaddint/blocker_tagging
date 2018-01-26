clear;clc
N=3;
shouldbe=N^4
% count(1)=nchoosek(N,1);
% count(2)=nchoosek(N,2)*nchoosek(4,1);
% count(3)=nchoosek(N,2)*nchoosek(4,2);
% count(4)=nchoosek(N,3)*3*4;
% count(5)=nchoosek(N,2)*nchoosek(4,1);
% count(6)=nchoosek(N,3)*3*4;
% count(7)=nchoosek(N,3)*3*4;
% count(8)=nchoosek(N,4)*factorial(4);
% sum(count);

test=zeros(N,N,N);
for a=0:N-1
    for b=0:N-1-a
        for c=0:N-1-a-b
            test(a+1,b+1,c+1)=count_fun(N,a,b,c);
        end
    end
end
sum(sum(sum(test)))

% test2=zeros(N,N,N);
% for x1=1:N
%     for x2=1:N
%         for x3=1:N
%             for x4=1:N
%                 y=sort([x1,x2,x3,x4]);
%                 test2(y(2)-y(1)+1,y(3)-y(2)+1,y(4)-y(3)+1)...
%                     =test2(y(2)-y(1)+1,y(3)-y(2)+1,y(4)-y(3)+1)+1;
%             end
%         end
%     end
% end
                
