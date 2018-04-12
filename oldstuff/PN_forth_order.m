clear;clc
N_range=[7 15 31 63];
for Nindex=1
N=N_range(Nindex);
cal0=PN_generator_downsample(N,N,1);
count=0;
cc=0;


for a=1:N
    for b=a:N
        for c=b:N
            for d=c:N
                count=count+1;
                result(count) = sum(cal0(a,:).*cal0(b,:).*cal0(c,:).*cal0(d,:));
                if result(count)==N
                    cc=cc+1;
                    PNcol(cc,:)=[a b c d];
                end
            end
        end
    end
end
%result=result/count;
figure
plot(sort(squeeze(result)))
% plot(1:count,ones(1,count)*N,'k--');
% plot(1:count,ones(1,count)*-N,'k--');
% plot(1:count,ones(1,count)*(-1/N),'k--');
%%
clearvars temp
temp=[PNcol(:,2)-PNcol(:,1),PNcol(:,3)-PNcol(:,2),PNcol(:,4)-PNcol(:,3)];
    if N==7
        keycomb7=unique(temp,'rows');
    elseif N==15
        keycomb15=unique(temp,'rows');
    elseif N==31
        keycomb31=unique(temp,'rows');
    elseif N==63
        keycomb63=unique(temp,'rows');
    end

end