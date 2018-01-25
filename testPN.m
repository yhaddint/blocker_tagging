clear;clc
%%

%cal=PNgenerator_v2(100,10);

% corrMat=zeros(10,10);
% for ii=1:10
%     for jj=1:10
%         corrMat(ii,jj)=sum(cal(:,ii).*cal(:,jj));
%     end
% end

%%
cal0=PNgenerator_v4(110);
blocker_num=10;

for ii=1:blocker_num
cal(ii,:)=cal0(ii:99+ii);
end

corrMat=zeros(blocker_num,blocker_num);
for ii=1:blocker_num
    for jj=ii:blocker_num
        corrMat(ii,jj)=sum(cal(ii,:).*cal(jj,:));
    end
end
corrMat