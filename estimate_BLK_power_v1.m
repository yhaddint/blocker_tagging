% 03/09/2015
% use optimization method to reconstrunct power of BLK
clear;clc
BLK_num=10;
N=2;
Power_BLK=ones(N,1)*0.5;
Power_BLK(2)=1;
mus=0.8288;
mun=0.006077;
M=10;
sigmas2=0.1036/M;
sigman2=0.0002854/M;
U=ones(N)*mun+diag(ones(1,N))*(mus-mun);
SIGMA_measured=randn(N)*sigman2+diag(randn(1,N)*sigmas2);
SIGMA=ones(N)*sigman2+diag(ones(1,N)*sigmas2);
R=(U+SIGMA_measured)*Power_BLK;



result=zeros(200,200);
for col=1:200
    for row=1:200
        P=[row/100;col/100];
        for ii=1:N
            result(row,col)=result(row,col)+(U(ii,:)*P-R(ii)).^2/(2*SIGMA(ii,:)*(P.*P))+...
                log((SIGMA(ii,:)*(P.*P)));
        end
    end
end
%%
% cvx_begin
%     variable P(2)
%     minimize (U(1,:)*P-R(1)).^2/(2*SIGMA(1,:)*P)+log(ones(1,2)*(SIGMA*P))+(U(2,:)*P-R(2)).^2/(2*SIGMA(2,:)*P)
% cvx_end



%% plot
figure
X=1:200;
Y=1:200;
%v=exp(0:1:50).*result(10,100);
contour(X,Y,result,100);

