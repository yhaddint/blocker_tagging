function [result] = corrDecMat(N,downSampleTo,option)
%CORRDECMAT_FUN_V1 Summary of this function goes here
%   Detailed explanation goes here

SamN=N*downSampleTo;
corrDecMat=zeros(SamN,SamN);
cal0=PN_generator_downsample(N,N,downSampleTo);
result=zeros(SamN,SamN);
if option==2
    for ii=1:N
        for jj=1:N
            if ii==jj
                result=result;
            else
                result=result+((cal0(:,ii).*cal0(:,jj))*(cal0(:,ii).*cal0(:,jj))')/(N*(N-1));
            end
        end
    end

elseif option==1
    for ii=1:N
        for jj=ii
            result=result+((cal0(:,ii).^2)*(cal0(:,jj).^2)')/N;
        end
    end
end
    

end

