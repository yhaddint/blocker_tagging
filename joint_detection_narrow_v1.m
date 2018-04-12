%  04/12/2015
% test for sequential detection.
% square waveform is considered here.
% this easiest version only generate their gaussian statistics rather than
% generate waveform.

%clear;clc;clf;close all
warning off
blk_num=15;


M=1e2;
mus=0.287;
mun=0.0116;
sigmas=0.109;
sigman=0.3e-3;
sig_pow_dB=linspace(0,60,blk_num)'-90;
sig_pow_dB=sort(sig_pow_dB,'descend');
sig_pow=10.^(sig_pow_dB/10);


mu_matrix=ones(blk_num,blk_num)*mun;
mu_matrix=mu_matrix+diag(ones(1,blk_num)*(mus-mun));

sig_get = mu_matrix*sig_pow;
runtimes=1e3;
sig_noproc=zeros(blk_num,runtimes);
for ii=1:blk_num
    sigma_matrix=diag(ones(1,blk_num))*sigman;
    sigma_matrix(ii,ii)=sigmas;
    for jj=1:blk_num
        for kk=1:blk_num
            if (jj~=kk) 
               if jj==ii || kk==ii
                   sigma_matrix(jj,kk)=mus*mun;
               else
                   sigma_matrix(jj,kk)=mun*mun;
               end
            end
        end
    end
    sig_var(ii)=sig_pow'*sigma_matrix*sig_pow/M;
    sig_noproc(ii,:)=(sig_get(ii)+randn(1,runtimes)*sqrt(sig_var(ii)))/mus;
end

%% joint sequential detection 2
TH=-50;

diag_comp=1/(mus-mun);
offdiag_comp=-mun/(mus-mun)/(mus+(blk_num-1)*mun);
temp=ones(ii)*offdiag_comp+diag(ones(1,ii)*(diag_comp-offdiag_comp));

for rr=1:runtimes
    sig_proc2(:,rr) = temp*sig_noproc(:,rr)*mus;
end
%% plot
pd2 = sum(10*log10(sig_proc2')>TH)/runtimes;
figure(1)
plot(sig_pow_dB,pd2,'r');hold on
xlabel('raw power (dBm)')
ylabel('prob being detected')
%sum(10*log10(randn(1,runtimes)*offdiag_comp*sum(sqrt(sig_var)))>TH)/runtimes
%%
% figure
% plot_setting;
% plot(10*log10(sig_pow),'bx');hold on
% plot(10*log10(sum(sig_noproc)./mus/runtimes),'rx');hold on
% plot(10*log10(sum(sig_proc)/runtimes),'kx');hold on
% legend('raw power','calculated power','processed power')
%%

% tempp=[sig_pow(1:2);zeros(28,1)];
% 10*log10(tempp'*diag(diag(sigma_matrix))*tempp/M/mus^2)