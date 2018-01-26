%  09/08/2014
%  This program is to compare performance of different CAL/BLK BW selection. 
%  In this one I don't want to consider LPF.
%  Maybe add it in next version.

clear;clc;clf;close all
warning('off')

BWratio_pool=10:2:100;

for bb=1:length(BWratio_pool)
upsam=BWratio_pool(bb);
N=upsam;

hdesign  = fdesign.pulseshaping(upsam,'Square Root Raised Cosine');
hpulse = design(hdesign);
rrc=hpulse.Numerator./sqrt(mean(hpulse.Numerator.^2));%./sqrt(length(hpulse.Numerator));
for ii=1:length(rrc)
    h(ii)=rrc(ii:end)*rrc(1:end-ii+1)';
end


for ii=1:upsam
    for jj=1:upsam
        Rb(ii,jj)=h(1+abs(ii-jj));
    end
end


power_benchmark=sum(sum(Rb));


% for findex=1:plotxNumber
deltaFreq=[0 1 2 4];
for ff=1:length(deltaFreq)
    cosCorr(:,ff)=cos(2*pi*deltaFreq(ff)/N*(0:N-1));
    for ii=1:N
        for jj=1:N
            Rc(ii,jj,ff)=cosCorr(abs(ii-jj)+1,ff);
        end
    end
end


Cmatrix(:,:,1)=-1/N*ones(N,N);
Cmatrix(:,:,1)=Cmatrix(:,:,1)+diag(ones(1,N)*((N+1)/N));
Cmatrix(:,:,2)=ones(N,N);

for suppressionChoice=1:2
    for ff=1:length(deltaFreq)
        result(bb,ff,suppressionChoice)=trace((Rc(:,:,ff).*Rb)*Cmatrix(:,:,suppressionChoice))/power_benchmark;   
    end
end
clearvars Cmatrix cosCorr Rc Rb
end


%% plot
color=['b  ';'r  ';'g  ';'c  ';'r-o';'g-o'];
figure
for ff=1:length(deltaFreq)
plot(BWratio_pool,10*log10(result(:,ff,1)),color(ff,:),'linewidth',2);
hold on
end
grid on
legend('0MHz','0.5MHz','1MHz','2MHz')