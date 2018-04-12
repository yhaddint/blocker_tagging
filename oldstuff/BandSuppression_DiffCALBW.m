%  10/01/2014
%  We use previous setting. For a given BLK, we tune BW ratio and parameter
%  N.

clear;clc;clf;close all

%  to compare fairly, we downsample to 5samples/cal period
downSampleTo=4;

N_pool=[7 15 31 63];%4:4:60;
% BLK_BW_pool=1;
BWratio_pool=4:4:60;

for nn=1:length(N_pool)
    N=N_pool(nn);
    SampleNumberAve=N*downSampleTo;
    %  Theoretically value of Rc: one for all diagonal component. -1/N for
    %  others,
%     Rc_initial=ones(N,N).*(-1/N)+diag(ones(N,1).*(1/N+1));
%     Cmatrix(:,:,1)=ones(SampleNumberAve,SampleNumberAve);
%     Cmatrix(:,:,2)=kron(Rc_initial,ones(downSampleTo,downSampleTo));
    Cmatrix=zeros(SampleNumberAve,SampleNumberAve,2);
    Cmatrix(:,:,1)=corrDecMat(N,downSampleTo,1);
    Cmatrix(:,:,2)=corrDecMat(N,downSampleTo,2);
    
    for ww=1:length(BWratio_pool)

        % envelop correlation matrix can be calculated in advance
        upsam=BWratio_pool(ww)*downSampleTo;
        hdesign  = fdesign.pulseshaping(fix(upsam),'Square Root Raised Cosine');
        hpulse = design(hdesign);
        rrc=hpulse.Numerator./sqrt(mean(hpulse.Numerator.^2));%./sqrt(length(hpulse.Numerator));
        for ii=1:length(rrc)
            h(ii)=rrc(ii:end)*rrc(1:end-ii+1)'/(length(rrc)-ii+1);
        end

        %%
        % Averaging length N and Samples to do such Averaging    

        cosCorr_1BW=cos(2*pi/upsam*(0:SampleNumberAve-1));
        %cosCorr_1MHz=cos(2*pi/fs*2*(0:SampleNumberAve-1));
        for ii=1:SampleNumberAve
            for jj=1:SampleNumberAve
                Rb(ii,jj)=h(1+abs(ii-jj));
                Rc_1BW(ii,jj)=cosCorr_1BW(abs(ii-jj)+1);
                %Rc_1MHz(ii,jj)=cosCorr_1MHz(abs(ii-jj)+1);
            end
        end

        powerBenchmark(ww)=SampleNumberAve^2;

        result_target(ww,nn)=sum(sum(Rb))/powerBenchmark(ww);
        result_1BW(ww,nn)=trace((Rc_1BW.*Rb)*Cmatrix(:,:,1))/powerBenchmark(ww);   
        %result_1MHz(nn,ss)=trace((Rc_1MHz.*Rb)*Cmatrix(:,:,1))/powerBenchmark(nn);   
        result_interf(ww,nn)=sum(sum(Rb))/trace((Rb)*Cmatrix(:,:,2));   
        result_interf_1BW(ww,nn)=trace((Rc_1BW.*Rb)*Cmatrix(:,:,2))/powerBenchmark(ww);   

    end
    clearvars Rb Rc_1BW Rc_1MHz
end
%%
clf
figure(18)
for nn=1:length(N_pool)
    plot(1./BWratio_pool*10,10*log10(result_target(:,nn)),'linewidth',2);
    hold on
end

title('CAL BW = 10MHz, gain for different BLK BW')
xlabel('BW of BLK (MHz)');
ylabel('Gain Ratio(dB)');
legend('N = 7','N = 15','N = 31','N = 63');
grid on
ylim([-10 0]);
xlim([0,1])
