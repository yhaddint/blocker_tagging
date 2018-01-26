%  09/02/2014
%  Tune CAL BW. BLK BW is set to be 1MHz
%  Different sampling time.


clear;clc
figIndex=99;

%  to compare fairly, we downsample to 5samples/cal period
downSampleTo=5;

% BW of blocker is a constant (1MHz)
BW_blk_pool=1;

% We change BW of calibration signals by chaning oversample rate
oversamp=[5 10 30];
fs=300;

% Averaging length N and Samples to do such Averaging
N=31;
SampleNumberAve=N*downSampleTo;

%% LPF
% How far a Cutoff is from BW of Calibration.
d=fdesign.lowpass('Fp,Fst,Ap,Ast',1/downSampleTo,1.25/downSampleTo,1,60);
Hd = design(d,'equiripple');
ll=(length(Hd.Numerator)+1)/2;

if length(Hd.Numerator)<=SampleNumberAve
    Lrow=[zeros(1,SampleNumberAve-length(Hd.Numerator)) Hd.Numerator];
else
    Lrow=[Hd.Numerator(1:SampleNumberAve)];
end

Lcol=[Hd.Numerator';zeros(SampleNumberAve-1,1)];
Lmatrix_full=toeplitz(Lcol,Lrow);
Lmatrix=Lmatrix_full(ll:end-ll+1,:);

%
%Spulse=cos((0:SampleNumberAve-1)./(SampleNumberAve-1)*2*pi*2)';
% Spulse=[1;zeros((SampleNumberAve-1),1)];
% figure;
% subplot(311)
% stem(Spulse);
% subplot(312)
% stem(Lmatrix*Spulse);
% subplot(313)
% stem(conv(Hd.Numerator',Spulse))


%% cal
cal0_temp=PN_generator(N);
cal0=zeros(N,N);
for ii=1:N
    cal0(:,ii)=[cal0_temp(ii:end);cal0_temp(1:ii-1)];
end

hdesign  = fdesign.pulseshaping(SampleNumberAve,'Square Root Raised Cosine');
hpulse = design(hdesign);
c=(length(hpulse.Numerator)+1)/2;
% hindex=zeros(7,oversamp);
% for kk=1:oversamp
%     hindex(:,kk)=((kk-1)*N+(c-3*N*oversamp:N*oversamp:c+3*N*oversamp))';
% end

cal=zeros(SampleNumberAve,N);
for ii=1:N
    cal(:,ii)=kron(cal0(:,ii),ones(downSampleTo,1));
end

% tempmat=zeros(N,oversamp);
% for ii=1:N
%     for kk=1:oversamp
%         temp=cconv(cal0(:,ii),hpulse.Numerator(hindex(:,kk))');
%         tempmat(:,kk)=temp(4:end-3)*31*oversamp;
%     end
%     cal_ni(:,ii)=reshape(tempmat',N*oversamp,1);
% end
clearvars temp

Cmatrix=zeros(SampleNumberAve,SampleNumberAve,4);
% 1. Neighbor without LPF
% Cmatrix(:,:,1)=ones(N*oversamp,N*oversamp);

% We want to consider averaging effect over ensemble of calibration
aveCalNum=(N-1)*N/2;

for ii=1:N
    for jj=ii+1:N
        % 2. Interference without LPF
        Cmatrix(:,:,2)=Cmatrix(:,:,2)+(cal(:,ii)*cal(:,ii)').*(cal(:,jj)*cal(:,jj)')/aveCalNum;
        % 4. Interference with LPF
        Cmatrix(:,:,4)=Cmatrix(:,:,4)+(Lmatrix'*cal(:,ii)*cal(:,ii)'*Lmatrix).*(cal(:,jj)*cal(:,jj)')/aveCalNum;
    end
    % 3. Neighbor with LPF
    Cmatrix(:,:,3)=Cmatrix(:,:,3)+(Lmatrix'*cal(:,ii)*cal(:,ii)'*Lmatrix).*(cal(:,ii)*cal(:,ii)')/N;
    Cmatrix(:,:,1)=Cmatrix(:,:,1)+(cal(:,ii)*cal(:,ii)').*(cal(:,ii)*cal(:,ii)')/N;
end

%%
for kk=1:length(oversamp)
BW_cal=fs/oversamp(kk);

% supression
%  1 Neighbor Suppression without LPF
%  2 Interference Suppression without LPF
%  3 Neighbor Suppression with LPF
%  4 Interference Suppression with LPF

%

upsam=fs/BW_blk_pool/(oversamp(kk)/downSampleTo);

hdesign  = fdesign.pulseshaping(fix(upsam),'Square Root Raised Cosine');
hpulse = design(hdesign);
rrc=hpulse.Numerator./sqrt(mean(hpulse.Numerator.^2));%./sqrt(length(hpulse.Numerator));
for ii=1:length(rrc)
    h(ii)=rrc(ii:end)*rrc(1:end-ii+1)';
end

% In plot x label is divided by plotxNumber
plotxNumber=30;
% In plot x label is up to plotxLim MHz
plotxLim=30;

df=linspace(0,plotxLim,plotxNumber);

for ii=1:SampleNumberAve
    for jj=1:SampleNumberAve
        Rb(ii,jj)=h(1+abs(ii-jj));
    end
end

power_benchmark=sum(sum(Rb));

for findex=1:plotxNumber
    cosCorr=cos(2*pi*df(findex)/fs*(oversamp(kk)/downSampleTo)*(0:SampleNumberAve-1));
    for ii=1:SampleNumberAve
        for jj=1:SampleNumberAve
            Rc(ii,jj)=cosCorr(abs(ii-jj)+1);
            Rb(ii,jj)=h(1+abs(ii-jj));
        end
    end



    for suppressionChoice=1:4
        result(findex,kk,suppressionChoice)=trace((Rc.*Rb)*Cmatrix(:,:,suppressionChoice))/power_benchmark;   
    end
end
end


%% plot
color=['b  ';'r  ';'g  ';'b-o';'r-o';'g-o'];
for suppressionChoice=1:4
for kk=1:length(oversamp)
    
    plotx=[-fliplr(df) df(2:end)];
    ploty=10*log10(result(:,kk,suppressionChoice))';
    ploty=[fliplr(ploty) ploty(2:end)];
   
    figure(figIndex)
    subplot(2,2,suppressionChoice)
    plot(plotx,ploty,color(kk,:),'linewidth',2);hold on
    
end
figure(figIndex)
legend('CAL BW = 60MHz','CAL BW = 30MHz','CAL BW =	10MHz');
ylim([-40,10]);  
xlabel('Distance to Center (MHz)');
ylabel('Suppression (dB)');
grid on
switch suppressionChoice
    case 1
        title('Neighbor Suppression without LPF');
    case 2
        title('Interference Suppression without LPF');
    case 3
        title('Neighbor Suppression with LPF');
    case 4
        title('Interference Suppression with LPF');
end

end

figIndex=figIndex+1;