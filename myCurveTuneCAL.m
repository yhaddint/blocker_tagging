%  09/02/2014
%  Tune CAL BW. BLK BW is set to be 1MHz
%  Different sampling time.


clear;clc;clf;close all
warning('off')

% BW of blocker is a constant (1MHz)
BW_blk_pool=1;

% We change BW of calibration signals by chaning oversample rate
cal_BW=[30 15 60];
fs=360;
downSampleTo=fs/cal_BW(3);

% Averaging length N and Samples to do such Averaging
N=63;
SampleNumberAve=N*downSampleTo;

for kk=1:3
clearvars -except result cal_BW fs downSampleTo N SampleNumberAve kk power_benchmark    
    
%% LPF
cutoff=cal_BW(kk)/fs;
% How far a Cutoff is from BW of Calibration.
d=fdesign.lowpass('Fp,Fst,Ap,Ast',cutoff,1.25*cutoff,1,60);
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
if kk==1
    blk_num=31;
elseif kk==2
    blk_num=15;
else
    blk_num=63;
end

cal0_temp=PN_generator(blk_num);
cal0=zeros(blk_num,blk_num);
for ii=1:blk_num
    cal0(:,ii)=[cal0_temp(ii:end);cal0_temp(1:ii-1)];
end

hdesign  = fdesign.pulseshaping(SampleNumberAve,'Square Root Raised Cosine');
hpulse = design(hdesign);
c=(length(hpulse.Numerator)+1)/2;
% hindex=zeros(7,oversamp);
% for kk=1:oversamp
%     hindex(:,kk)=((kk-1)*N+(c-3*N*oversamp:N*oversamp:c+3*N*oversamp))';
% end

cal=zeros(SampleNumberAve,blk_num);
for ii=1:blk_num
    temp=kron(cal0(:,ii),ones(downSampleTo*cal_BW(3)/cal_BW(kk),1));
    cal(:,ii)=[temp;temp(1:SampleNumberAve-length(temp))];
    clearvars temp
end

% tempmat=zeros(N,oversamp);
% for ii=1:N
%     for kk=1:oversamp
%         temp=cconv(cal0(:,ii),hpulse.Numerator(hindex(:,kk))');
%         tempmat(:,kk)=temp(4:end-3)*31*oversamp;
%     end
%     cal_ni(:,ii)=reshape(tempmat',N*oversamp,1);
% end


Cmatrix=zeros(SampleNumberAve,SampleNumberAve,4);
% 1. Neighbor without LPF
% Cmatrix(:,:,1)=ones(N*oversamp,N*oversamp);

% We want to consider averaging effect over ensemble of calibration
aveCalNum=(blk_num-1)*blk_num/2;

for ii=1:blk_num
    for jj=ii+1:blk_num
        % 2. Interference without LPF
        Cmatrix(:,:,2)=Cmatrix(:,:,2)+(cal(:,ii)*cal(:,ii)').*(cal(:,jj)*cal(:,jj)')/aveCalNum;
        % 4. Interference with LPF
        Cmatrix(:,:,4)=Cmatrix(:,:,4)+(Lmatrix'*cal(:,ii)*cal(:,ii)'*Lmatrix).*(cal(:,jj)*cal(:,jj)')/aveCalNum;
    end
    % 3. Neighbor with LPF
    Cmatrix(:,:,3)=Cmatrix(:,:,3)+(Lmatrix'*cal(:,ii)*cal(:,ii)'*Lmatrix).*(cal(:,ii)*cal(:,ii)')/N;
    Cmatrix(:,:,1)=Cmatrix(:,:,1)+(cal(:,ii)*cal(:,ii)').*(cal(:,ii)*cal(:,ii)')/N;
end


% supression
%  1 Neighbor Suppression without LPF
%  2 Interference Suppression without LPF
%  3 Neighbor Suppression with LPF
%  4 Interference Suppression with LPF

%

upsam=fs*cal_BW(kk)/cal_BW(3);

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

if kk==1
power_benchmark=sum(sum(Rb));
end

    for findex=1:plotxNumber
        cosCorr=cos(2*pi*df(findex)/fs*(cal_BW(3)/cal_BW(kk))*(0:SampleNumberAve-1));
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
for kk=1:3
    
    plotx=[-fliplr(df) df(2:end)];
    ploty=10*log10(result(:,kk,suppressionChoice))';
    ploty=[fliplr(ploty) ploty(2:end)];
   
    figure(99)
    subplot(2,2,suppressionChoice)
    plot(plotx,ploty,color(kk,:),'linewidth',2);hold on
    
end
figure(99)
legend('CAL BW = 30MHz','CAL BW = 15MHz','CAL BW = 60MHz');
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
