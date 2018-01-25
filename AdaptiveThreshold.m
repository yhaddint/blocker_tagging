%
%  07/23/2014
%  Stonger blockers will leak into weak one while vice versa is not that
%  important. So we tagging strongest one first, then modifies threshold
%  for other. 


clear;clc;clf;close all
warning off

% determined whether to apply estimation of power and detection
do_est=0;
do_det=0;
plot_intermediate=1;

symbols=5e4;
upsam=500;


N=400;
M=40;

%int_power_pool=0.5';

blocker_pow_pool=[1 2 4 8 16 32 64]/1024;


  

L=M*N*60;

for ii=1:3
    clear data
    hmod = modem.pskmod('M', 4, 'InputType', 'integer');
    hdesign  = fdesign.pulseshaping(upsam,'Square Root Raised Cosine');
    hpulse = design(hdesign);
    data = randi(4,symbols,1)-1;
    data = modulate(hmod, data);
    data = upsample(data,upsam);
    sig_bb(:,ii) = conv(data,hpulse.Numerator);
end

%% standard
%  We define Blocker with certain power, let's say 10, as standard.
%  It can be regarded as training signal to determine TH
for jj=1:length(blocker_pow_pool)
    
    blocker(1,:)=real(sig_bb(1:L,1))./sqrt(mean(real(sig_bb(1:L,1)).^2))*blocker_pow_pool(jj);
    cal0=PNgenerator_v1(M*N*4);

    % Equivalently lower rate of PN signal
    Dsampling=4;
    clearvars cal
    for ii=1:3
        % reuse PN sequences of length M*N
        cal(ii,:)=kron(ones(L/M/N/Dsampling,1),cal0(M*N*ii+1:M*N*(ii+1)));
        % Rate of PN is lower than Fs, which is equivalent to upsamling PN sequence
        cal_temp(ii,:)=kron(cal(ii,:),ones(1,Dsampling));
    end
    cal=cal_temp;
    clearvars cal_temp
        temp3=0.5*blocker(1,:).^2.*cal(1,:);
        temp4=blocker(1,:);
    raw_standard=temp3+temp4;
    
    for ii=1:fix(L/N/2)
        Nave_standard(ii)=abs(sum(raw_standard((ii-1)*N+1:ii*N)/N));
    end
    
    l_need=fix(length(Nave_standard)/M);
    MNave_standard=mean(reshape(Nave_standard(1:M*l_need),M,l_need));
    
    figure
    stem(sort(MNave_standard))
    title(['power =',num2str(blocker_pow_pool(jj))])
    mean(MNave_standar);

end