% 05/15/2015
% mean and variance of both S^2 and N^2 has been compared between anal. and
% sim. 
% Square waveform is considered here.
% 4-QAM, 16QAM, 64QAM are compared here.
% here random codes are used instead of PN codes.

% specially we care about
% E(b(n1)b(n2)b(n3)b(n4))-E(b(n1)b(n2))*E(b(n3)b(n4))
% for arbitary n1,n2,n3,n4

clear;
clc;%clf;close all
warning off

% Independent repeating experiment.
runtimes=5e4;

%N_range=[4];
%N_num = length(N_range);
M=50;
rrcPulseShape=0;
% BLK
R_range = 40;
R_num=length(R_range);
N_range=[15,31,63];
N_num=length(N_range);


var_typical = zeros(R_num,N_num);

SampleperSymbol=R_range;
SampleperFrame=M*SampleperSymbol;
%  How many BLK we want to tag simutaneously
tag_num=20;

%  How many active BLK existing. To get reasonable number of different
%  envelope, we generate 2x of them.
BLK_num=20;
BLKcan_num=20;

%  Define length of BLK we want.
stat_num=40;
L=SampleperFrame*stat_num;
clearvars sig_bb_r sig_bb_i

if rrcPulseShape==0
    clearvars unnormailzed
    for ii=1:BLKcan_num
        unnormalized_rec= kron(randi(4,M*stat_num,1)*2-5,ones(SampleperSymbol,1));
        sig_bb_r(:,ii) = unnormalized_rec./sqrt(mean(unnormalized_rec.^2));
        unnormalized_rec = kron(randi(4,M*stat_num,1)*2-5,ones(SampleperSymbol,1));
        sig_bb_i(:,ii) = unnormalized_rec./sqrt(mean(unnormalized_rec.^2));
        sig(:,ii) = (sig_bb_r(:,ii)+sig_bb_i(:,ii))/sqrt(2);
    end
else
    for ii=1:BLKcan_num
        upsam = SampleperSymbol;
        symbols = fix(L*2/upsam);

        clearvars data temp_data
        hmod = modem.qammod('M', 16, 'InputType', 'integer');
        %hmod = modem.pskmod('M', 16, 'InputType', 'integer');
        hdesign  = fdesign.pulseshaping(upsam,'Square Root Raised Cosine');
        hpulse = design(hdesign);
        data = randi(16,symbols,1)-1;
        data = modulate(hmod, data);
        data = upsample(data,upsam);
        temp_data = conv(data,hpulse.Numerator);
        unnormalized = real(temp_data(end-L+1:end));
        sig_bb_r(:,ii) = unnormalized./sqrt(sum(unnormalized.^2)/L);
        unnormalized = imag(temp_data(end-L+1:end));
        sig_bb_i(:,ii) = unnormalized./sqrt(sum(unnormalized.^2)/L);
        sig(:,ii)=(sig_bb_r(:,ii)+sig_bb_i(:,ii))/sqrt(2);
    end
end
%
for runindex=1:runtimes
    si=randi(L-SampleperFrame-50)+50;
    bi= randi(BLKcan_num);
    for n3=1:121
        n1=0;
        n2=n1;
        n4=n3;
        E_record_1(n3,runindex) = sig(si+n1,bi)*sig(si+n2,bi)*sig(si+n3-1,bi)*sig(si+n4-1,bi);
    end
end

for n3=1:121
    E4(n3) = mean(E_record_1(n3,:));
end
%%
% E4_smooth = zeros(1,length(E4));
% E4_smooth(1:2) = E4(1:2);
% E4_smooth(3:end) = (E4(1:end-2)+E4(2:end-1)+E4(3:end))/3;
xdata=(1:121)-1;
figure
plot_setting
plot(xdata(1:2:end),E4(1:2:end))
grid on
%% plot
% figure
% plot(var_typical);hold all