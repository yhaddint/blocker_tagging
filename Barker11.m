clear all;
close all;
clc;

fs=32e9;
ts=1/fs;
XDelta=ts;


%%%%%%%%%%%%%%Start here to define x and y in tems of a time variable with
%%%%%%%%%%%%%%steps ts

% parameters for PN generation
rng(2);                                % Random seed
T_win      = 500e-6;                   % Time window to capture signal in [s], i.e., 500us 
sig_length = T_win * fs;               % Number of sample (conceptual sample in RF domain)
PN_num     = 2;                        % Num of inserted PNs
fchip      = 2e6;                      % PN rate in [Hz]
PN_OS      = fs/fchip;                 % Oversampling ratio for PN signal
PN_CFO     = 10e6;                     % offset b/w PNs and BLKs
fc         = [0.2,0.5]*1e9 + PN_CFO;   % PN center freq
t          = (0:sig_length-1).' * ts;  % s

% prepare PN per band
for bb = 1:PN_num
    Code(:,bb) = randi(2,sig_length/PN_OS,1) * 2 - 3; % PN code
    PN = kron( Code(:,bb),ones(PN_OS,1) );            % interpolate to requried BW
    x_real(:,bb) = PN .* ( cos(2*pi*t*fc(bb)) );      % upconvert to carrier
    y_imag(:,bb) = zeros(sig_length,1);               % imag part not used

end

% summation over all PNs
x = sum(x_real,2);
y = sum(y_imag,2);


%%%%%%%%%%%%%%%%%%%% End of changes 




Y=x+i*y;
save('data1.mat','Y','XDelta')



%%%%%Fig verification 

t=t/1e-9;

figure(1)
set(gcf,'color','w');   
subplot(2,1,1)
plot(t,x,'b','linewidth',2)
axis([0 220 -2 2])
xticks([0:20:220])

AX = gca;
set(AX,'fontsize',40);
set(AX,'fontweight','bold');
set(AX,'fontname','Times');
xlabel('Time (ns)','fontsize',50,'fontweight','bold')
ylabel(' Amplitude (V)','fontsize',50,'fontweight','bold','color','k')
set(gca,'LineWidth',5,'GridLineStyle',':')
grid on
box on



subplot(2,1,2)
plot(t,y,'b','linewidth',10)
axis([0 220 -2 2])
xticks([0:20:220])
AX = gca;
set(AX,'fontsize',40);
set(AX,'fontweight','bold');
set(AX,'fontname','Times');
xlabel('Time (ns)','fontsize',50,'fontweight','bold')
ylabel(' Amplitude (V)','fontsize',50,'fontweight','bold','color','k')
set(gca,'LineWidth',5,'GridLineStyle',':')
grid on
box on
