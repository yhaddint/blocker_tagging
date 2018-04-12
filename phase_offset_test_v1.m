clear;clc
N=5;
for ii=1:1e3
symbol_re = randi(2,1,N)*2-3;
symbol_im = randi(2,1,N)*2-3;
phi = rand*2*pi;
%phase_offset(ii) = 0.8478;
amp1(ii)=norm(symbol_re*(cos(phi)+sin(phi))+symbol_im*(-sin(phi)+cos(phi)));
amp2(ii)=norm(symbol_re+symbol_im);
end
figure
plot(amp1,'r');hold on
plot(amp2,'b')
%%
% figure
% plot(symbol_re*cos(phase_offset(ii))-symbol_im*sin(phase_offset(ii)));hold on
% plot(symbol_re+symbol_im);hold on