function [input_V] = dBm2V(input_dBm,resistance)
%DBM2MV Summary of this function goes here
%   Detailed explanation goes here

if nargin==1
    res = 50;
elseif nargin==2
    res = resistance;
else
    error('[input_dBm, Resistance] as input');
end
    
input_dBW = input_dBm - 30;
input_W = 10^((input_dBW)/10);
input_V = sqrt(input_W*2*res);

% sig_pow = 10.^((inputArg1-17.5)./10); % the fixed offset is due to 50 om resistance
end

