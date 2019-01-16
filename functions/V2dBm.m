function [output_dBm] = V2dBm(inputV)
%V2DBM Summary of this function goes here
%   Detailed explanation goes here
if nargin==1
    res = 50;
elseif nargin==2
    res = resistance;
else
    error('[input_dBm, Resistance] as input');
end
    
input_W = inputV^2/(2*res);

output_dBm = 10*log10(input_W*1000);
end

