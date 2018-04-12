function [ sig_out ] = get_rapp_square( sig_in,Vsat,P )
%GET_RAPP_NL Summary of this function goes here
%   Detailed explanation goes here
    
    sig_ratio = abs(sig_in)/Vsat;
    sig_out = (sig_in).^2./((1+sig_ratio.^(2*P)).^(1/P));

end

