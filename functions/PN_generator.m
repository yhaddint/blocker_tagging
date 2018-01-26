function [ cal ] = PN_generator( N )
%PN_GENERATOR Summary of this function goes here
% PN_generator using LFSR. It is named PN_generator_v3 in the previous
% folder
% Input: N, The intended length of PN sequence. It can be from <15, 31, 63>
% Output: cal, the PN sequence with length N

    if N==31
        poly=[5 3 0];
    elseif N==15
        poly=[4 3 0];
    elseif N==63
        poly=[6 5 0];
    end
    order=poly(1);
    initial_cond=[zeros(1,order-1) 1];

    hpn = comm.PNSequence('Polynomial',poly, ...
               'SamplesPerFrame', N,'InitialConditions',initial_cond);
    x1 = step(hpn);
    cal=x1*2-1;

end

