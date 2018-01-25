function cal=PNgenerator_v1(N)


poly=[35 33 0];
order=poly(1);
initial_cond=[zeros(1,order-1) 1];

hpn = comm.PNSequence('Polynomial',poly, ...
           'SamplesPerFrame', N,'InitialConditions',initial_cond);
x1 = step(hpn);
cal=x1*2-1;
