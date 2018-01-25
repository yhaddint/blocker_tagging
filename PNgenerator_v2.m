function cal=PNgenerator_v2(N,blocker_num)


poly1=[7 3 0];
poly2=[7 3 2 1 0];
order=poly1(1);
initial_cond=[zeros(1,order-1) 1];

hgld = comm.GoldSequence('FirstPolynomial',poly1,...
             'SecondPolynomial', poly2,...
             'FirstInitialConditions', initial_cond,...
             'SecondInitialConditions', initial_cond,...
             'Index', 4, 'SamplesPerFrame', N);
for ii=1:blocker_num
cal(:,ii) = step(hgld);
cal(:,ii)=cal(:,ii).*2-1;
end

%cal=x1*2-1;
