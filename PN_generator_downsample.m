function cal_downsample=PN_generator_downsample(blk_num,N,downSampleTo)

cal=zeros(N,blk_num);
cal_initial = zeros(blk_num,blk_num);
cal_downsample=zeros(blk_num*downSampleTo,blk_num);

switch blk_num
    case{7}
        poly=[3 2 0];
    case{15}
        poly=[4 3 0];
    case{31}
        poly=[5 3 0];
    case{63}
        poly=[6 5 0];
    case{127}
        poly=[7 6 0];
    case{255}
        poly=[8 6 5 4 0];
end

order=poly(1);
initial_cond=[zeros(1,order-1) 1];

hpn = comm.PNSequence('Polynomial',poly, ...
           'SamplesPerFrame', blk_num,'InitialConditions',initial_cond);
x1 = step(hpn);
cal_single=x1*2-1;

for ii=1:blk_num
    cal_initial(:,ii)=[cal_single(ii:end);cal_single(1:ii-1)];
end

replicate_num = fix(N/blk_num)+1;
cal_temp = kron(ones(replicate_num,1),cal_initial);


cal=cal_temp(1:N,:);
%  Now consider samples not exactly on the peak of calibration signals.

if downSampleTo==1
    cal_downsample=cal;
else
    hmod = modem.pskmod('M', 2, 'InputType', 'integer');
    hdesign  = fdesign.pulseshaping(downSampleTo,'Square Root Raised Cosine');
    hpulse = design(hdesign);
    for nn=1:blk_num
        data = upsample(cal(:,nn),downSampleTo);
        hpulse.Numerator=hpulse.Numerator./max(hpulse.Numerator);
        temp_data = conv(data,hpulse.Numerator);
        ll=length(temp_data);
        sl=length(cal(:,1))*downSampleTo;
        cal_downsample(:,nn)=temp_data(ll/2-sl/2+1:ll/2+sl/2)';
    end
end