%  05/14/2015
% test for joint estimation.
% RRC waveform is considered here.(square waveform is in version 2)
% this version uses real waveform.

runtimes_sim = 5e2;
runtimes_analy = 1e4;

blk_num = 1;
M = 1e2;
N = 43;
SampleNumberAve = 20;
stat_num = 3e3;
P = M*N;
L = M*SampleNumberAve*stat_num;

for ii = 1:blk_num
    upsam = SampleNumberAve;
    symbols = fix(L*2/upsam);
    
    clearvars data temp_data
    hmod = modem.pskmod('M', 4, 'InputType', 'integer');
    hdesign  = fdesign.pulseshaping(upsam,'Square Root Raised Cosine');
    hpulse = design(hdesign);
    data = randi(4,symbols,1)-1;
    data = modulate(hmod, data);
    data = upsample(data,upsam);
    temp_data = conv(data,hpulse.Numerator);
    unnormalized = real(temp_data(end-L+1:end));
    sig_r(:,ii) = unnormalized./sqrt(sum(unnormalized.^2)/L);
    unnormalized = imag(temp_data(end-L+1:end));
    sig_i(:,ii) = unnormalized./sqrt(sum(unnormalized.^2)/L);
    sig_i(:,ii) = zeros(size(sig_r(:,ii)));
    sig(:,ii) = (sig_r(:,ii)+sig_i(:,ii));

end
%%
E2 = mean(sig.^2);
E4 = mean(sig.^4);
value_to_record = zeros(N,N,runtimes_sim);

for runindex=1:runtimes_sim
    start_point=randi(length(sig(:,1)))-N-1;
    for nn=1:N
        for kk=1:N
            value_to_record(runindex,nn,kk) = sig(start_point+nn).^2*sig(start_point+kk).^4;
        end
    end
end
E_matrix=zeros(N,N);
for nn=1:N
    for kk=1:N
        E_matrix(nn,kk) = sum(squeeze(value_to_record(:,nn,kk)))/runtimes_sim-E2*E4;
    end
end

sum(sum(E_matrix>0))