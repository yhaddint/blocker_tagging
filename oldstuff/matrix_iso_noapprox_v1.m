%  04/14/2015
% test for sequential detection.
% RRC waveform is considered here.(square waveform is in version 2)
% this version uses real waveform.
clear;clc;clf;close all


runtimes_sim = 1e4;
runtimes_analy = 1e4;

blk_num = 63;
M = 1;
N = 63;
SampleNumberAve = 20;
stat_num = 3e3;
P = M*N;
L = M*SampleNumberAve*stat_num;

sig_pow_dB = zeros(1,blk_num);
sig_pow_dB = sort(sig_pow_dB,'descend');
sig_pow = 10.^(sig_pow_dB/10);


upsam = SampleNumberAve;
symbols = fix(L*2/upsam);
hmod = modem.pskmod('M', 4, 'InputType', 'integer');
hdesign  = fdesign.pulseshaping(upsam,'Square Root Raised Cosine');
hpulse = design(hdesign);

for ii = 1:blk_num
    
    clearvars data temp_data
    data = randi(4,symbols,1)-1;
    data = modulate(hmod, data);
    data = upsample(data,upsam);
    temp_data = conv(data,hpulse.Numerator);
    unnormalized = real(temp_data(end-L+1:end));
    sig_r(:,ii) = unnormalized./sqrt(sum(unnormalized.^2)/L);
    unnormalized = imag(temp_data(end-L+1:end));
    sig_i(:,ii) = unnormalized./sqrt(sum(unnormalized.^2)/L);
    sig_i(:,ii) = zeros(size(sig_r(:,ii)));
    sig(:,ii) = (sig_r(:,ii)+sig_i(:,ii))*sqrt(sig_pow(ii));

end

%% calibration sequences
cal0=PNgenerator_v5(N,N,1);
CAL=LowerRate_v2(cal0,P);

%% test begin
result_noproc_sim=zeros(runtimes_sim,N,N);
for runindex=1:runtimes_sim
    
    if mod(runindex,100) ==0
        runindex/100
    end
    
    start_point = randi(L-P-1,blk_num);
    for ii=1:blk_num
        sig_cache(:,ii) = sig(start_point(ii):start_point(ii)+P-1,ii);
    end
    
    for rr=1:N
        sig_mix = sig_cache(:,rr).*CAL(:,rr);
        for cc=rr:N
            result_noproc_sim(runindex,rr,cc) = tagging_v3(sig_mix,CAL(:,cc),N,M);
        end
    end
end
%% organize coeeficient
for rr=1:N
    for cc=rr:N
        mean_coeff(rr,cc) = mean(squeeze(result_noproc_sim(:,rr,cc)));
        var_coeff(rr,cc) = var(squeeze(result_noproc_sim(:,rr,cc)));
    end
end
%%
matrix_iso_mean = mean_coeff+mean_coeff'-diag(diag(mean_coeff));
Z = 10*log10(matrix_iso_mean);
X=1:N;
Y=1:N;
figure
imagesc(X,Y,Z)
colorbar
matrix_proc_sep = inv(matrix_iso_mean);
%%
matrix_iso_var = var_coeff+var_coeff'-diag(diag(var_coeff));
Z = 10*log10(matrix_iso_var);
X=1:N;
Y=1:N;
figure
imagesc(X,Y,Z)
colorbar
%%
mean_comp = sort(reshape(matrix_iso_mean,1,63*63));
mean_offdiag_comp = mean_comp(1:end-63);
figure;
plot(mean_offdiag_comp)
%%
A = diag(ones(1,63)*(mean(diag(matrix_iso_mean))-mean(mean_offdiag_comp)));
B = ones(63,63)*mean(mean_offdiag_comp);
Z = 10*log10(A+B);
X=1:N;
Y=1:N;
figure
imagesc(X,Y,Z)
colorbar