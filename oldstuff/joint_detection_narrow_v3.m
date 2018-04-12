%  04/14/2015
% test for sequential detection.
% RRC waveform is considered here.(square waveform is in version 2)
% this version uses real waveform
% 
% Sensitivity of BW mismatch.
clear;clc;clf;close all


runtimes_sim = 2e4;
runtimes_analy = 2e4;

blk_num = 3;
M = 1;
N = 63;
SampleNumberAve = 20;
stat_num = 1e4;
P = M*N;
L = M*SampleNumberAve*stat_num;

sig_pow_dB = [-30,-40,-90];
sig_pow_dB = sort(sig_pow_dB,'descend');
sig_pow = 10.^(sig_pow_dB/10);

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
    sig(:,ii) = (sig_r(:,ii)+sig_i(:,ii))*sqrt(sig_pow(ii));

end

%% calibration sequences
cal0=PNgenerator_v5(N,N,1);
CAL=LowerRate_v2(cal0,P);

%% process matrix
% N = 63, R = 20
mus=0.3039;
mun=0.0113;
sigmas=0.2256;
sigman=0.235e-3;

diag_comp=1/(mus-mun);
offdiag_comp=-mun/(mus-mun)/(mus+(blk_num-1)*mun);
matrix_proc=ones(ii)*offdiag_comp+diag(ones(1,ii)*(diag_comp-offdiag_comp));

%% tagging (sim.1 results)
for runindex=1:runtimes_sim
    CAL=CAL(:,randperm(N));
    start_point = randi(L-P-1,blk_num);
    for ii=1:blk_num
        sig_cache(:,ii) = sig(start_point(ii):start_point(ii)+P-1,ii);
    end
    sig_mix = zeros(P,1);

    for ii=1:blk_num
        sig_mix = sig_mix + sig_cache(:,ii).*CAL(:,ii);
    end
    
    for ii=1:blk_num
        result_noproc_sim(runindex,ii) = tagging_v3(sig_mix,CAL(:,ii),N,1);
    end
    
    result_noproc_sim_lookin(runindex,1) = tagging_v3(sig_cache(:,1),ones(P,1),N,1);
    result_noproc_sim_lookin(runindex,2) = tagging_v3(sig_cache(:,1).*CAL(:,1),CAL(:,3),N,1);
    result_noproc_sim_lookin(runindex,3) = tagging_v3(sig_cache(:,2),ones(P,1),N,1);
    result_noproc_sim_lookin(runindex,4) = tagging_v3(sig_cache(:,2).*CAL(:,2),CAL(:,3),N,1);
    result_noproc_sim_lookin(runindex,5) = tagging_v3(sig_cache(:,3),ones(P,1),N,1);
    result_noproc_sim_lookin(runindex,6) = tagging_v3(sig_cache(:,3).*CAL(:,3),CAL(:,1),N,1);
 
    result_noproc_sim(runindex,3) = result_noproc_sim_lookin(runindex,2)+result_noproc_sim_lookin(runindex,4)+result_noproc_sim_lookin(runindex,5);
    
    result_proc_sim(runindex,:) = (matrix_proc*result_noproc_sim(runindex,:)')';
end

%%
term1 = mean(result_noproc_sim_lookin(:,1));
term2 = mean(result_noproc_sim_lookin(:,2));
term3 = mean(result_noproc_sim_lookin(:,3));
term4 = mean(result_noproc_sim_lookin(:,4));
term5 = mean(result_noproc_sim_lookin(:,5));
term6 = mean(result_noproc_sim_lookin(:,6));

sr(1) = term1/sig_pow(1);
or(1) = term2/sig_pow(1);
sr(2) = term3/sig_pow(2);
or(2) = term4/sig_pow(2);
sr(3) = term5/sig_pow(3);
or(3) = term6/sig_pow(3);

matrix_iso = [sr(1),or(2),or(3);or(1),sr(2),or(3);or(1),or(2),sr(3)];
step_num = 11;

result_robust=zeros(step_num,blk_num);
%%
for bb=1:3
    matrix_iso = [sr(1),or(2),or(3);or(1),sr(2),or(3);or(1),or(2),sr(3)];
    sr_range=linspace(1.2*sr(bb),0.8*sr(bb),step_num);
    or_range=linspace(0.8*or(bb),1.2*or(bb),step_num);
    for ss=1:step_num
        for kk=1:blk_num
            if kk==bb
                matrix_iso(kk,bb)=sr_range(ss);
            else
                matrix_iso(kk,bb)=or_range(ss);
            end
        end
        matrix_proc2 = inv(matrix_iso);
        for runindex=1:runtimes_sim
            temp = (matrix_proc2*result_noproc_sim(runindex,:)')';
            result_proc_sim(runindex,ss) = temp(blk_num);
        end
        result_robust(ss,bb) = mean(result_proc_sim(:,ss));
    end
end

%%
figure
plot_setting
plot(linspace(-0.2,0.2,step_num),10*log10(result_robust),'o-');hold on
plot(linspace(-0.2,0.2,step_num),ones(1,step_num)*10*log10(mean(result_noproc_sim(:,3))),'k--');
xlabel('BW Mismatch')
ylabel('Mean of Esimated Power (dBm)')
legend('Mismatch with Interferer 1 (-30 dBm)','Mismatch with Interferer 2 (-40 dBm)','Mismatch with Interferer 3 (-90 dBm)','No Calibartion (Ref.)')
