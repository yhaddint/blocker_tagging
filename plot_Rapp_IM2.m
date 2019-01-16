clear;clc
input_range_dBm = -20:10;
resistance = 50;
fs = 1e9;
fc = 100e6;
T_window = 1e-3;
t_sample = transpose( linspace(0, T_window, (T_window * fs)) );


% ---- Analog Anti-aliasing Filter -------
fcutoff = 0.015;                          % Fcutoff is around 15MHz
fcutoff_DSP = fcutoff/fs*2;               % Effective "digital" frequency in RF sim.
AA_taps = 2000;                           % Num of Tap (just for RF sim purpose)
bhi_AA = fir1(AA_taps,fcutoff_DSP,'low'); % Fcutoff is around 15MHz


Asat = 0.5; % Parameter to control saturation input [volt]
P_rapp = 1; % Parameter to control IIP2

for IP_pow_idx = 1:length(input_range_dBm)
    
    % Print status
    clc
    fprintf( 'Input power trial %d out of %d\n',IP_pow_idx, length(input_range_dBm) );
    
    % Input tone power scaling
    input_V = dBm2V( input_range_dBm(IP_pow_idx), resistance );
    input_tone = sin( 2 * pi * fc * t_sample ) * input_V;
    
    % Rapp model for square-law LNA
    output_tone = get_rapp_square( input_tone, Asat, P_rapp );
    
    % Filter tone in DC
    AA_shift = AA_taps/2+1; % Delay after AA filter
    temp_out = filter(bhi_AA,1,output_tone);
    output_tone = temp_out(AA_shift:end);
    
    % Convert waveform to dBm
    output_rms = sqrt( mean(abs(output_tone).^2) );
    output_dBm(IP_pow_idx) = V2dBm( output_rms );
    
end

%%
figure
plot(input_range_dBm, output_dBm,'o-');hold on
xlabel('Input Power [dBm]')
ylabel('Output Power [dBm]')
grid on
xlim([-20,10])
