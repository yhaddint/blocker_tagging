function [ x_freq, PSD_actual ] = get_PSD( input_waveform, FFT_size, channel_num )
%GET_PSD Summary of this function goes here
%   Detailed explanation goes here
if nargin==2
    channel_num = 100;
end

PSD_flip = 20*log10(abs(fft(input_waveform(1:FFT_size).*hann(FFT_size)))/FFT_size);
PSD_actual  = [PSD_flip(FFT_size/2+1:end);PSD_flip(1:FFT_size/2)];
x_freq = linspace(-channel_num/2,channel_num/2,FFT_size).';

end

