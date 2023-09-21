function [data_ft,frequency_out,magnitude_out] = make_freq_spect(data,fs,bool_plot)
% Make_freq_spect performs a fourier analysis with the matlab function fft
% and plots the frequency spectrum of the signal
%
%   Input Arguments:
%           - Data = vector with input signal
%           - fs = sampling frequency of the input signal
%           - bool_plot= A Boolean, when 1 = figure of the freq. spectrum
%                                   when 0 = no figure
%
%   Output Arguments:
%           - data_ft = output from the fourier analysis
%           - frequency_out = x-axis of the freq. plot.(freq. harmonics)
%           - magnitude_out = y-axis of the freq. plot (amplitude harmonics)
%
%   --------Author: Maarten Afschrift (14/03/2014)-------

% Discrete Fourier transform.
data_ft=fft(data); 
sze = length(data); 
ff= fix(sze/2) + 1; 
f = [0:ff-1]*fs/sze; 

% get the frequency 
frequency_out=f(1:ff)';
% get the magnitude
magnitude_out=abs(data_ft(1:ff)/sze*2);

if bool_plot
    figure
    plot(frequency_out, magnitude_out);  
    xlabel('Frequency in Hz');
    ylabel('Magnitude');
    title('Frequency spectrum');
    axis tight;
end

end

