function [ampl,freq] = freq_spectrum(ts,fs,mode)
if nargin<3
    mode = 'abs';
end
if nargin==1
    t = ts(:,1);
    y = ts(:,2);
    samplingRate = (length(t)-1) / (t(end) - t(1)); %Hz
elseif nargin>=2
    y = ts(:,1);
    samplingRate = fs;
end

% y = y - mean(y);

% nsamples = length(y);
% nfft = nsamples; %2^nextpow2(length(y));
% % w = hamming(length(y),'periodic');
% w = ones(size(y));
% 
% Et = sum(y.^2) / nsamples;
% 
% y_fft = fft(y.*w,nfft);
% A_rms = abs(y_fft(:)) / (nsamples); %two-sided rms amplitude
% A_pwr = abs(y_fft(:).^2) / (nsamples); %two-sided rms power
% 
% freq = linspace(0,samplingRate/2,nfft/2+1); 
% freq = freq(:);
% A_rms = A_rms(1:length(freq));
% A_rms(2:nfft/2) = sqrt(2.0) * A_rms(2:nfft/2); %频率 0~fs/2
% resp = A_rms(:);
% 
% Ef = sum(resp.^2);

signal = y(:);
window = rectwin(length(signal));
nfft = length(signal); %2^nextpow2(length(signal));
[pows,freq] = periodogram(signal,window,nfft,samplingRate,'onesided','power'); 
% [pows,freq] = pwelch(signal,length(signal),0,length(signal),samplingRate,'onesided','power'); ampl = (pows*2).^(1/2);
ampl_rms = pows.^(1/2); %rms amplitude
ampl_abs = (pows).^(1/2) * sqrt(2); %peak amplitude

% freq = freq(2:end-1);
% resp = resp(2:end-1);
if strcmp(mode,'abs')
    ampl = ampl_abs;
elseif strcmp(mode,'rms')
    ampl = ampl_rms;
elseif strcmp(mode,'dB')
    ampl = 20*log10(ampl_abs);
elseif strcmp(mode,'power')
    ampl = pows;
elseif strcmp(mode,'dBW')
    ampl = 10*log10(pows);
elseif strcmp(mode,'norm-power')
    ampl = pows / sum(pows);
elseif strcmp(mode,'norm-dBW')
    ampl = 10*log10(pows/sum(pows));
elseif strcmp(mode,'norm-abs')
    ampl = ampl_abs / sum(ampl_abs);
end