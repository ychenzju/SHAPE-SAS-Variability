function [Ts,Ys] = gnrtPulseSignal_x(fc,n,repeats,samplingRate,base,mode,shape,lead_n,cease_n)
% generate the gaussian signal
% fc : the central frequency; 4MHz
% n:   the number of periods of central frequency; 8 periods
% samplingRate: the sampling frequency, 50MHz
% Author: Yao Chen, 2021-12-27
% email: ychen.zju@outlook.com
if nargin>7
    augmented = 1;
end
% fc = 2.5e6;
% 
% n = [1, 2, 1; %magnitude
%     1,0.5,1;  %frequency
%     1,0,1];   %phase
% repeats = 1;
% samplingRate = 100e6;
% base = 'sin';
% mode = 'bipolar';
% shape = 'step';
% lead_n = 1;
% cease_n = 12;

global plot_on

[rows, cols] = size(n);
if rows==1 && cols==1
    period = 1/fc * n;
    ts = (0:1/samplingRate:period)';
    ts = ts(1:end-1);
    if strcmp('sin',base)
        xs = (sin(2*pi*fc*ts));
    elseif strcmp('cos',base)
        xs = (cos(2*pi*fc*ts));
    elseif strcmp('rect',base)
        xs = sign(sin(2*pi*fc*ts));
    elseif strcmp('tri',base)
        xs = sin(2*pi*(2*1-1)*fc*ts) .* (8 * (-1)^(1+1) /(2*1-1)^2/pi^2) ...
            +sin(2*pi*(2*2-1)*fc*ts) .* (8 * (-1)^(2+1) /(2*2-1)^2/pi^2) ...
            +sin(2*pi*(2*3-1)*fc*ts) .* (8 * (-1)^(3+1) /(2*3-1)^2/pi^2) ...
            +sin(2*pi*(2*4-1)*fc*ts) .* (8 * (-1)^(4+1) /(2*4-1)^2/pi^2) ...
            +sin(2*pi*(2*5-1)*fc*ts) .* (8 * (-1)^(5+1) /(2*5-1)^2/pi^2) ...
            +sin(2*pi*(2*6-1)*fc*ts) .* (8 * (-1)^(5+1) /(2*6-1)^2/pi^2) ...
            +sin(2*pi*(2*7-1)*fc*ts) .* (8 * (-1)^(5+1) /(2*7-1)^2/pi^2);
    elseif strcmp('general',base)
        xs = (1-cos(2*pi*fc*ts/n)).*...
            cos(2*pi*fc*ts) / 2;
    end
    zs = (xs + 1)/2;
elseif rows==3 && cols>1
    ts = [];
    xs = [];
    zs = [];
    for col_ndx=1:cols
        mag_times = n(1,col_ndx);
        frq_times = n(2,col_ndx);
        phs = n(3,col_ndx) * pi;
        
        frq = fc * frq_times;
        period = 1 / frq * 0.5; %half period
        ti = (0:1/samplingRate:period)';
        ti = ti(2:end);
        
        if strcmp('sin',base)
            xi = (sin(2*pi*frq*ti + phs));
        elseif strcmp('cos',base)
            xi = (cos(2*pi*frq*ti + phs));
        elseif strcmp('rect',base)
            xi = sign(sin(2*pi*fc*ti + phs));
        elseif strcmp('tri',base)
            xi = sin(2*pi*(2*1-1)*frq*ti + phs) .* (8 * (-1)^(1+1) /(2*1-1)^2/pi^2) ...
                +sin(2*pi*(2*2-1)*frq*ti + phs) .* (8 * (-1)^(2+1) /(2*2-1)^2/pi^2) ...
                +sin(2*pi*(2*3-1)*frq*ti + phs) .* (8 * (-1)^(3+1) /(2*3-1)^2/pi^2) ...
                +sin(2*pi*(2*4-1)*frq*ti + phs) .* (8 * (-1)^(4+1) /(2*4-1)^2/pi^2) ...
                +sin(2*pi*(2*5-1)*frq*ti + phs) .* (8 * (-1)^(5+1) /(2*5-1)^2/pi^2) ...
                +sin(2*pi*(2*6-1)*frq*ti + phs) .* (8 * (-1)^(5+1) /(2*6-1)^2/pi^2) ...
                +sin(2*pi*(2*7-1)*frq*ti + phs) .* (8 * (-1)^(5+1) /(2*7-1)^2/pi^2);
        end
        zi = (xi + 1)/2;
        xi = xi * mag_times;
        zi = zi * mag_times;
        
        if ~isempty(ts)
            ti = ti + ts(end) ;
        end
        
        ts = cat(1,ts,ti);
        xs = cat(1,xs,xi);
        zs = cat(1,zs,zi);
    end
    
end

switch shape
    case 'gauss'
        if strcmp('unipolar',mode)
            ys = gausswin(length(ts)) .* zs;
        elseif strcmp('bipolar',mode)
            ys = gausswin(length(ts)) .* xs;
        end
    case 'tukey_p25'
        if strcmp('unipolar',mode)
            ys = tukeywin(length(ts)) .* zs;
        elseif strcmp('bipolar',mode)
            ys = tukeywin(length(ts),0.25) .* xs;
        end
    case 'tukey_p50'
        if strcmp('unipolar',mode)
            ys = tukeywin(length(ts)) .* zs;
        elseif strcmp('bipolar',mode)
            ys = tukeywin(length(ts),0.50) .* xs;
        end
    case 'tukey_p75'
        if strcmp('unipolar',mode)
            ys = tukeywin(length(ts)) .* zs;
        elseif strcmp('bipolar',mode)
            ys = tukeywin(length(ts),0.75) .* xs;
        end
    case 'hann'
        if strcmp('unipolar',mode)
            ys = tukeywin(length(ts),1) .* zs;
        elseif strcmp('bipolar',mode)
            ys = tukeywin(length(ts),1) .* xs;
        end
    case 'sine'
        if strcmp('unipolar',mode)
            ys = sin(2*pi*fc/n/2*ts) .* zs;
        elseif strcmp('bipolar',mode)
            ys = sin(2*pi*fc/n/2*ts) .* xs;
        end
    case 'step'
        if strcmp('unipolar',mode)
            ys = zs;
        elseif strcmp('bipolar',mode)
            ys = xs;
        end
end

if augmented==1 && (lead_n>=1 || cease_n>=1)
    if lead_n>=1
        lead_ts = fliplr(0:-1/samplingRate:-lead_n/fc);
        lead_ts = lead_ts';
        ts = [lead_ts(1:end-1); ts];
        ys = [zeros(length(lead_ts)-1,1); ys];
    end
    if cease_n>=1
        cease_ts = (0:1/samplingRate:cease_n/fc)';
        ts = [ts; cease_ts(2:end)+ts(end)];
        ys = [ys; zeros(length(cease_ts)-1,1)];
    end

    if lead_n>=1
        ts = ts - lead_ts(1);
    end
    
%     if strcmp('unipolar',mode)
%         switch shape
%             case 'gauss'
%                 ys = [zeros(length(lead_ts)-1,1)+0.5; ys; zeros(length(cease_ts)-1,1)+0.5];
%             case 'sin'
%                 ys = [zeros(length(lead_ts)-1,1)+0.5; ys; zeros(length(cease_ts)-1,1)+0.5];
%             case 'step'
%                 ys = [zeros(length(lead_ts)-1,1); ys; zeros(length(cease_ts)-1,1)];
%         end
%     else
%         ys = [zeros(length(lead_ts)-1,1); ys; zeros(length(cease_ts)-1,1)];
%     end
end

Ts = ts;
Ys = ys;

if repeats>1
    for rep = 2:repeats
        Ts = [Ts; ts(2:end)+ts(end)*(rep-1)];
        Ys = [Ys; ys(2:end)];
    end
end

if plot_on==1
    nsamples = length(ys);
    ys_fft = fft(ys); ys_fft=ys_fft(:);
    ys_fft_fspan = linspace(0,samplingRate/2,nsamples/2+1); ys_fft_fspan=ys_fft_fspan(:);
    ys_ssfft = ys_fft(1:length(ys_fft_fspan));
    ys_ssfft(2:end-1) = 2.0 * ys_ssfft(2:end-1);
    
    figure; 
    abs_ys_ssfft = abs(ys_ssfft);
    subplot(211), plot(Ts*1e6,Ys); xlabel('Time (us)'); ylabel('Signal Amplitude'); title('Generated Pulse Signal'); grid on
    subplot(212), plot(ys_fft_fspan*1e-6,abs_ys_ssfft); xlabel('Freq (MHz)'); ylabel('Magnitude (AU)'); grid on
    xlim([0,3*fc/1e6])
end

