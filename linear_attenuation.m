function [output_signal,impulse_response,h,f] = linear_attenuation(beta,depth,Ts,input_signal)

%% input arguments
Ts = Ts * 1e6; %microseconds
NFFT = 11; NPT = 2^NFFT; N2 = NPT / 2;
%% generate attenuation frequency responses
[impulse_response,h,f] = frqrsp(beta,depth,Ts,NFFT);

if nargin>3 && ~isempty(input_signal)
    output_signal = conv(impulse_response(:),input_signal(:));
    output_signal = output_signal(1:length(input_signal));
    output_signal = output_signal(:);
else
    output_signal = [];
end

% %% attenuation model
% frd0 = idfrd(h(1:N2),f(1:N2)*1e6,Ts*1e-6,'FrequencyUnit','Hz','TimeUnit','seconds'); 
% sys0 = tfest(frd0,5,4,'Ts',Ts*1e-6);
% %% model simulation
% if nargin>3
%     tsig = input_signal(:,1);
%     usig = input_signal(:,2);
%     [yout,tout] = lsim(sys0,usig,tsig);
%     output_signal = [tout,yout];
% else
%     output_signal = [];
% end
end

function [impulse_response,h,f] = frqrsp(beta,depth,Ts,NFFT)
% firstly generate the log-magnitude transfer function from mag_slope and 
% sampling rate; secondly generate the phase characteristic from hilbert 
% transform
NPT = 2^NFFT; %NFFT=10;
N2 = NPT / 2;
N2P1 = N2 + 1;
NPTP2 = NPT + 2;

f = (0:NPT-1) / NPT / Ts; %MHz
xmag = zeros(1,NPT); %log(H(w))
for ndx=2:N2P1
    frqMhz = (ndx-1) / NPT / Ts; %MHz
    xmag(ndx) = -beta * frqMhz * depth / 20; %log-mag function
    xmag(NPTP2-ndx) = xmag(ndx); % even function of f
end
xphs = dht(xmag,NPT); %hilbert transform

h = zeros(1,NPT);
for ndx=2:N2P1
    xm = 10.^(xmag(ndx)); %mag from log-mag
    h(ndx) = xm * (cos(xphs(ndx)) + 1i * sin(xphs(ndx))); 
    h(NPTP2-ndx) = xm * (cos(xphs(ndx)) - 1i * sin(xphs(ndx))); 
end

impulse_response = real(ifft(h,NPT)); %fft for inverse
end


function [xphs] = dht(xmag,NPT)

NPTP2 = NPT+2;
xphs = zeros(1,NPT);

for ndx=2:fix(NPT/2)
    omega = (ndx - 1) * 2.0 * pi / NPT;
    % calculate the [-pi, pi] integral of log(|H(w)|)cot((theta-omega)/2) d(theta)
    % d(theta) = 2 * pi / NPT
    sum = 0;
    for k=1:NPT
        if k==ndx
            continue;
        else
            theta = (k - 1) * 2.0 * pi / NPT;
            cotan = cos((theta - omega) / 2.0) / sin((theta - omega) / 2.0);
            xm = log(10.^(xmag(k)));
            sum = sum + xm * cotan * (2.0 * pi / NPT);
        end
    end
    
    xphs(ndx) = sum / (2.0 * pi);
    xphs(NPTP2-ndx) = -xphs(ndx); % phase is odd fcn
end

end