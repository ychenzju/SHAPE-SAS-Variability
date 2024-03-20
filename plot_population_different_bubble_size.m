close all;
clearvars;
%% 显示参数
colors = {'b','r','k','g','m','c'};
styles = {'-','-.','--',':','-','-.'};
mrkers = {'s','^','v','o','*','d'};
widths = {1.0,1.0,1.0,1.0};
total_number = 4096;
min_number = 0;
while min_number<total_number
    colors = cat(2,colors,colors);
    styles = cat(2,styles,styles);
    mrkers = cat(2,mrkers,mrkers);
    widths = cat(2,widths,widths);
    min_number = min([length(colors),length(styles),length(mrkers),length(widths)]);
end
colors = colors(1:total_number);
styles = styles(1:total_number);
mrkers = mrkers(1:total_number);
widths = widths(1:total_number);

%% plot signals
bubble_sizes = [0,1,2,3,4,5];
caselab = {};
SHpsd = zeros(6,length(bubble_sizes));
SHdem = zeros(6,length(bubble_sizes));
sh = SHdem;
sh_pi = sh;
sh_am = sh;
sh_cps = sh;
for bsiz_ndx = 1:length(bubble_sizes)
    rad = bubble_sizes(bsiz_ndx);
    fname = 'population_003_Povs';
    if abs(rad-0)<eps
        caselab{bsiz_ndx} = ['Random Radii'];
    else
        fname = [fname,'_R',strrep(num2str(rad),'.','p'),'.mat'];
        caselab{bsiz_ndx} = ['Radii = ',num2str(rad), '\mum'];
    end
    load(fname);
    
    fig = figure(1001); fig.Position = [700,100,1000,1000];
    for pov_ndx = 1 : length(povs)
        if pov_ndx == 1
            figure(1001);
            subplot(length(bubble_sizes),3,(bsiz_ndx-1)*3+1);
            X = Treceived{pov_ndx} * 1e6;
            Y = Preceived{pov_ndx};
            plot(X,Y); grid on; hold on;
            ylim([-150,150]);
            YL = ylim();
            rectangle('Position',[sg_window_beg*1e6,YL(1),(sg_window_end-sg_window_beg)*1e6,YL(2)-YL(1)],'LineStyle',':','EdgeColor','g','LineWidth',1.4);
            xlabel('Time (\mus)');
            ylabel({caselab{bsiz_ndx},'Pressure (Pa)'});
            title('Receiving Acoustic Signal');
            xlim([X(1),X(end)]);
            
            ax = subplot(length(bubble_sizes),3,(bsiz_ndx-1)*3+2);
            Tt = Treceived{pov_ndx};
            Pt = Preceived{pov_ndx};
            ROI_LIST = (Tt>=sg_window_beg) & (Tt<=sg_window_end);
            signal = Pt(ROI_LIST);
            [resp,freq] = freq_spectrum(Pt(ROI_LIST),probe_fs,'abs');
            
            % [pows,freq] = periodogram(signal,rectwin(length(signal)),length(signal),probe_fs,'onesided','power'); resp = (pows*2).^(1/2);
            % [pows,freq] = pwelch(signal,length(signal),0,length(signal),probe_fs,'onesided','power'); resp = (pows*2).^(1/2);
            % relation between power spectrum density and power spectrum: psd = pow / (fs * N);
            X = freq(2:end) / 1e6; Y = 20*log10(resp(2:end));
            plot(X,Y,'-o','MarkerSize',3); grid on; hold on;
            ylim([-40,40]);
            YL = ylim();
            %rectangle('Position',[frqs(frq_ndx)/2/1e6-0.2,YL(1),0.4,YL(2)-YL(1)],'LineStyle',':','EdgeColor','r','LineWidth',1.4);
            xlabel('Frequency (MHz)');
            ylabel('Amplitude(dB)');
            title('PSD-derived Amplitude Spectrum');
            xlim([0,7.5]);
            ax.XTick = 0:frqs(frq_ndx)/2/1e6:7.5;
            
            subplot(length(bubble_sizes),3,(bsiz_ndx-1)*3+3);
            SHA = rf2iq_filter(Pt,probe_fs,Frq/2);
            ROI_LIST = (Tt>=sg_window_beg) & (Tt<=sg_window_end);
            sh_meanvalue = mean(SHA(ROI_LIST));
            X = Tt * 1e6; Y = 20*log10(SHA);
            plot(X,Y); grid on; hold on;
            ylim([-40,40]);
            YL = ylim();
            rectangle('Position',[sg_window_beg*1e6,YL(1),(sg_window_end-sg_window_beg)*1e6,YL(2)-YL(1)],'LineStyle',':','EdgeColor','g','LineWidth',1.4);
            %line([sg_window_beg;sg_window_end]*1e6,[sh_meanvalue;sh_meanvalue],'Color','r','LineStyle',':','LineWidth',1.4);
            xlabel('Time (\mus)');
            ylabel('Amplitude (dB)');
            title('IQ-demodulated Subharmonic Signal');
            xlim([X(1),X(end)]);
        end
        
        %% calculate the subharmonic amplitude vs. ambient pressure
        Tt = Treceived{pov_ndx};
        Pt = Preceived{pov_ndx};
        
        ROI_LIST = (Tt>=sg_window_beg) & (Tt<=sg_window_end);
        signal = Pt(ROI_LIST);
        [resp,freq] = freq_spectrum(signal,probe_fs,'dBW');
        sh_fc = frqs(frq_ndx)/2;
        sh_bw = 0.2e6;
        FOCUS_LIST = freq>=(sh_fc-sh_bw) & freq<=(sh_fc+sh_bw);
        SHpsd(pov_ndx,bsiz_ndx) = mean(resp(FOCUS_LIST));

        SHA = rf2iq_filter(Pt,probe_fs,Frq/2);
        ROI_LIST = (Tt>=sg_window_beg) & (Tt<=sg_window_end);
        SHdem(pov_ndx,bsiz_ndx) = mean(20*log10(SHA(ROI_LIST)));
        sh(pov_ndx,bsiz_ndx) = mean(20*log10(SHA(ROI_LIST)));

        Pt = PreceivedPI{pov_ndx};
        SHA = rf2iq_filter(Pt,probe_fs,Frq/2);
        ROI_LIST = (Tt>=sg_window_beg) & (Tt<=sg_window_end);
        sh_pi(pov_ndx,bsiz_ndx) = mean(20*log10(SHA(ROI_LIST)));

        Pt = PreceivedAM{pov_ndx};
        SHA = rf2iq_filter(Pt,probe_fs,Frq/2);
        ROI_LIST = (Tt>=sg_window_beg) & (Tt<=sg_window_end);
        sh_am(pov_ndx,bsiz_ndx) = mean(20*log10(SHA(ROI_LIST)));

        Pt = PreceivedCPS{pov_ndx};
        SHA = rf2iq_filter(Pt,probe_fs,Frq/2);
        ROI_LIST = (Tt>=sg_window_beg) & (Tt<=sg_window_end);
        sh_cps(pov_ndx,bsiz_ndx) = mean(20*log10(SHA(ROI_LIST)));
        
    end
    
end

%% 显示PSD和IQ-Demodulation分析的次谐波响应-环境压力的关系
fig = figure(1002); fig.Position = [200,200,1000,800];
subplot(2,2,1);
X = povs / 1e3;
Y = SHpsd;
plt = plot(X,Y); grid on;
legend(caselab,'Location','southwest');
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
end
xlabel('Ambient Pressure (kPa)');
ylabel('Amplitude (dB)');
title('PSD-derived Subharmonic Response');
subplot(2,2,2);
X = povs / 1e3;
Y = SHpsd-ones(size(SHpsd,1),1)*SHpsd(1,:);
plt = plot(X,Y); grid on;
legend(caselab,'Location','southwest');
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
end
xlabel('Ambient Pressure (kPa)');
ylabel('Amplitude (dB)');
title('PSD-derived Subharmonic Response (Relative)');
subplot(2,2,3);
X = povs / 1e3;
Y = SHdem;
plt = plot(X,Y); grid on;
legend(caselab,'Location','southwest');
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
end
xlabel('Ambient Pressure (kPa)');
ylabel('Amplitude (dB)');
title('IQ-demodulated Subharmonic Response');
subplot(2,2,4);
X = povs / 1e3;
Y = SHdem-ones(size(SHdem,1),1)*SHdem(1,:);
plt = plot(X,Y); grid on;
legend(caselab,'Location','southwest');
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
end
xlabel('Ambient Pressure (kPa)');
ylabel('Amplitude (dB)');
title('IQ-demodulated Subharmonic Response (Relative)');

%% 显示不同造影模式下IQ-Demodulation分析的次谐波响应（相对变化）
fig = figure(1004); fig.Position = [100,100,1000,800];
subplot(2,2,1);
X = povs / 1e3;
Y = sh - ones(size(sh,1),1)*sh(1,:);
plt = plot(X,Y); grid on;
legend(caselab,'Location','southwest');
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
end
xlabel('Ambient Pressure (kPa)');
ylabel('Amplitude (dB)');
title('Subharmonic Response in Linear Imaging');
subplot(2,2,2);
X = povs / 1e3;
Y = sh_pi - ones(size(sh_pi,1),1)*sh_pi(1,:);
plt = plot(X,Y); grid on;
legend(caselab,'Location','southwest');
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
end
xlabel('Ambient Pressure (kPa)');
ylabel('Amplitude (dB)');
title('Subharmonic Response in PI Imaging');
subplot(2,2,3)
X = povs / 1e3;
Y = sh_am - ones(size(sh_am,1),1)*sh_am(1,:);
plt = plot(X,Y); grid on;
legend(caselab,'Location','southwest');
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
end
xlabel('Ambient Pressure (kPa)');
ylabel('Amplitude (dB)');
title('Subharmonic Response in AM Imaging');
subplot(2,2,4)
X = povs / 1e3;
Y = sh_cps - ones(size(sh_cps,1),1)*sh_cps(1,:);
plt = plot(X,Y); grid on;
legend(caselab,'Location','southwest');
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
end
xlabel('Ambient Pressure (kPa)');
ylabel('Amplitude (dB)');
title('Subharmonic Response in CPS Imaging');

%% 显示不同造影模式下IQ-Demodulation分析的次谐波响应（绝对变化）
fig = figure(1005); fig.Position = [200,200,1000,800];
subplot(2,2,1);
X = povs / 1e3;
Y = sh;
plt = plot(X,Y); grid on;
legend(caselab,'Location','southwest');
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
end
xlabel('Ambient Pressure (kPa)');
ylabel('Amplitude (dB)');
title('Subharmonic Response in Linear Imaging');
subplot(2,2,2);
X = povs / 1e3;
Y = sh_pi;
plt = plot(X,Y); grid on;
legend(caselab,'Location','southwest');
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
end
xlabel('Ambient Pressure (kPa)');
ylabel('Amplitude (dB)');
title('Subharmonic Response in PI Imaging');
subplot(2,2,3)
X = povs / 1e3;
Y = sh_am;
plt = plot(X,Y); grid on;
legend(caselab,'Location','southwest');
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
end
xlabel('Ambient Pressure (kPa)');
ylabel('Amplitude (dB)');
title('Subharmonic Response in AM Imaging');
subplot(2,2,4)
X = povs / 1e3;
Y = sh_cps;
plt = plot(X,Y); grid on;
legend(caselab,'Location','southwest');
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
end
xlabel('Ambient Pressure (kPa)');
ylabel('Amplitude (dB)');
title('Subharmonic Response in CPS Imaging');

%% 显示不同造影模式下IQ-Demodulation分析的次谐波响应（相对变化）

%% 显示不同造影模式下IQ-Demodulation分析的次谐波响应（绝对+相对变化）
fig = figure(1006); fig.Position = [200,40,1000,1080];
subplot(4,2,1);
X = povs / 1e3;
Y = sh;
plt = plot(X,Y); grid on;
legend(caselab,'Location','best');
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
end
xlabel('Ambient Pressure (kPa)');
ylabel('Amplitude (dB)');
title({'Pressure-Depedent Subharmonic Response','Linear Imaging'});
subplot(4,2,3);
X = povs / 1e3;
Y = sh_pi;
plt = plot(X,Y); grid on;
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
end
xlabel('Ambient Pressure (kPa)');
ylabel('Amplitude (dB)');
title('PI Imaging');
subplot(4,2,5)
X = povs / 1e3;
Y = sh_am;
plt = plot(X,Y); grid on;
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
end
xlabel('Ambient Pressure (kPa)');
ylabel('Amplitude (dB)');
title('AM Imaging');
subplot(4,2,7)
X = povs / 1e3;
Y = sh_cps;
plt = plot(X,Y); grid on;
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
end
xlabel('Ambient Pressure (kPa)');
ylabel('Amplitude (dB)');
title('CPS Imaging');

subplot(4,2,2);
X = povs / 1e3;
Y = sh - ones(size(sh,1),1)*sh(1,:);
plt = plot(X,Y); grid on;
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
end
xlabel('Ambient Pressure (kPa)');
ylabel('Amplitude (dB)');
title({'Pressure-Dependent Subharmonic Variation', 'Linear Imaging'});
subplot(4,2,4);
X = povs / 1e3;
Y = sh_pi - ones(size(sh_pi,1),1)*sh_pi(1,:);
plt = plot(X,Y); grid on;
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
end
xlabel('Ambient Pressure (kPa)');
ylabel('Amplitude (dB)');
title('PI Imaging');
subplot(4,2,6)
X = povs / 1e3;
Y = sh_am - ones(size(sh_am,1),1)*sh_am(1,:);
plt = plot(X,Y); grid on;
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
end
xlabel('Ambient Pressure (kPa)');
ylabel('Amplitude (dB)');
title('AM Imaging');
subplot(4,2,8)
X = povs / 1e3;
Y = sh_cps - ones(size(sh_cps,1),1)*sh_cps(1,:);
plt = plot(X,Y); grid on;
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
end
xlabel('Ambient Pressure (kPa)');
ylabel('Amplitude (dB)');
title('CPS Imaging');

