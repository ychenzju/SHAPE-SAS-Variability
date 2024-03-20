close all;
clearvars;
%% 显示参数
colors = {'b','r','k','g','m','c'};
styles = {'-','-.','--',':'};
mrkers = {'s','^','v','o','*','d','+','p','<','>','h'};
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

if ~exist('plot_population_random_bubble_sizes.mat','file')
    %% plot signals
    Seeds = [1999,2099,2199,2299,2399,2499,2599,2699,2799,2899];
    caselab = {};
    SHpsd = zeros(6,length(Seeds));
    SHdem = zeros(6,length(Seeds));
    for seed_ndx = 1:length(Seeds)
        seednum = Seeds(seed_ndx);
        fname = ['population_003_R',num2str(seednum)];
        fname = [fname,'_Povs.mat'];
        caselab{seed_ndx} = ['R-seed = ',num2str(seednum)];
        load(fname);
        
        for pov_ndx = 1 : length(povs)
            [seed_ndx, pov_ndx]
            %         if pov_ndx == 1
            %             figure(3001); fig.Position = [700,100,1000,1000];
            %             subplot(length(Rseed_cases),3,(seed_ndx-1)*3+1);
            %             X = Treceived{pov_ndx} * 1e6;
            %             Y = Preceived{pov_ndx};
            %             plot(X,Y); grid on; hold on;
            %             ylim([-150,150]);
            %             YL = ylim();
            %             rectangle('Position',[sg_window_beg*1e6,YL(1),(sg_window_end-sg_window_beg)*1e6,YL(2)-YL(1)],'LineStyle',':','EdgeColor','g','LineWidth',1.4);
            %             xlabel('Time (\mus)');
            %             ylabel({caselab{seed_ndx},'Pressure (Pa)'});
            %             title('Receiving Acoustic Signal');
            %             xlim([X(1),X(end)]);
            %
            %             ax = subplot(length(Rseed_cases),3,(seed_ndx-1)*3+2);
            %             Tt = Treceived{pov_ndx};
            %             Pt = Preceived{pov_ndx};
            %             ROI_LIST = (Tt>=sg_window_beg) & (Tt<=sg_window_end);
            %             signal = Pt(ROI_LIST);
            %             [resp,freq] = freq_spectrum(Pt(ROI_LIST),probe_fs,'abs');
            %
            %             % [pows,freq] = periodogram(signal,rectwin(length(signal)),length(signal),probe_fs,'onesided','power'); resp = (pows*2).^(1/2);
            %             % [pows,freq] = pwelch(signal,length(signal),0,length(signal),probe_fs,'onesided','power'); resp = (pows*2).^(1/2);
            %             % relation between power spectrum density and power spectrum: psd = pow / (fs * N);
            %             X = freq(2:end) / 1e6; Y = 20*log10(resp(2:end));
            %             plot(X,Y,'-o','MarkerSize',3); grid on; hold on;
            %             ylim([-40,40]);
            %             YL = ylim();
            %             %rectangle('Position',[frqs(frq_ndx)/2/1e6-0.2,YL(1),0.4,YL(2)-YL(1)],'LineStyle',':','EdgeColor','r','LineWidth',1.4);
            %             xlabel('Frequency (MHz)');
            %             ylabel('Amplitude(dBa)');
            %             title('PSD-based Amplitude Spectrum');
            %             xlim([0,7.5]);
            %             ax.XTick = 0:frqs(frq_ndx)/2/1e6:7.5;
            %
            %             subplot(length(Rseed_cases),3,(seed_ndx-1)*3+3);
            %             SHA = rf2iq_filter(Pt,probe_fs,Frq/2);
            %             ROI_LIST = (Tt>=sg_window_beg) & (Tt<=sg_window_end);
            %             sh_meanvalue = mean(SHA(ROI_LIST));
            %             X = Tt * 1e6; Y = 20*log10(SHA);
            %             plot(X,Y); grid on; hold on;
            %             ylim([-40,40]);
            %             YL = ylim();
            %             rectangle('Position',[sg_window_beg*1e6,YL(1),(sg_window_end-sg_window_beg)*1e6,YL(2)-YL(1)],'LineStyle',':','EdgeColor','g','LineWidth',1.4);
            %             %line([sg_window_beg;sg_window_end]*1e6,[sh_meanvalue;sh_meanvalue],'Color','r','LineStyle',':','LineWidth',1.4);
            %             xlabel('Time (\mus)');
            %             ylabel('Amplitude (dBa)');
            %             title('IQ-demodulated Subharmonic Signal');
            %             xlim([X(1),X(end)]);
            %         end
            
            %% calculate the subharmonic amplitude vs. ambient pressure
            Tt = Treceived{pov_ndx};
            Pt = Preceived{pov_ndx};
            
            ROI_LIST = (Tt>=sg_window_beg) & (Tt<=sg_window_end);
            signal = Pt(ROI_LIST);
            [resp,freq] = freq_spectrum(signal,probe_fs,'dBW');
            sh_fc = frqs(frq_ndx)/2;
            sh_bw = 0.2e6;
            FOCUS_LIST = freq>=(sh_fc-sh_bw) & freq<=(sh_fc+sh_bw);
            SHpsd(pov_ndx,seed_ndx) = mean(resp(FOCUS_LIST));
            
            SHA = rf2iq_filter(Pt,probe_fs,Frq/2);
            ROI_LIST = (Tt>=sg_window_beg) & (Tt<=sg_window_end);
            SHdem(pov_ndx,seed_ndx) = mean(20*log10(SHA(ROI_LIST)));
            
            Tt = Treceived{pov_ndx};
            Pt = Preceived{pov_ndx};
            SHA = rf2iq_filter(Pt,probe_fs,Frq/2);
            FDA = rf2iq_filter(Pt,probe_fs,Frq);
            UHA = rf2iq_filter(Pt,probe_fs,Frq*3/2);
            HMA = rf2iq_filter(Pt,probe_fs,Frq*2);
            ROI_LIST = (Tt>=sg_window_beg) & (Tt<=sg_window_end);
            xsh(pov_ndx,seed_ndx) = mean(20*log10(SHA(ROI_LIST)));
            xfd(pov_ndx,seed_ndx) = mean(20*log10(FDA(ROI_LIST)));
            xuh(pov_ndx,seed_ndx) = mean(20*log10(UHA(ROI_LIST)));
            xhm(pov_ndx,seed_ndx) = mean(20*log10(HMA(ROI_LIST)));
            
            Tt = Treceived{pov_ndx};
            Pt = PreceivedPI{pov_ndx};
            SHA = rf2iq_filter(Pt,probe_fs,Frq/2);
            FDA = rf2iq_filter(Pt,probe_fs,Frq);
            UHA = rf2iq_filter(Pt,probe_fs,Frq*3/2);
            HMA = rf2iq_filter(Pt,probe_fs,Frq*2);
            ROI_LIST = (Tt>=sg_window_beg) & (Tt<=sg_window_end);
            xsh_pi(pov_ndx,seed_ndx) = mean(20*log10(SHA(ROI_LIST)));
            xfd_pi(pov_ndx,seed_ndx) = mean(20*log10(FDA(ROI_LIST)));
            xuh_pi(pov_ndx,seed_ndx) = mean(20*log10(UHA(ROI_LIST)));
            xhm_pi(pov_ndx,seed_ndx) = mean(20*log10(HMA(ROI_LIST)));
            
            Tt = Treceived{pov_ndx};
            Pt = PreceivedAM{pov_ndx};
            SHA = rf2iq_filter(Pt,probe_fs,Frq/2);
            FDA = rf2iq_filter(Pt,probe_fs,Frq);
            UHA = rf2iq_filter(Pt,probe_fs,Frq*3/2);
            HMA = rf2iq_filter(Pt,probe_fs,Frq*2);
            ROI_LIST = (Tt>=sg_window_beg) & (Tt<=sg_window_end);
            xsh_am(pov_ndx,seed_ndx) = mean(20*log10(SHA(ROI_LIST)));
            xfd_am(pov_ndx,seed_ndx) = mean(20*log10(FDA(ROI_LIST)));
            xuh_am(pov_ndx,seed_ndx) = mean(20*log10(UHA(ROI_LIST)));
            xhm_am(pov_ndx,seed_ndx) = mean(20*log10(HMA(ROI_LIST)));
            
            Tt = Treceived{pov_ndx};
            Pt = PreceivedCPS{pov_ndx};
            SHA = rf2iq_filter(Pt,probe_fs,Frq/2);
            FDA = rf2iq_filter(Pt,probe_fs,Frq);
            UHA = rf2iq_filter(Pt,probe_fs,Frq*3/2);
            HMA = rf2iq_filter(Pt,probe_fs,Frq*2);
            ROI_LIST = (Tt>=sg_window_beg) & (Tt<=sg_window_end);
            xsh_cps(pov_ndx,seed_ndx) = mean(20*log10(SHA(ROI_LIST)));
            xfd_cps(pov_ndx,seed_ndx) = mean(20*log10(FDA(ROI_LIST)));
            xuh_cps(pov_ndx,seed_ndx) = mean(20*log10(UHA(ROI_LIST)));
            xhm_cps(pov_ndx,seed_ndx) = mean(20*log10(HMA(ROI_LIST)));
            
        end
    end
    save('plot_population_random_bubble_sizes.mat','*');
else
    load('plot_population_random_bubble_sizes.mat');
end

sh = xsh; 
sh_pi = xsh_pi;
sh_am = xsh_am;
sh_cps = xsh_cps;

%% 显示PSD和IQ-Demodulation方法计算的次谐波-环境压力的关系
fig = figure(3002); fig.Position = [200,40,1000,1080];
subplot(2,2,1);
X = povs / 1e3;
Y = SHpsd;
plt = plot(X,Y); grid on;
legend(caselab,'Location','southwest');
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
end
xlabel('Ambient Pressure (kPa)');
ylabel('Amplitude (dBa)');
title({'Subharmonic Response','PSD Method'});
subplot(2,2,2);
X = povs / 1e3;
Y = SHpsd-ones(size(SHpsd,1),1)*SHpsd(1,:);
plt = plot(X,Y); grid on;
legend(caselab,'Location','southwest');
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
end
xlabel('Ambient Pressure (kPa)');
ylabel('Amplitude (dBa)');
title({'Relative Subharmonic Response','PSD Method'});

subplot(2,2,3);
X = povs / 1e3;
Y = SHdem;
plt = plot(X,Y); grid on;
legend(caselab,'Location','southwest');
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
end
xlabel('Ambient Pressure (kPa)');
ylabel('Amplitude (dBa)');
title({'Subharmonic Response','IQ Demodulation Method'});
subplot(2,2,4);
X = povs / 1e3;
Y = SHdem-ones(size(SHdem,1),1)*SHdem(1,:);
plt = plot(X,Y); grid on;
legend(caselab,'Location','southwest');
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
end
xlabel('Ambient Pressure (kPa)');
ylabel('Amplitude (dBa)');
title({'Relative Subharmonic Response','IQ Demodulation Method'});


%% 显示四种成像模式下的次谐波-环境压力响应
fig=figure(3004); fig.Position = [100,40,1000,1080];

subplot(4,2,1);
X = povs / 1e3;
Y = sh;
plt = plot(X,Y); grid on;
legend(caselab,'Location','southwest');
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
end
xlabel('Ambient Pressure (kPa)');
ylabel('Amplitude (dBa)');
radii = unique(R00);
if length(radii)>1
    title({'Polydisperse Microbubbles with Random Radii','Subharmonic Response in Linear Imaging'});
else
    title({['Monodisperse Microbubbles with Radii = ',num2str(radii*1e6),' \mum'],'Subharmonic Response in Linear Imaging'});
end
subplot(4,2,2);
X = povs / 1e3;
Y = sh-ones(size(sh,1),1)*sh(1,:);
plt = plot(X,Y); grid on;
legend(caselab,'Location','southwest');
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
end
xlabel('Ambient Pressure (kPa)');
ylabel('Amplitude (dBa)');
title({'Subharmonic Response (Relative) in Linear Imaging'});

subplot(4,2,3);
X = povs / 1e3;
Y = sh_pi;
plt = plot(X,Y); grid on;
legend(caselab,'Location','southwest');
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
end
xlabel('Ambient Pressure (kPa)');
ylabel('Amplitude (dBa)');
title({'Subharmonic Response in PI Imaging'});
subplot(4,2,4);
X = povs / 1e3;
Y = sh_pi-ones(size(sh_pi,1),1)*sh_pi(1,:);
plt = plot(X,Y); grid on;
legend(caselab,'Location','southwest');
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
end
xlabel('Ambient Pressure (kPa)');
ylabel('Amplitude (dBa)');
title({'Subharmonic Response (Relative) in PI Imaging'});

subplot(4,2,5);
X = povs / 1e3;
Y = sh_am;
plt = plot(X,Y); grid on;
legend(caselab,'Location','southwest');
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
end
xlabel('Ambient Pressure (kPa)');
ylabel('Amplitude (dBa)');
title({'Subharmonic Response in AM Imaging'});
subplot(4,2,6);
X = povs / 1e3;
Y = sh_am-ones(size(sh_am,1),1)*sh_am(1,:);
plt = plot(X,Y); grid on;
legend(caselab,'Location','southwest');
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
end
xlabel('Ambient Pressure (kPa)');
ylabel('Amplitude (dBa)');
title({'Subharmonic Response (Relative) in AM Imaging'});

subplot(4,2,7);
X = povs / 1e3;
Y = sh_cps;
plt = plot(X,Y); grid on;
legend(caselab,'Location','southwest');
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
end
xlabel('Ambient Pressure (kPa)');
ylabel('Amplitude (dBa)');
title({'Subharmonic Response in CPS Imaging'});
subplot(4,2,8);
X = povs / 1e3;
Y = sh_cps-ones(size(sh_cps,1),1)*sh_cps(1,:);
plt = plot(X,Y); grid on;
legend(caselab,'Location','southwest');
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
end
xlabel('Ambient Pressure (kPa)');
ylabel('Amplitude (dBa)');
title({'Subharmonic Response (Relative) in CPS Imaging'});

%% 显示次谐波信号的变动性
fig = figure(3006); fig.Position = [200,40,1000,1080];
subplot(4,1,1);
X = (1:length(Seeds))';
Y = sh';
plt = plot(X,Y);
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k}; plt(k).MarkerSize=4;
end
xlabel('Case Number (#)');
ylabel('Amplitude (dBa)');
legend(compose("Pov = %2.0f kPa", povs/1e3),'Location','best');
title('Variability of Subharmonic Amplitudes in Linear Imaging')

subplot(4,1,2);
X = (1:length(Seeds))';
Y = sh_pi';
plt = plot(X,Y);
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k}; plt(k).MarkerSize=4;
end
xlabel('Case Number (#)');
ylabel('Amplitude (dBa)');
legend(compose("Pov = %2.0f kPa", povs/1e3),'Location','best');
title('Variability of Subharmonic Amplitudes in PI Imaging')

subplot(4,1,3);
X = (1:length(Seeds))';
Y = sh_am';
plt = plot(X,Y);
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k}; plt(k).MarkerSize=4;
end
xlabel('Case Number (#)');
ylabel('Amplitude (dBa)');
legend(compose("Pov = %2.0f kPa", povs/1e3),'Location','best');
title('Variability of Subharmonic Amplitudes in AM Imaging')

subplot(4,1,4);
X = (1:length(Seeds))';
Y = sh_cps';
plt = plot(X,Y);
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k}; plt(k).MarkerSize=4;
end
xlabel('Case Number (#)');
ylabel('Amplitude (dBa)');
legend(compose("Pov = %2.0f kPa", povs/1e3),'Location','best');
title('Variability of Subharmonic Amplitudes in CPS Imaging')

