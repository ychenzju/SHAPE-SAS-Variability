close all;
clearvars;

%% 显示参数
colors = {'b','r','k','g','m','c'};
styles = {'-','-.','--',':','-','-.'};
mrkers = {'s','^','v','o','*','d'};
widths = {0.8,1.2,1.4,1.0};
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

filename = 'plot_population_driving_frequencies.mat';
reload = 1;
if ~exist(filename,'file') || reload==1
    load population_003_Povs_Frqs
    
    xsh = zeros(length(povs),length(frqs),length(pacs));
    xfd = xsh; xuh = xsh; xhm = xsh;
    xsh_pi = xsh; xfd_pi = xsh; xuh_pi = xsh; xhm_pi = xsh;
    xsh_am = xsh; xfd_am = xsh; xuh_am = xsh; xhm_am = xsh;
    xsh_cps = xsh; xfd_cps = xsh; xuh_cps = xsh; xhm_cps = xsh;
    
    for pov_ndx = 1:length(povs)
        Pov = povs(pov_ndx);
        for frq_ndx = 1:length(frqs)
            Frq = frqs(frq_ndx);
            for pac_ndx = 1:length(pacs)
                Pac = pacs(pac_ndx);
                [Pov/1e3, Frq/1e6, Pac/1e3]
                
                %             Tt   = Treceived{pov_ndx,frq_ndx,pac_ndx};
                %             Pt   = Preceived{pov_ndx,frq_ndx,pac_ndx};
                %             sigl_list = Tt>=sg_window_beg-eps & Tt<sg_window_end+eps;
                %             [pr_roi_fft,pr_roi_frq] = freq_spectrum(Pt(sigl_list),probe_fs,'abs');
                %             SH_fc = frqs(frq_ndx) * 0.5;
                %             FD_fc = frqs(frq_ndx) * 1.0;
                %             UH_fc = frqs(frq_ndx) * 1.5;
                %             HM_fc = frqs(frq_ndx) * 2.0;
                %             fr_wd = 0.2e6; %Hz
                %
                %             xsh(pov_ndx,frq_ndx,pac_ndx) = mean(20*log10(pr_roi_fft(pr_roi_frq<SH_fc+fr_wd & pr_roi_frq>=SH_fc-fr_wd)));
                %             xfd(pov_ndx,frq_ndx,pac_ndx) = mean(20*log10(pr_roi_fft(pr_roi_frq<FD_fc+fr_wd & pr_roi_frq>=FD_fc-fr_wd)));
                %             xuh(pov_ndx,frq_ndx,pac_ndx) = mean(20*log10(pr_roi_fft(pr_roi_frq<UH_fc+fr_wd & pr_roi_frq>=UH_fc-fr_wd)));
                %             xhm(pov_ndx,frq_ndx,pac_ndx) = mean(20*log10(pr_roi_fft(pr_roi_frq<HM_fc+fr_wd & pr_roi_frq>=HM_fc-fr_wd)));
                
                Tt = Treceived{pov_ndx,frq_ndx,pac_ndx};
                Pt = Preceived{pov_ndx,frq_ndx,pac_ndx};
                SHA = rf2iq_filter(Pt,probe_fs,Frq/2);
                FDA = rf2iq_filter(Pt,probe_fs,Frq);
                UHA = rf2iq_filter(Pt,probe_fs,Frq*3/2);
                HMA = rf2iq_filter(Pt,probe_fs,Frq*2);
                ROI_LIST = (Tt>=sg_window_beg) & (Tt<=sg_window_end);
                xsh(pov_ndx,frq_ndx,pac_ndx) = mean(20*log10(SHA(ROI_LIST)));
                xfd(pov_ndx,frq_ndx,pac_ndx) = mean(20*log10(FDA(ROI_LIST)));
                xuh(pov_ndx,frq_ndx,pac_ndx) = mean(20*log10(UHA(ROI_LIST)));
                xhm(pov_ndx,frq_ndx,pac_ndx) = mean(20*log10(HMA(ROI_LIST)));
                
                Tt = Treceived{pov_ndx,frq_ndx,pac_ndx};
                Pt = PreceivedPI{pov_ndx,frq_ndx,pac_ndx};
                SHA = rf2iq_filter(Pt,probe_fs,Frq/2);
                FDA = rf2iq_filter(Pt,probe_fs,Frq);
                UHA = rf2iq_filter(Pt,probe_fs,Frq*3/2);
                HMA = rf2iq_filter(Pt,probe_fs,Frq*2);
                ROI_LIST = (Tt>=sg_window_beg) & (Tt<=sg_window_end);
                xsh_pi(pov_ndx,frq_ndx,pac_ndx) = mean(20*log10(SHA(ROI_LIST)));
                xfd_pi(pov_ndx,frq_ndx,pac_ndx) = mean(20*log10(FDA(ROI_LIST)));
                xuh_pi(pov_ndx,frq_ndx,pac_ndx) = mean(20*log10(UHA(ROI_LIST)));
                xhm_pi(pov_ndx,frq_ndx,pac_ndx) = mean(20*log10(HMA(ROI_LIST)));
                
                Tt = Treceived{pov_ndx,frq_ndx,pac_ndx};
                Pt = PreceivedAM{pov_ndx,frq_ndx,pac_ndx};
                SHA = rf2iq_filter(Pt,probe_fs,Frq/2);
                FDA = rf2iq_filter(Pt,probe_fs,Frq);
                UHA = rf2iq_filter(Pt,probe_fs,Frq*3/2);
                HMA = rf2iq_filter(Pt,probe_fs,Frq*2);
                ROI_LIST = (Tt>=sg_window_beg) & (Tt<=sg_window_end);
                xsh_am(pov_ndx,frq_ndx,pac_ndx) = mean(20*log10(SHA(ROI_LIST)));
                xfd_am(pov_ndx,frq_ndx,pac_ndx) = mean(20*log10(FDA(ROI_LIST)));
                xuh_am(pov_ndx,frq_ndx,pac_ndx) = mean(20*log10(UHA(ROI_LIST)));
                xhm_am(pov_ndx,frq_ndx,pac_ndx) = mean(20*log10(HMA(ROI_LIST)));
                
                Tt = Treceived{pov_ndx,frq_ndx,pac_ndx};
                Pt = PreceivedCPS{pov_ndx,frq_ndx,pac_ndx};
                SHA = rf2iq_filter(Pt,probe_fs,Frq/2);
                FDA = rf2iq_filter(Pt,probe_fs,Frq);
                UHA = rf2iq_filter(Pt,probe_fs,Frq*3/2);
                HMA = rf2iq_filter(Pt,probe_fs,Frq*2);
                ROI_LIST = (Tt>=sg_window_beg) & (Tt<=sg_window_end);
                xsh_cps(pov_ndx,frq_ndx,pac_ndx) = mean(20*log10(SHA(ROI_LIST)));
                xfd_cps(pov_ndx,frq_ndx,pac_ndx) = mean(20*log10(FDA(ROI_LIST)));
                xuh_cps(pov_ndx,frq_ndx,pac_ndx) = mean(20*log10(UHA(ROI_LIST)));
                xhm_cps(pov_ndx,frq_ndx,pac_ndx) = mean(20*log10(HMA(ROI_LIST)));
                
            end
        end
        save(filename,'*');
    end
else
    load(filename);
end

rads = Rs; population_figures;

sh = reshape(xsh, length(povs), length(frqs));
sh_pi = reshape(xsh_pi, length(povs), length(frqs));
sh_am = reshape(xsh_am, length(povs), length(frqs));
sh_cps = reshape(xsh_cps, length(povs), length(frqs));

%% 显示不同造影模式下IQ-Demodulation分析的次谐波响应
if size(sh,1)==6
    povs_selectedIdx = [1,3,6];
elseif size(sh,1)==3
    povs_selectedIdx = [1,2,3];
end

% frqs_selected = [2.2,2.4]*1e6;
% [~,ia,ib] = intersect(frqs_selected,frqs);
% frqs_selectedIdx = ib;
frqs_selectedIdx = [5,6,10,11,12];

povs_labels = compose("Pov = %2.0f kPa",(povs(povs_selectedIdx)/1e3));
frqs_labels = compose("f = %2.1f MHz",(frqs(frqs_selectedIdx)/1e6));

fig = figure(4006); fig.Position = [200,40,1000,1080];
subplot(4,2,1);
X = frqs(:) / 1e6;
Y = sh(povs_selectedIdx,:)';
plt = plot(X,Y); grid on;
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k}; plt(k).MarkerSize=4;
end
xlabel('Driving Frequency (MHz)');
ylabel('Amplitude (dB)');
title({'Frequency-dependent Subharmonic Response','Linear Imaging'});
legend(povs_labels,'Location','southwest');
xlim([frqs(1)/1e6,frqs(end)/1e6]);

subplot(4,2,2);
X = povs(:) / 1e3;
Y = sh(:,frqs_selectedIdx) - sh(1,frqs_selectedIdx);
plt = plot(X,Y); grid on;
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k}; plt(k).MarkerSize=4;
end
xlabel('Ambient Pressure (kPa)');
ylabel('Amplitude (dB)');
title({'Pressure-dependent Subharmonic Variation','Linear Imaging'});
legend(frqs_labels,'Location','southwest');

subplot(4,2,3);
X = frqs(:) / 1e6;
Y = sh_pi(povs_selectedIdx,:)';
plt = plot(X,Y); grid on;
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k}; plt(k).MarkerSize=4;
end
xlabel('Driving Frequency (MHz)');
ylabel('Amplitude (dB)');
title({'PI Imaging'});
legend(povs_labels,'Location','southwest');
xlim([frqs(1)/1e6,frqs(end)/1e6]);

subplot(4,2,4);
X = povs(:) / 1e3;
Y = sh_pi(:,frqs_selectedIdx) - sh_pi(1,frqs_selectedIdx);
plt = plot(X,Y); grid on;
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k}; plt(k).MarkerSize=4;
end
xlabel('Ambient Pressure (kPa)');
ylabel('Amplitude (dB)');
title({'PI Imaging'});
legend(frqs_labels,'Location','southwest');

subplot(4,2,5);
X = frqs(:) / 1e6;
Y = sh_am(povs_selectedIdx,:)';
plt = plot(X,Y); grid on;
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k}; plt(k).MarkerSize=4;
end
xlabel('Driving Frequency (MHz)');
ylabel('Amplitude (dB)');
title({'AM Imaging'});
legend(povs_labels,'Location','southwest');
xlim([frqs(1)/1e6,frqs(end)/1e6]);

subplot(4,2,6);
X = povs(:) / 1e3;
Y = sh_am(:,frqs_selectedIdx) - sh_am(1,frqs_selectedIdx);
plt = plot(X,Y); grid on;
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k}; plt(k).MarkerSize=4;
end
xlabel('Ambient Pressure (kPa)');
ylabel('Amplitude (dB)');
title({'AM Imaging'});
legend(frqs_labels,'Location','southwest');

subplot(4,2,7);
X = frqs(:) / 1e6;
Y = sh_cps(povs_selectedIdx,:)';
plt = plot(X,Y); grid on;
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k}; plt(k).MarkerSize=4;
end
xlabel('Driving Frequency (MHz)');
ylabel('Amplitude (dB)');
title({'CPS Imaging'});
legend(povs_labels,'Location','southwest');
xlim([frqs(1)/1e6,frqs(end)/1e6]);

subplot(4,2,8);
X = povs(:) / 1e3;
Y = sh_cps(:,frqs_selectedIdx) - sh_cps(1,frqs_selectedIdx);
plt = plot(X,Y); grid on;
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k}; plt(k).MarkerSize=4;
end
xlabel('Ambient Pressure (kPa)');
ylabel('Amplitude (dB)');
title({'CPS Imaging'});
legend(frqs_labels,'Location','southwest');
