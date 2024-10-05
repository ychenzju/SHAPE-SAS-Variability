close all;
clearvars;


%% plot Monte-Carlo simulation results
%% The user must choose the population firstly.
%% one is for polydisperse populations
%% the other is for monodisper populations of 3um
%% This can be done by commenting and uncommenting the following codes

% MCcaseLists = dir('population_003_MC_*'); %for polydisperse populations
% filename = 'plot_population_random_bubbles.mat'; %for polydisper populations
MCcaseLists = dir('population_003_R3MC_*'); %for monodisperse populations
filename = 'plot_population_random_bubbles_R3.mat'; %for monodisperse populations

reload = 0;
if ~exist(filename,'file') || reload==1
    caselab = {};
    existingLargeBubbles = zeros(length(MCcaseLists),1);
    
    SHpsd = zeros(6,length(MCcaseLists));
    SHdem = zeros(6,length(MCcaseLists));
    SHnrm = zeros(6,length(MCcaseLists));
    for seed_ndx = 1:length(MCcaseLists)
        
        fname = MCcaseLists(seed_ndx).name;
        disp(fname);
        pattern1 = "R" + digitsPattern + "_L";
        pattern2 = "_L" + digitsPattern + "_Povs";
        pattern3 = "_Povs.mat";
        pos1 = strfind(fname,pattern1);
        pos2 = strfind(fname,pattern2);
        pos3 = strfind(fname,pattern3);
        
        Rseed = str2num(fname(pos1+1:pos2-1));
        Lseed = str2num(fname(pos2+2:pos3-1));
        
        caselab{seed_ndx} = ['R-seed = ',num2str(Rseed),', L-seed = ',num2str(Lseed)];
        load(fname);
        
        roi_list = Zs>=-4.6e-2;
        existingLargeBubbles(seed_ndx) = sum(Rs(roi_list)>=5e-6);
        
        Pac = pacs(1);
        Frq = frqs(1);
        for pov_ndx = 1 : length(povs)
            Pov = povs(pov_ndx);
            %         if pov_ndx == 6
            %             fig = figure(3001); fig.Position = [700,100,1000,1000];
            %             subplot(length(caseLists),3,(case_ndx-1)*3+1);
            %             X = Treceived{pov_ndx} * 1e6;
            %             Y = Preceived{pov_ndx};
            %             plot(X,Y); grid on; hold on;
            %             ylim([-150,150]);
            %             YL = ylim();
            %             rectangle('Position',[sg_window_beg*1e6,YL(1),(sg_window_end-sg_window_beg)*1e6,YL(2)-YL(1)],'LineStyle',':','EdgeColor','g','LineWidth',1.4);
            %             xlabel('Time (\mus)');
            %             ylabel({caselab{case_ndx},'Pressure (Pa)'});
            %             title('Receiving Acoustic Signal');
            %             xlim([X(1),X(end)]);
            %
            %             ax = subplot(length(caseLists),3,(case_ndx-1)*3+2);
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
            %             subplot(length(caseLists),3,(case_ndx-1)*3+3);
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
    save(filename,'*');
else
    load(filename);
end

%% 显示参数
colors = {'b','r','g','k','m','c'};
styles = {':',':',':',':'};
mrkers = {'s','^','v','o','*','d','+','p','<','>','h'};
widths = {0.8,0.8,0.8,0.8};
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

%% 显示PSD方法的次谐波响应
fig = figure(3002); fig.Position = [500,100,1000,1000];
subplot(2,2,1);
X = povs / 1e3;
Y = SHpsd;
plt = plot(X,Y); grid on;
% legend(caselab,'Location','southwest');
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
% legend(caselab,'Location','southwest');
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
end
xlabel('Ambient Pressure (kPa)');
ylabel('Amplitude (dBa)');
title({'Relative Subharmonic Response','PSD Method'});
subplot(2,2,3);
X = (1:size(SHpsd,2))';
Y = SHpsd';
plt = plot(X,Y,'o','MarkerSize',3); grid on;
legend(compose("Pov = %2.0f kPa",povs/1e3),'Location','southwest');
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).Marker=mrkers{k}; plt(k).MarkerFaceColor=colors{k};
end
xlabel('Monte-Carlo Test Number (#)');
ylabel('Amplitude (dBa)');
title({'Subharmonic Response','PSD Method'});
subplot(2,2,4);
X = (1:size(SHpsd,2))';
Y = SHpsd(end,:)'-SHpsd(1,:)';
plt = plot(X,Y,'-o','MarkerSize',3); grid on;
% legend(compose("Pov = %2.0f kPa",povs/1e3),'Location','southwest');
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).Marker=mrkers{k}; plt(k).MarkerFaceColor=colors{k};
end
XL = xlim();
line(XL,[mean(Y),mean(Y)],'Color','r','LineStyle',':','LineWidth',1.0);
xlabel('Monte-Carlo Test Number (#)');
ylabel('Amplitude (dBa)');
title({'Subharmonic Sensitivity','PSD Method'});

%% 显示IQ-Demodulation方法的次谐波响应
fig = figure(3003); fig.Position = [100,100,1000,1000];
subplot(2,2,1);
X = povs / 1e3;
Y = SHdem;
plt = plot(X,Y); grid on;
% legend(caselab,'Location','southwest');
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
end
xlabel('Ambient Pressure (kPa)');
ylabel('Amplitude (dBa)');
title({'Subharmonic Response','IQ Demodulation Method'});
subplot(2,2,2);
X = povs / 1e3;
Y = SHdem-ones(size(SHdem,1),1)*SHdem(1,:);
plt = plot(X,Y); grid on;
% legend(caselab,'Location','southwest');
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
end
xlabel('Ambient Pressure (kPa)');
ylabel('Amplitude (dBa)');
title({'Relative Subharmonic Response','IQ Demodulation Method'});
subplot(2,2,3);
X = (1:size(SHdem,2))';
Y = SHdem';
plt = plot(X,Y,'o','MarkerSize',3); grid on;
legend(compose("Pov = %2.0f kPa",povs/1e3),'Location','southwest');
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).Marker=mrkers{k}; plt(k).MarkerFaceColor=colors{k};
end
xlabel('Monte-Carlo Test Number (#)');
ylabel('Amplitude (dBa)');
title({'Subharmonic Response','IQ Demodulation Method'});
subplot(2,2,4);
X = (1:size(SHdem,2))';
Y = SHdem(end,:)'-SHdem(1,:)';
plt = plot(X,Y,'-o','MarkerSize',3); grid on;
% legend(compose("Pov = %2.0f kPa",povs/1e3),'Location','southwest');
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).Marker=mrkers{k}; plt(k).MarkerFaceColor=colors{k};
end
XL = xlim();
line(XL,[mean(Y),mean(Y)],'Color','r','LineStyle',':','LineWidth',1.0);
xlabel('Monte-Carlo Test Number (#)');
ylabel('Amplitude (dBa)');
title({'Subharmonic Sensitivity','IQ Demodulation Method'});

%% 显示四种成像模式下的次谐波-环境压力响应
sh = xsh; 
sh_pi = xsh_pi;
sh_am = xsh_am;
sh_cps = xsh_cps;

fig=figure(3004); fig.Position = [100,40,1000,1080];

subplot(4,2,1);
X = povs / 1e3;
Y = sh;
plt = plot(X,Y); grid on;
legend(caselab,'Location','southwest');
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k}; plt(k).MarkerSize=3;
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
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k}; plt(k).MarkerSize=3;
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
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k}; plt(k).MarkerSize=3;
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
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k}; plt(k).MarkerSize=3;
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
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k}; plt(k).MarkerSize=3;
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
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k}; plt(k).MarkerSize=3;
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
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k}; plt(k).MarkerSize=3;
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
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k}; plt(k).MarkerSize=3;
end
xlabel('Ambient Pressure (kPa)');
ylabel('Amplitude (dBa)');
title({'Subharmonic Response (Relative) in CPS Imaging'});

%% 显示次谐波信号的变动性
pov_list = [1,6];
fig = figure(3006); fig.Position = [200,40,1000,1080];
subplot(4,2,1); grid on;
X = (1:length(MCcaseLists))';
Y = sh(pov_list,:)';
plt = plot(X,Y); grid on; hold on;
%line([1,length(caseLists)],[mean(Y),mean(Y)],'Color','r','LineStyle',':','LineWidth',1.0);
XL = xlim(); YL = ylim();
legend(compose("Pov = %2.0f kPa, %2.1f \x00B1 %2.1f",(povs(pov_list)/1e3)',(mean(Y))',(std(Y))'),'Location','best');
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k}; plt(k).MarkerSize=3;
end
xlabel('Case Number (#)');
ylabel('Amplitude (dB)');
%title({'Monte-Carlo Simulated Subharmonic Amplitudes','Linear Imaging'});
title({'Monte-Carlo Simulated Subharmonic Amplitudes',''});

subplot(4,2,2);
X = (1:length(MCcaseLists))';
Y = (sh(end,:)' - sh(1,:)')/187;
plt = plot(X,Y); grid on;hold on;
line([1,length(MCcaseLists)],[mean(Y),mean(Y)],'Color','r','LineStyle','--','LineWidth',1.0);
for k=1:length(plt)
    plt(k).Color='k'; plt(k).LineStyle=styles{k}; plt(k).Marker='*'; plt(k).LineWidth=widths{k}; plt(k).MarkerSize=3;
end
% plot(X(existingLargeBubbles>eps),Y(existingLargeBubbles>eps),'ro');
legend(sprintf('\x0394ShA / \x0394Pa = %2.3f \x00B1 %2.3f',mean(Y),std(Y)),'Location','southeast');
xlabel('Case Number (#)');
ylabel('Sensitivity (dB / mmHg)');
%title({'\DeltaShA / \DeltaPa over 0 - 25 kPa Ambient Pressure','Linear Imaging'});
title({'\DeltaShA / \DeltaPa over 0 - 25 kPa Ambient Pressure',''});

subplot(4,2,3);
X = (1:length(MCcaseLists))';
Y = sh_pi(pov_list,:)';
plt = plot(X,Y); grid on;hold on;
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k}; plt(k).MarkerSize=3;
end
%line([1,length(caseLists)],[mean(Y),mean(Y)],'Color','r','LineStyle',':','LineWidth',1.0);
XL = xlim(); YL = ylim();
legend(compose("Pov = %2.0f kPa, %2.1f \x00B1 %2.1f",(povs(pov_list)/1e3)',(mean(Y))',(std(Y))'),'Location','best');
xlabel('Case Number (#)');
ylabel('Amplitude (dB)');
%title({'PI Imaging'});

subplot(4,2,4); grid on;
X = (1:length(MCcaseLists))';
Y = (sh_pi(end,:)' - sh_pi(1,:)')/187;
plt = plot(X,Y); grid on;hold on;
line([1,length(MCcaseLists)],[mean(Y),mean(Y)],'Color','r','LineStyle','--','LineWidth',1.0);
for k=1:length(plt)
    plt(k).Color='k'; plt(k).LineStyle=styles{k}; plt(k).Marker='*'; plt(k).LineWidth=widths{k}; plt(k).MarkerSize=3;
end
% plot(X(existingLargeBubbles>eps),Y(existingLargeBubbles>eps),'ro');
legend(sprintf('\x0394ShA / \x0394Pa = %2.3f \x00B1 %2.3f',mean(Y),std(Y)),'Location','southeast');
xlabel('Case Number (#)');
ylabel('Sensitivity (dB / mmHg)');
%title({'PI Imaging'});

subplot(4,2,5);
X = (1:length(MCcaseLists))';
Y = sh_am(pov_list,:)';
plt = plot(X,Y); grid on;hold on;
%line([1,length(caseLists)],[mean(Y),mean(Y)],'Color','r','LineStyle',':','LineWidth',1.0);
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k}; plt(k).MarkerSize=3;
end
XL = xlim(); YL = ylim();
legend(compose("Pov = %2.0f kPa, %2.1f \x00B1 %2.1f",(povs(pov_list)/1e3)',(mean(Y))',(std(Y))'),'Location','best');
xlabel('Case Number (#)');
ylabel('Amplitude (dB)');
%title({'AM Imaging'});

subplot(4,2,6);
X = (1:length(MCcaseLists))';
Y = (sh_am(end,:)' - sh_am(1,:)')/187;
plt = plot(X,Y); grid on;hold on;
line([1,length(MCcaseLists)],[mean(Y),mean(Y)],'Color','r','LineStyle','--','LineWidth',1.0);
for k=1:length(plt)
    plt(k).Color='k'; plt(k).LineStyle=styles{k}; plt(k).Marker='*'; plt(k).LineWidth=widths{k}; plt(k).MarkerSize=3;
end
% plot(X(existingLargeBubbles>eps),Y(existingLargeBubbles>eps),'ro');
legend(sprintf('\x0394ShA / \x0394Pa = %2.3f \x00B1 %2.3f',mean(Y),std(Y)),'Location','southeast');
xlabel('Case Number (#)');
ylabel('Sensitivity (dB / mmHg)');
%title({'AM Imaging'});

subplot(4,2,7);
X = (1:length(MCcaseLists))';
Y = sh_cps(pov_list,:)';
plt = plot(X,Y); grid on;hold on;
%line([1,length(caseLists)],[mean(Y),mean(Y)],'Color','r','LineStyle',':','LineWidth',1.0);
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k}; plt(k).MarkerSize=3;
end
XL = xlim(); YL = ylim();
legend(compose("Pov = %2.0f kPa, %2.1f \x00B1 %2.1f",(povs(pov_list)/1e3)',(mean(Y))',(std(Y))'),'Location','best');
xlabel('Case Number (#)');
ylabel('Amplitude (dB)');
%title({'CPS Imaging'});

subplot(4,2,8);
X = (1:length(MCcaseLists))';
Y = (sh_cps(end,:)' - sh_cps(1,:)')/187;
plt = plot(X,Y); grid on;hold on;
line([1,length(MCcaseLists)],[mean(Y),mean(Y)],'Color','r','LineStyle','--','LineWidth',1.0);
for k=1:length(plt)
    plt(k).Color='k'; plt(k).LineStyle=styles{k}; plt(k).Marker='*'; plt(k).LineWidth=widths{k}; plt(k).MarkerSize=3;
end
% plot(X(existingLargeBubbles>eps),Y(existingLargeBubbles>eps),'ro');
legend(sprintf('\x0394ShA / \x0394Pa = %2.3f \x00B1 %2.3f',mean(Y),std(Y)),'Location','southeast');
xlabel('Case Number (#)');
ylabel('Sensitivity (dB / mmHg)');
%title({'CPS Imaging'});

