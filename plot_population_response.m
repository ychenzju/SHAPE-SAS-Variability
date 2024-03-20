close all; 
clearvars; 
load population_003_Povs.mat;
close all;

% Preceived = PreceivedPI;

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

%% 显示微泡位置和随机分布
% plot the randomly distributed bubbles
figure(1);
subplot(1,2,1);
histogram(Rs*1e6,12);
subplot(1,2,2);
scatter3(Xs,Ys,Zs,Rs*1e6); xlabel('X'); ylabel('Y'); zlabel('Z'); title('Randomly generated bubbles');
hold on;

%% 显示发射和接收信号
pov_ndx = 1;
frq_ndx = 1;
r = lognrnd(lognorm_mean,lognorm_std,100000,1);
xbins = 0.25:0.5:9;
counts = hist(r,xbins);
percents = counts/sum(counts);
if percents(end)>percents(end-1)
    percents(end) = percents(end-1)/2;
    percents = percents / sum(percents);
end
 
fig=figure(2); fig.Position = [300,300,810,810];
subplot(3,2,1);
histogram(Rs*1e6,'BinEdges',(0:0.5:9)); hold on;
plot(xbins,percents*length(Rs),'-*');
xlabel('Bubble Radius (\mum)')
ylabel('Counts');
title('Bubble size distribution');
legend('Randomly generated','Ideally distributed')

subplot(3,2,2);
scatter(Xs*1e2,Zs*1e2,Rs*1e6,'b'); hold on; 
xlabel('X (cm)'); ylabel('Z (cm)'); title('1D randomly generated bubbles');
ylim([-0.05,-0.04]*1e2);
XL = xlim();
line(XL,[-0.04,-0.04]*1e2,'Color','r'); 
text(XL(1),-0.04*1e2,' Vessel Top','HorizontalAlignment','left','VerticalAlignment','top');
line(XL,[-0.043,-0.043]*1e2,'Color','r'); 
text(XL(1),-0.043*1e2,' ROI Top','HorizontalAlignment','left','VerticalAlignment','top');
line(XL,[-0.046,-0.046]*1e2,'Color','r');
text(XL(1),-0.046*1e2,' ROI Bottom','HorizontalAlignment','left','VerticalAlignment','top');
line(XL,[-0.049,-0.049]*1e2,'Color','r');
text(XL(1),-0.049*1e2,' Vessel Bottom','HorizontalAlignment','left','VerticalAlignment','top');

subplot(3,2,3);
X = tex * 1e6;
Y = xex / 1e3;
Y2 = uex * Pac / 1e3;
plot(X,Y,'LineWidth',1.0); grid on; hold on
plot(X,Y2,'r--');
xlabel('Time (\mus)');
ylabel('Pressure (kPa)');
title('Excitation signals');
legend('Transmitting acoustic signal','Driving signal')
subplot(3,2,4);
X = Treceived{pov_ndx} * 1e6;
Y = Preceived{pov_ndx};
plot(X,Y); grid on; hold on;
YL = ylim();
rectangle('Position',[sg_window_beg*1e6,YL(1),(sg_window_end-sg_window_beg)*1e6,YL(2)-YL(1)],'LineStyle',':','EdgeColor','g','LineWidth',1.4);
xlabel('Time (\mus)');
ylabel('Pressure (Pa)');
title('Receiving acoustic signal');
xlim([X(1),X(end)]);

ax=subplot(3,2,5);
Tt = Treceived{pov_ndx};
Pt = Preceived{pov_ndx};
ROI_LIST = (Tt>=sg_window_beg) & (Tt<=sg_window_end);
signal = Pt(ROI_LIST);
% [resp,freq] = freq_spectrum(Pt(ROI_LIST),probe_fs,'abs');
[pows,freq,powsc] = periodogram(signal,rectwin(length(signal)),length(signal),probe_fs,'onesided','power','ConfidenceLevel',0.95); 
resp = 10*log10(pows*2);
respc = 10*log10(powsc*2);
% [pows,freq] = pwelch(signal,length(signal),0,length(signal),probe_fs,'onesided','power'); resp = (pows*2).^(1/2);
% relation between power spectrum density and power spectrum: psd = pow / (fs * N);
X = freq(2:end) / 1e6; Y = (resp(2:end));
sh_index = find(abs(X-1.25)<0.1);
plot(X,Y,'-o','MarkerSize',4); grid on; hold on;
plot(X(sh_index),Y(sh_index),'r*');
plot(X,respc(2:end,:),'c:','LineWidth',1.2);
ylim([-10,70]);
YL = ylim();
%rectangle('Position',[frqs(frq_ndx)/2/1e6-0.2,YL(1),0.4,YL(2)-YL(1)],'LineStyle',':','EdgeColor','r','LineWidth',1.4);
xlabel('Frequency (MHz)');
ylabel('Amplitude (dB)');
title('Power Spectrum of receiving signal');
xlim([0,7.5]); 
ax.XTick = 0:frqs(frq_ndx)/2/1e6:7.5;
legend('Estimated spectrum','Subharmonic amplitude','95% Confidence bounds');

subplot(3,2,6);
SHA = rf2iq_filter(Pt,probe_fs,Frq/2);
ROI_LIST = (Tt>=sg_window_beg) & (Tt<=sg_window_end);
X = Tt * 1e6; Y = 20*log10(SHA);
sh_meanvalue = mean(Y(ROI_LIST));
plot(X,Y); grid on; hold on;
ylim([0,50]);
YL = ylim();
rectangle('Position',[sg_window_beg*1e6,YL(1),(sg_window_end-sg_window_beg)*1e6,YL(2)-YL(1)],'LineStyle',':','EdgeColor','g','LineWidth',1.4);
line([sg_window_beg;sg_window_end]*1e6,[sh_meanvalue;sh_meanvalue],'Color','r','LineStyle',':','LineWidth',1.4);
xlabel('Time (\mus)');
ylabel('Amplitude (dB)');
title('IQ-demodulated subharmonic signal');
xlim([X(1),X(end)]);
legend('Subharmonic signal','Average amplitude over ROI');

%% 显示接收信号及其频谱
for pov_ndx = 1:length(povs)
    for frq_ndx = 1:length(frqs)
        for pac_ndx = 1:length(pacs)
            
            Tt   = Treceived{pov_ndx,frq_ndx,pac_ndx};
            Pt   = Preceived{pov_ndx,frq_ndx,pac_ndx};
            
            sigl_list = Tt>=sg_window_beg-eps & Tt<sg_window_end+eps;
            [pr_roi_fft,pr_roi_frq] = freq_spectrum(Pt(sigl_list),probe_fs,'abs');
            SH_fc = frqs(frq_ndx) * 0.5;
            FD_fc = frqs(frq_ndx) * 1.0;
            UH_fc = frqs(frq_ndx) * 1.5;
            HM_fc = frqs(frq_ndx) * 2.0;
            fr_wd = 0.2e6; %Hz
            
            xsh(pov_ndx,frq_ndx,pac_ndx) = mean(20*log10(pr_roi_fft(pr_roi_frq<SH_fc+fr_wd & pr_roi_frq>=SH_fc-fr_wd)));
            xfd(pov_ndx,frq_ndx,pac_ndx) = mean(20*log10(pr_roi_fft(pr_roi_frq<FD_fc+fr_wd & pr_roi_frq>=FD_fc-fr_wd)));
            xuh(pov_ndx,frq_ndx,pac_ndx) = mean(20*log10(pr_roi_fft(pr_roi_frq<UH_fc+fr_wd & pr_roi_frq>=UH_fc-fr_wd)));
            xhm(pov_ndx,frq_ndx,pac_ndx) = mean(20*log10(pr_roi_fft(pr_roi_frq<HM_fc+fr_wd & pr_roi_frq>=HM_fc-fr_wd)));
            
            fig=figure(3); fig.Position = [200,200,1680,240];
            subplot(121), plot(recv_time*1e6,Pt); grid on; hold on; YL = ylim();
            rectangle('Position',[sg_window_beg*1e6,YL(1),(sg_window_end-sg_window_beg)*1e6,YL(2)-YL(1)],'EdgeColor','g');
            title('Received Acoustic Signal'); xlabel('Time (\mus)'); ylabel('Pressure (Pa)');
            subplot(122), plot(pr_roi_frq(2:end)/1e6,20*log10(pr_roi_fft(2:end))); grid on; hold on;
            title('Spectrum of Received Signal'); xlabel('Freq (MHz)'); ylabel('Amplitude (dB)');
            xlim([0,7.5]);
            
        end
    end
end
disp(xsh(:)');

%% 显示谐波响应vs环境压力
% rads = Rs; population_figures;

%% 显示每个微泡的振动和散射信号
ROIset = find((Zs <= -(Probe2Vessel_Depth + Vessel2ROI_Depth)) & (Zs >= -(Probe2Vessel_Depth + Vessel2ROI_Depth + ROI_Zwidth)));
TOPset = find((Zs <= -(Probe2Vessel_Depth)) & (Zs >= -(Probe2Vessel_Depth + Vessel2ROI_Depth)));
BOTset = find((Zs <= -(Probe2Vessel_Depth + Vessel2ROI_Depth + ROI_Zwidth)) & (Zs >= -(Probe2Vessel_Depth + Vessel2ROI_Depth + ROI_Zwidth + ROI2Bottom_Depth)));

RROI = R00(ROIset);
RTOP = R00(TOPset);
RBOT = R00(BOTset);

pov_roi = [0,10,25] * 1e3;
[~,pov_lst] = min(abs(pov_roi - povs(:)));

pov_ndx = 1;
for rad_ndx = 1:length(R00)
    rad = Rinit(rad_ndx,pov_ndx);
    fig=figure(200); fig.Position = [180,80,980,1024];
    subplot(1,2,1);
    X = Tincident{pov_ndx}*1e6;
    Y = (Rincident{pov_ndx}(:,rad_ndx)-rad)*1e6-rad_ndx;
    ln = plot(X,Y,'Color',colors{rad_ndx});
    if ismember(rad_ndx,ROIset)
        ln.LineWidth = 1.4;
    end
    text(X(1),Y(1),['R = ',num2str(rad*1e6,'%1.2f'),' '],'HorizontalAlignment','right','VerticalAlignment','middle');
    hold on;
    subplot(1,2,2);
    X = Treceived{pov_ndx}*1e6;
    Y = (Psimulated{pov_ndx}(:,rad_ndx))-rad_ndx*20;
    ln = plot(X,Y,'Color',colors{rad_ndx});
    if ismember(rad_ndx,ROIset)
        ln.LineWidth = 1.4;
    end
     text(X(1),Y(1),['#',num2str(rad_ndx),' '],'HorizontalAlignment','right','VerticalAlignment','middle');
   hold on;
end
fig=figure(200);
ax = subplot(121);
XL = xlim();
xlim([XL(1)-2, XL(2)]);
title('Bubble Radius Dynamics');
xlabel('Time (\mus)');
ylim([-56,1]);
ax.YTickLabel = [];
%YL = ylim();
%rectangle('Position',[ix_window_beg*1e6,YL(1),(ix_window_end-ix_window_beg)*1e6,YL(2)-YL(1)],'EdgeColor','k');
ax = subplot(122);
XL = xlim();
xlim([XL(1), XL(2)]);
ylim([-1120,20]);
YL = ylim();
rectangle('Position',[sg_window_beg*1e6,YL(1),(sg_window_end-sg_window_beg)*1e6,YL(2)-YL(1)],...
    'EdgeColor','g','LineWidth',1.4,'LineStyle',':');
title('Receiving Signals from Each Bubble')
xlabel('Time (\mus)');
ax.YTickLabel = [];

%% RF2IQ混频法次谐波信号处理
SH = zeros(length(povs),1);
FD = zeros(length(povs),1);
UH = zeros(length(povs),1);
HM = zeros(length(povs),1);
Frq = frqs;
Pac = pacs;
for pov_ndx=1:length(povs)
    Pov = povs(pov_ndx);
    
    Tt = Treceived{pov_ndx};
    Pt = Preceived{pov_ndx};
    
    SHA = rf2iq_filter(Pt,probe_fs,Frq/2);
    FDA = rf2iq_filter(Pt,probe_fs,Frq);
    UHA = rf2iq_filter(Pt,probe_fs,Frq*3/2);
    HMA = rf2iq_filter(Pt,probe_fs,Frq*2);
    
    %显示混频法的波形
%     fig = figure(400); fig.Position = [100,500,1600,200];
%     subplot(1,length(povs),pov_ndx);
%     plot(Tt*1e6,[SHA,FDA,UHA,HMA]);
%     grid on;
%     if pov_ndx==1
%         legend('SH','FD','UH','HM');
%     end
%     title(['Pov = ',num2str(Pov/1e3,'%2.0f'),' kPa']);
%     xlabel('Time (\mus)');
%     ylabel('Amplitude (AU)');
%     xlim([Tt(1),Tt(end)]*1e6);
% %     ylim([0,80]);
%     YL = ylim();
%     rectangle('Position',[sg_window_beg*1e6,YL(1),(sg_window_end-sg_window_beg)*1e6,YL(2)-YL(1)],...
%     'EdgeColor','g','LineWidth',1.4,'LineStyle',':');
    
    ROI_LIST = (Tt>=sg_window_beg) & (Tt<=sg_window_end);
    SH(pov_ndx) = mean(20*log10(SHA(ROI_LIST)));
    FD(pov_ndx) = mean(20*log10(FDA(ROI_LIST)));
    UH(pov_ndx) = mean(20*log10(UHA(ROI_LIST)));
    HM(pov_ndx) = mean(20*log10(HMA(ROI_LIST)));
end
disp(SH');
fig=figure(301); fig.Position = [100,200,1600,400]; subplot(1,3,1);
plt = plot(povs(:)/1e3,([SH,FD,UH,HM]),'LineWidth',1.2);
plt(1).Marker = 's';
plt(2).Marker = '^';
plt(3).Marker = 'v';
plt(4).Marker = 'o';
grid on;
lgd_labels = {['Subharmonic \Delta = ',num2str(SH(end)-SH(1),'%2.2f'), 'dB'],...
    ['Fundamental \Delta = ',num2str(FD(end)-FD(1),'%2.2f'), 'dB'],...
    ['Ultraharmonic \Delta = ',num2str(UH(end)-UH(1),'%2.2f'), 'dB'],...
    ['2nd-harmonic \Delta = ',num2str(HM(end)-HM(1),'%2.2f'), 'dB']};
legend(lgd_labels);
xlabel('Ambient Pressure (kPa)');
ylabel('Amplitude (dB)');
title('Harmonic Response (IQ-Demod) vs. Ambient Pressure')

%% FFT法次谐波信号处理 -- 谐波幅值
SH = zeros(length(povs),1);
FD = zeros(length(povs),1);
UH = zeros(length(povs),1);
HM = zeros(length(povs),1);
Frq = frqs;
Pac = pacs;
for pov_ndx=1:length(povs)
    Pov = povs(pov_ndx);
    
    Tt = Treceived{pov_ndx};
    Pt = Preceived{pov_ndx};
    
    ROI_LIST = (Tt>=sg_window_beg) & (Tt<=sg_window_end);
    [resp,freq] = freq_spectrum(Pt(ROI_LIST),probe_fs,'abs');
    SH_fc = frqs(frq_ndx) * 0.5;
    FD_fc = frqs(frq_ndx) * 1.0;
    UH_fc = frqs(frq_ndx) * 1.5;
    HM_fc = frqs(frq_ndx) * 2.0;
    fr_wd = 0.2e6; %Hz
    
    SH(pov_ndx) = mean(20*log10(resp(freq<SH_fc+fr_wd & freq>=SH_fc-fr_wd)));
    FD(pov_ndx) = mean(20*log10(resp(freq<FD_fc+fr_wd & freq>=FD_fc-fr_wd)));
    UH(pov_ndx) = mean(20*log10(resp(freq<UH_fc+fr_wd & freq>=UH_fc-fr_wd)));
    HM(pov_ndx) = mean(20*log10(resp(freq<HM_fc+fr_wd & freq>=HM_fc-fr_wd)));

end
disp(SH');
figure(301); subplot(1,3,2);
plt = plot(povs(:)/1e3,([SH,FD,UH,HM]),'LineWidth',1.2);
plt(1).Marker = 's';
plt(2).Marker = '^';
plt(3).Marker = 'v';
plt(4).Marker = 'o';
grid on;
lgd_labels = {['Subharmonic \Delta = ',num2str(SH(end)-SH(1),'%2.2f'), 'dB'],...
    ['Fundamental \Delta = ',num2str(FD(end)-FD(1),'%2.2f'), 'dB'],...
    ['Ultraharmonic \Delta = ',num2str(UH(end)-UH(1),'%2.2f'), 'dB'],...
    ['2nd-harmonic \Delta = ',num2str(HM(end)-HM(1),'%2.2f'), 'dB']};
legend(lgd_labels);
xlabel('Ambient Pressure (kPa)');
ylabel('Amplitude (dB)');
title('Harmonic Response (DSP) vs. Ambient Pressure')

%% FFT法次谐波信号处理 -- 归一化幅值
SH = zeros(length(povs),1);
FD = zeros(length(povs),1);
UH = zeros(length(povs),1);
HM = zeros(length(povs),1);
Frq = frqs;
Pac = pacs;
for pov_ndx=1:length(povs)
    Pov = povs(pov_ndx);
    
    Tt = Treceived{pov_ndx};
    Pt = Preceived{pov_ndx};
    
    ROI_LIST = (Tt>=sg_window_beg) & (Tt<=sg_window_end);
    [resp,freq] = freq_spectrum(Pt(ROI_LIST),probe_fs,'norm-dBW');
    SH_fc = frqs(frq_ndx) * 0.5;
    FD_fc = frqs(frq_ndx) * 1.0;
    UH_fc = frqs(frq_ndx) * 1.5;
    HM_fc = frqs(frq_ndx) * 2.0;
    fr_wd = 0.2e6; %Hz
    
    SH(pov_ndx) = (mean(resp(freq<SH_fc+fr_wd & freq>=SH_fc-fr_wd)));
    FD(pov_ndx) = (mean(resp(freq<FD_fc+fr_wd & freq>=FD_fc-fr_wd)));
    UH(pov_ndx) = (mean(resp(freq<UH_fc+fr_wd & freq>=UH_fc-fr_wd)));
    HM(pov_ndx) = (mean(resp(freq<HM_fc+fr_wd & freq>=HM_fc-fr_wd)));

end
figure(301); subplot(1,3,3);
plt = plot(povs(:)/1e3,([SH,FD,UH,HM]),'LineWidth',1.2);
plt(1).Marker = 's';
plt(2).Marker = '^';
plt(3).Marker = 'v';
plt(4).Marker = 'o';
grid on;
lgd_labels = {['Subharmonic \Delta = ',num2str(SH(end)-SH(1),'%2.2f'), 'dB'],...
    ['Fundamental \Delta = ',num2str(FD(end)-FD(1),'%2.2f'), 'dB'],...
    ['Ultraharmonic \Delta = ',num2str(UH(end)-UH(1),'%2.2f'), 'dB'],...
    ['2nd-harmonic \Delta = ',num2str(HM(end)-HM(1),'%2.2f'), 'dB']};
legend(lgd_labels);
xlabel('Ambient Pressure (kPa)');
ylabel('Amplitude (dB)');
title('Harmonic Response (Power Normalized) vs. Ambient Pressure')

return

%% 显示每个微泡的次谐波响应 -- ROI区域
for r_ndx = 1:length(ROIset)
    rad_ndx = ROIset(r_ndx);
    rad = R00(rad_ndx);
    for p_ndx = 1:length(pov_lst)
        pov_ndx = pov_lst(p_ndx);
        pov = povs(pov_ndx);
        
        T = Tincident{pov_ndx};
        R = Rincident{pov_ndx}(:,rad_ndx);
        TN = Treceived{pov_ndx};
        PN = Psimulated{pov_ndx}(:,rad_ndx);
        
        ROI_LIST = find(TN>=sg_window_beg & TN<=sg_window_end);
        [resp,freq] = freq_spectrum(PN(ROI_LIST)/Pac,probe_fs,'abs');
        
        fig=figure(100); fig.Position=[680 80 980 1024];
        subplot(length(ROIset),3,(r_ndx-1)*3+1),
        plot(T*1e6,R*1e6,'Color',colors{p_ndx});
        grid on; hold on;
        xlim([ix_window_beg,ix_window_end+6e-6]*1e6);
        ylim([rad*0.5,rad+2e-6]*1e6);
        YL = ylim();
        rectangle('Position',[ix_window_beg*1e6,YL(1),(ix_window_end-ix_window_beg)*1e6,YL(2)-YL(1)],'EdgeColor','g');
%         xlabel('Time (\mus)');
%         ylabel('Radius (\mum)');
        
        subplot(length(ROIset),3,(r_ndx-1)*3+2),
        plot(TN*1e6,PN,'Color',colors{p_ndx});
        grid on; hold on;
        ylim([-60,60])
        YL = ylim();
        rectangle('Position',[sg_window_beg*1e6,YL(1),(sg_window_end-sg_window_beg)*1e6,YL(2)-YL(1)],'EdgeColor','g');
%         xlabel('Time (\mus)');
%         ylabel('Pressure (Pa)');
        
        ax = subplot(length(ROIset),3,(r_ndx-1)*3+3);
        plot(freq(2:end)/1e6,20*log10(resp(2:end)),'Color',colors{p_ndx});
        grid on; hold on;
        xlim([0,7]);
        ax.XTick = (0:Frq/1e6/2:7);
%         xlabel('Freq (MHz)');
%         ylabel('Amplitude (dB)');
        
    end
end

%% 显示每个微泡的次谐波响应 -- TOP区域
for r_ndx = 1:length(TOPset)
    rad_ndx = TOPset(r_ndx);
    rad = R00(rad_ndx);
    for p_ndx = 1:length(pov_lst)
        pov_ndx = pov_lst(p_ndx);
        pov = povs(pov_ndx);
        
        T = Tincident{pov_ndx};
        R = Rincident{pov_ndx}(:,rad_ndx);
        TN = Treceived{pov_ndx};
        PN = Psimulated{pov_ndx}(:,rad_ndx);
        
        ROI_LIST = find(TN>=sg_window_beg & TN<=sg_window_end);
        [resp,freq] = freq_spectrum(PN(ROI_LIST)/Pac,probe_fs,'abs');
        
        fig=figure(101); fig.Position=[680 80 980 1024];
        subplot(length(TOPset),3,(r_ndx-1)*3+1),
        plot(T*1e6,R*1e6,'Color',colors{p_ndx});
        grid on; hold on;
        xlim([ix_window_beg,ix_window_end+6e-6]*1e6);
        ylim([rad*0.5,rad+2e-6]*1e6);
        YL = ylim();
        rectangle('Position',[ix_window_beg*1e6,YL(1),(ix_window_end-ix_window_beg)*1e6,YL(2)-YL(1)],'EdgeColor','g');
%         xlabel('Time (\mus)');
%         ylabel('Radius (\mum)');
        
        subplot(length(TOPset),3,(r_ndx-1)*3+2),
        plot(TN*1e6,PN,'Color',colors{p_ndx});
        grid on; hold on;
        ylim([-60,60])
        YL = ylim();
        rectangle('Position',[sg_window_beg*1e6,YL(1),(sg_window_end-sg_window_beg)*1e6,YL(2)-YL(1)],'EdgeColor','g');
%         xlabel('Time (\mus)');
%         ylabel('Pressure (Pa)');
        
        ax = subplot(length(TOPset),3,(r_ndx-1)*3+3);
        plot(freq(2:end)/1e6,20*log10(resp(2:end)),'Color',colors{p_ndx});
        grid on; hold on;
        xlim([0,7]);
        ax.XTick = (0:Frq/1e6/2:7);
%         xlabel('Freq (MHz)');
%         ylabel('Amplitude (dB)');
        
    end
end

%% 显示每个微泡的次谐波响应 -- BOTTOM区域
for r_ndx = 1:length(BOTset)
    rad_ndx = BOTset(r_ndx);
    rad = R00(rad_ndx);
    for p_ndx = 1:length(pov_lst)
        pov_ndx = pov_lst(p_ndx);
        pov = povs(pov_ndx);
        
        T = Tincident{pov_ndx};
        R = Rincident{pov_ndx}(:,rad_ndx);
        TN = Treceived{pov_ndx};
        PN = Psimulated{pov_ndx}(:,rad_ndx);
        
        ROI_LIST = find(TN>=sg_window_beg & TN<=sg_window_end);
        [resp,freq] = freq_spectrum(PN(ROI_LIST)/Pac,probe_fs,'abs');
        
        fig=figure(102); fig.Position=[680 80 980 1024];
        subplot(length(BOTset),3,(r_ndx-1)*3+1),
        plot(T*1e6,R*1e6,'Color',colors{p_ndx});
        grid on; hold on;
        xlim([ix_window_beg,ix_window_end+6e-6]*1e6);
        ylim([rad*0.5,rad+2e-6]*1e6);
        YL = ylim();
        rectangle('Position',[ix_window_beg*1e6,YL(1),(ix_window_end-ix_window_beg)*1e6,YL(2)-YL(1)],'EdgeColor','g');
%         xlabel('Time (\mus)');
%         ylabel('Radius (\mum)');
        
        subplot(length(BOTset),3,(r_ndx-1)*3+2),
        plot(TN*1e6,PN,'Color',colors{p_ndx});
        grid on; hold on;
        ylim([-60,60])
        YL = ylim();
        rectangle('Position',[sg_window_beg*1e6,YL(1),(sg_window_end-sg_window_beg)*1e6,YL(2)-YL(1)],'EdgeColor','g');
%         xlabel('Time (\mus)');
%         ylabel('Pressure (Pa)');
        
        ax = subplot(length(BOTset),3,(r_ndx-1)*3+3);
        plot(freq(2:end)/1e6,20*log10(resp(2:end)),'Color',colors{p_ndx});
        grid on; hold on;
        xlim([0,7]);
        ax.XTick = (0:Frq/1e6/2:7);
%         xlabel('Freq (MHz)');
%         ylabel('Amplitude (dB)');
        
    end
end
