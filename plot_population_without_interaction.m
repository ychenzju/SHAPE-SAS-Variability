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

%% plot signals
caselab = {};
SHpsd = zeros(6,2);
SHdem = zeros(6,2);
sh = SHdem;
sh_pi = sh;
sh_am = sh;
sh_cps = sh;

%% Plot the comparison between existing and non-existing interactions
%% the user must choose the objects to be compared
%% this can be performed by commenting and uncommenting the following codes
fname1 = 'population_003_Povs_noX.mat'; %polydisperse population without interactions
fname2 = 'population_003_Povs.mat'; %polydisperse population with interactions
% fname1 = 'population_003_Povs_R2_noX.mat'; %monodisperse population of 2um w/o interactions
% fname2 = 'population_003_Povs_R2.mat'; %monodisperse population of 2um w/ interactions
% fname1 = 'population_003_Povs_R3_noX.mat'; %monodisperse population of 3um w/o interactions
% fname2 = 'population_003_Povs_R3.mat'; %monodisperse population of 3um w/ interactions

%% -------------------------------------------------------------------
%% Followed are simulation results for SonoVue microbubbles using Marmottant model
%% Data filenames align with the Marmottant parameter settings in population2.m 
% Marmottant1,   constant Kappa_s = 15e-9, Pac = 200kPa
% Marmottant2,   Size-dependent effect modeled, with Kappa_s = 10^(-9 + 0.6*R0); Pac = 200kPa
% Marmottant2_1, Shear-thinning effect modeled, with Kappa_s_0 = 0.32e-9; Pac = 200kPa
% Marmottant2_2, Shear-thinning effect modeled, with Kappa_s_0 = 75e-9; Pac = 200kPa
% Marmottant2_3, Shear-thinning effect modeled, with Kappa_s_0 = 23e-9; Pac = 200kPa
%% --------------------------------------------------------------------
% fname1 = 'population_003_Povs_noX_Marmottant5_3.mat'; %polydisperse population without interactions
% fname2 = 'population_003_Povs_Marmottant5_3.mat'; %polydisperse population with interactions


caselab{1} = ['Without Interaction'];
caselab{2} = ['With Interaction'];    
for case_ndx = 1:2
    switch case_ndx
        case 1
            load(fname1,'R*','Zs','probe_fs','frqs','povs','sg_window*','*received*','*incident*','Psimulated');
        case 2
            load(fname2,'R*','Zs','probe_fs','frqs','povs','sg_window*','*received*','*incident*','Psimulated');
    end
    
    pov_ndx = 1;
    
    %% 显示每个微泡的振动和散射信号
    for rad_ndx = 1:length(R00)
        rad = Rinit(rad_ndx,pov_ndx);
        fig=figure(2006); fig.Position = [200,40,1000,1080];
        subplot(1,2,1);
        X = Tincident{pov_ndx}*1e6;
        Y = (Rincident{pov_ndx}(:,rad_ndx)-rad)*1e6+rad_ndx;
        if case_ndx==1, lineStyle='-'; else lineStyle=':'; end
        ln = plot(X,Y,'Color',colors{rad_ndx},'LineStyle',lineStyle);
        text(X(1),Y(1),['R = ',num2str(rad*1e6,'%1.2f'),' '],'HorizontalAlignment','right','VerticalAlignment','middle');
        hold on;
        subplot(1,2,2);
        X = Treceived{pov_ndx}*1e6;
        Y = (Psimulated{pov_ndx}(:,rad_ndx))+rad_ndx*20;
        if case_ndx==1, lineStyle='-'; else lineStyle=':'; end
        ln = plot(X,Y,'Color',colors{rad_ndx},'LineStyle',lineStyle);
        text(X(1),Y(1),['#',num2str(rad_ndx),' '],'HorizontalAlignment','right','VerticalAlignment','middle');
        hold on;
    end
    ax = subplot(121);
    XL = xlim();
    xlim([XL(1)-2, XL(2)]);
    title('Radius Dynamics');
    xlabel('Time (\mus)');
    ax.YTickLabel = [];
    %YL = ylim();
    %rectangle('Position',[ix_window_beg*1e6,YL(1),(ix_window_end-ix_window_beg)*1e6,YL(2)-YL(1)],'EdgeColor','k');
    ax = subplot(122);
    XL = xlim();
    xlim([XL(1), XL(2)]);
    YL = ylim();
    rectangle('Position',[sg_window_beg*1e6,YL(1),(sg_window_end-sg_window_beg)*1e6,YL(2)-YL(1)],...
        'EdgeColor','g','LineWidth',1.4,'LineStyle','--');
    title('Received Acoustic Pressures')
    xlabel('Time (\mus)');
    ylim([0,YL(2)]);
    ax.YTickLabel = [];
    
    %% 显示#14 ~ #20微泡振动模式在有无相互作用下的区别
    if case_ndx==1
        bub_indices = (15:19)'; %(23:27)'; %
        fig = figure(2005); fig.Position = [100,100,1000,1000];
        subplot(3,2,1);
        XX = zeros(length(bub_indices),1);
        ZZ = Zs(bub_indices)*1e2;
        RR = R00(bub_indices)*1e6;
        scatter([0;0],Zs([bub_indices(1)-1,bub_indices(end)+1])*1e2,R00([bub_indices(1)-1,bub_indices(end)+1])*1e6*4,'b'); 
        hold on;
        scatter(XX,ZZ,RR*4,'r');
        grid on;
        title('Bubbles To Be Studied');
        xlabel('X');
        ylabel('Z (cm)');
        text(XX,ZZ,...
            compose("    #%2.0f, R = %1.2f \x03BCm",bub_indices,RR),...
            'FontSize',10,'HorizontalAlignment','left','VerticalAlignment','middle');
    end
    
    pov_ndx = 6;
    figure(2005); 
    for ndx = 1:length(bub_indices)
        bub_ndx = bub_indices(ndx);
        X = Tincident{pov_ndx} * 1e6;
        Y = Rincident{pov_ndx}(:,bub_ndx)*1e6;
        subplot(3,2,ndx+1)
        plot(X,Y); hold on; grid on;
        xlabel('Time (\mus)');
        ylabel('Radius Dynamics (\mum)');
        title(['Bubble #',num2str(bub_ndx)]);
        if case_ndx==2
            legend(caselab);
        end
        xlim([29,36]);
    end
    
    
    frq_ndx = 1; Frq = frqs(frq_ndx);
    
    for pov_ndx = 1 : length(povs)
        fig = figure(2001); fig.Position = [700,100,1000,1000];
        subplot(length(povs),3,[(pov_ndx-1)*3+1,(pov_ndx-1)*3+2]);
        X = Treceived{pov_ndx} * 1e6;
        Y = Preceived{pov_ndx};
        plot(X,Y); grid on; hold on;
        ylim([-150,150]);
        YL = ylim();
        rectangle('Position',[sg_window_beg*1e6,YL(1),(sg_window_end-sg_window_beg)*1e6,YL(2)-YL(1)],'LineStyle',':','EdgeColor','g','LineWidth',1.4);
        xlabel('Time (\mus)');
        ylabel({['\bfPov = ',num2str(povs(pov_ndx)/1e3),' kPa\rm'],'Pressure (Pa)'});
        if pov_ndx==1
            if contains(fname1,'R2')
                title_str = 'Received Acoustic Signal of 2\mum Monodisperse Population';
            elseif contains(fname1,'R3')
                title_str = 'Received Acoustic Signal of 3\mum Monodisperse Population';
            else
                title_str = 'Received Acoustic Signal of Reference Polydisperse Population';
            end
            title(title_str);
        end
        xlim([X(1),X(end)]);
        
        Tt = Treceived{pov_ndx};
        Pt = Preceived{pov_ndx};
        ROI_LIST = (Tt>=sg_window_beg) & (Tt<=sg_window_end);
        signal = Pt(ROI_LIST);
%         [resp,freq] = freq_spectrum(signal,probe_fs,'abs');
%         
%         % [pows,freq] = periodogram(signal,rectwin(length(signal)),length(signal),probe_fs,'onesided','power'); resp = (pows*2).^(1/2);
%         % [pows,freq] = pwelch(signal,length(signal),0,length(signal),probe_fs,'onesided','power'); resp = (pows*2).^(1/2);
%         % relation between power spectrum density and power spectrum: psd = pow / (fs * N);
%         ax = subplot(length(povs),3,(pov_ndx-1)*3+2);
%         X = freq(2:end) / 1e6; Y = 20*log10(resp(2:end));
%         plot(X,Y,'-o','MarkerSize',3); grid on; hold on;
%         ylim([-40,40]);
%         YL = ylim();
%         %rectangle('Position',[frqs(frq_ndx)/2/1e6-0.2,YL(1),0.4,YL(2)-YL(1)],'LineStyle',':','EdgeColor','r','LineWidth',1.4);
%         xlabel('Frequency (MHz)');
%         ylabel('Amplitude(dB)');
%         if pov_ndx==1
%             title('PSD-derived Amplitude Spectrum');
%         end
%         xlim([0,7.5]);
%         ax.XTick = 0:frqs(frq_ndx)/2/1e6:7.5;
        
        subplot(length(povs),3,(pov_ndx-1)*3+3);
        SHA = rf2iq_filter(Pt,probe_fs,Frq/2);
        ROI_LIST = (Tt>=sg_window_beg) & (Tt<=sg_window_end);
        sh_meanvalue = mean(SHA(ROI_LIST));
        X = Tt * 1e6; Y = 20*log10(SHA);
        plot(X,Y); grid on; hold on;
        ylim([-20,40]);
        YL = ylim();
        rectangle('Position',[sg_window_beg*1e6,YL(1),(sg_window_end-sg_window_beg)*1e6,YL(2)-YL(1)],'LineStyle',':','EdgeColor','g','LineWidth',1.4);
        %line([sg_window_beg;sg_window_end]*1e6,[sh_meanvalue;sh_meanvalue],'Color','r','LineStyle',':','LineWidth',1.4);
        xlabel('Time (\mus)');
        ylabel('Amplitude (dB)');
        if pov_ndx==1
            title('HT-demodulated Subharmonic Signal');
        end
        xlim([X(1),X(end)]);
        
        %% calculate the subharmonic amplitude vs. ambient pressure
        Tt = Treceived{pov_ndx};
        Pt = Preceived{pov_ndx};
        
        ROI_LIST = (Tt>=sg_window_beg) & (Tt<=sg_window_end);
        signal = Pt(ROI_LIST);
        [resp,freq] = freq_spectrum(signal,probe_fs,'dBW');
        sh_fc = frqs(frq_ndx)/2;
        sh_bw = 0.2e6;
        FOCUS_LIST = freq>=(sh_fc-sh_bw) & freq<=(sh_fc+sh_bw);
        SHpsd(pov_ndx,case_ndx) = mean(resp(FOCUS_LIST));

        SHA = rf2iq_filter(Pt,probe_fs,Frq/2);
        ROI_LIST = (Tt>=sg_window_beg) & (Tt<=sg_window_end);
        SHdem(pov_ndx,case_ndx) = mean(20*log10(SHA(ROI_LIST)));

        sh(pov_ndx,case_ndx) = mean(20*log10(SHA(ROI_LIST)));

        Pt = PreceivedPI{pov_ndx};
        SHA = rf2iq_filter(Pt,probe_fs,Frq/2);
        ROI_LIST = (Tt>=sg_window_beg) & (Tt<=sg_window_end);
        sh_pi(pov_ndx,case_ndx) = mean(20*log10(SHA(ROI_LIST)));

        Pt = PreceivedAM{pov_ndx};
        SHA = rf2iq_filter(Pt,probe_fs,Frq/2);
        ROI_LIST = (Tt>=sg_window_beg) & (Tt<=sg_window_end);
        sh_am(pov_ndx,case_ndx) = mean(20*log10(SHA(ROI_LIST)));

        Pt = PreceivedCPS{pov_ndx};
        SHA = rf2iq_filter(Pt,probe_fs,Frq/2);
        ROI_LIST = (Tt>=sg_window_beg) & (Tt<=sg_window_end);
        sh_cps(pov_ndx,case_ndx) = mean(20*log10(SHA(ROI_LIST)));
         
    end
    
    figure(2002);
    pov_ndx = 1;
    bub_ndx = 17;
    subplot(2,2,1);
    %plot the #16 bubble for the difference between interaction and non-interaction
    X = Tincident{pov_ndx} * 1e6;
    Y = Rincident{pov_ndx}(:,bub_ndx) * 1e6;
    plt = plot(X,Y); grid on; hold on;
    xlabel('Time (\mus)');
    ylabel('Radius (\mum)');
    title(['Radius Dynamics of Bubble #',num2str(bub_ndx)]);
    
    subplot(2,2,2);
    X = Treceived{pov_ndx} * 1e6;
    Y = Preceived{pov_ndx};
    plot(X,Y); grid on; hold on;
    xlabel('Time (\mus)');
    ylabel('Pressure (Pa)');
    title('ROI-windowed Signal From Population');
    xlim([sg_window_beg,sg_window_end]*1e6);

    figure(2003);
    pov_ndx = 1;
    subplot(2,2,1);
    %plot the #16 bubble for the difference between interaction and non-interaction
    X = Tincident{pov_ndx} * 1e6;
    Y = Rincident{pov_ndx}(:,bub_ndx) * 1e6;
    plt = plot(X,Y); grid on; hold on;
    xlabel('Time (\mus)');
    ylabel('Radius (\mum)');
    title(['Radius Dynamics of Bubble #',num2str(bub_ndx)]);
    xlim([X(1),X(end)]);
    
    subplot(2,2,2);
    X = Treceived{pov_ndx} * 1e6;
    Y = Preceived{pov_ndx};
    plot(X,Y); grid on; hold on;
    xlabel('Time (\mus)');
    ylabel('Pressure (Pa)');
    title('ROI-Windowed Signal From Population');
    xlim([sg_window_beg,sg_window_end]*1e6);

end

%% plot the subharmonic amplitude vs. ambient pressure
fig = figure(2001); fig.Position = [700,40,1000,1080];
subplot(length(povs),3,[1,2]); legend(caselab,'Location','southeast');
subplot(length(povs),3,3); legend(caselab,'Location','southeast');

fig = figure(2002); fig.Position = [100,100,1000,800];
subplot(2,2,1);
legend(caselab);
subplot(2,2,2);
legend(caselab);

subplot(2,2,3);
X = povs / 1e3;
Y = SHpsd;
plt = plot(X,Y); grid on;
legend(caselab,'Location','southwest');
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
end
xlabel('Ambient Pressure (kPa)');
ylabel('Amplitude (dB)');
title({'Subharmonic Response','PSD Analysis Method'});
subplot(2,2,4);
X = povs / 1e3;
Y = SHpsd-ones(size(SHpsd,1),1)*SHpsd(1,:);
plt = plot(X,Y); grid on;
legend(caselab,'Location','southwest');
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
end
xlabel('Ambient Pressure (kPa)');
ylabel('Amplitude (dB)');
title({'Subharmonic Response (Relative Variation)','PSD Analysis Method'});

fig = figure(2003); fig.Position = [600,200,1000,800];
subplot(2,2,1);
legend(caselab);
subplot(2,2,2);
legend(caselab);


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
title({'Subharmonic Response from Population','HT Demodulation Method'});
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
title({'Relative Subharmonic Response from Population','HT Demodulation Method'});

%% 显示四种成像模式下的次谐波-环境压力响应
fig=figure(2004); fig.Position = [100,40,1000,1080];

subplot(4,2,1);
X = povs / 1e3;
Y = sh;
plt = plot(X,Y); grid on;
legend(caselab,'Location','southwest','FontSize',9.5);
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
end
xlabel('Ambient Pressure (kPa)');
ylabel('Amplitude (dB)');
radii = unique(R00);
if length(radii)>1
    %title({'Polydisperse Microbubble Population','Subharmonic Response in Linear Imaging'});
    title({'Polydisperse Microbubble Population','Subharmonic Response'})
else
%     title({['Monodisperse Microbubble Population of ',num2str(radii*1e6),' \mum'],'Subharmonic Response in Linear Imaging'});
    title({['Monodisperse Microbubble Population of ',num2str(radii*1e6),' \mum'],'Subharmonic Response'});
end
subplot(4,2,2);
X = povs / 1e3;
Y = sh-ones(size(sh,1),1)*sh(1,:);
plt = plot(X,Y); grid on;
legend(caselab,'Location','southwest','FontSize',9.5);
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
end
xlabel('Ambient Pressure (kPa)');
ylabel('Amplitude (dB)');
%title({'Subharmonic Variation in Linear Imaging'});

subplot(4,2,3);
X = povs / 1e3;
Y = sh_pi;
plt = plot(X,Y); grid on;
legend(caselab,'Location','southwest','FontSize',9.5);
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
end
xlabel('Ambient Pressure (kPa)');
ylabel('Amplitude (dB)');
% title({'Subharmonic Response in PI Imaging'});
subplot(4,2,4);
X = povs / 1e3;
Y = sh_pi-ones(size(sh_pi,1),1)*sh_pi(1,:);
plt = plot(X,Y); grid on;
legend(caselab,'Location','southwest','FontSize',9.5);
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
end
xlabel('Ambient Pressure (kPa)');
ylabel('Amplitude (dB)');
%title({'Subharmonic Variation in PI Imaging'});

subplot(4,2,5);
X = povs / 1e3;
Y = sh_am;
plt = plot(X,Y); grid on;
legend(caselab,'Location','southwest','FontSize',9.5);
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
end
xlabel('Ambient Pressure (kPa)');
ylabel('Amplitude (dB)');
%title({'Subharmonic Response in AM Imaging'});
subplot(4,2,6);
X = povs / 1e3;
Y = sh_am-ones(size(sh_am,1),1)*sh_am(1,:);
plt = plot(X,Y); grid on;
legend(caselab,'Location','southwest','FontSize',9.5);
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
end
xlabel('Ambient Pressure (kPa)');
ylabel('Amplitude (dB)');
%title({'Subharmonic Variation in AM Imaging'});

subplot(4,2,7);
X = povs / 1e3;
Y = sh_cps;
plt = plot(X,Y); grid on;
legend(caselab,'Location','southwest','FontSize',9.5);
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
end
xlabel('Ambient Pressure (kPa)');
ylabel('Amplitude (dB)');
%title({'Subharmonic Response in CPS Imaging'});

subplot(4,2,8);
X = povs / 1e3;
Y = sh_cps-ones(size(sh_cps,1),1)*sh_cps(1,:);
plt = plot(X,Y); grid on;
legend(caselab,'Location','southwest','FontSize',9.5);
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
end
xlabel('Ambient Pressure (kPa)');
ylabel('Amplitude (dB)');
%title({'Subharmonic Variation in CPS Imaging'});
