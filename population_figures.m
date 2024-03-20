%% plot figures in bubble.m
% close all

%% plot out for quick check
colors = {'b','r','k','m','g','c'};
styles = {'-','-.',':','--'};
mrkers = {'*','o','s','d','v','^','+','<','>','x','p','none'};
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

%% plot population response --- Amplitude vs. overpressure
if length(rads)>1 && length(povs)>1 && length(frqs)==1 && length(pacs)==1
    if exist('Param','var') && isfield(Param,'use_energy_spectrum') && Param.use_energy_spectrum==1
        xsh = reshape(ysh,length(povs),1);
        xfd = reshape(yfd,length(povs),1);
        xuh = reshape(yuh,length(povs),1);
        xhm = reshape(yhm,length(povs),1);
%         xfq = reshape(yfq,length(povs),1)/1e6;
%         xmg = 20*log10(reshape(ymg,length(povs),1));
    else
        xsh = reshape(xsh,length(povs),1);
        xfd = reshape(xfd,length(povs),1);
        xuh = reshape(xuh,length(povs),1);
        xhm = reshape(xhm,length(povs),1);
%         xfq = reshape(xfq,length(povs),1)/1e6;
%         xmg = 20*log10(reshape(xmg,length(povs),1));
    end
    
    fig=figure(40019); fig.Position = [680 610 560 480];
    plt = plot(povs/1e3,[xsh(:),xfd(:),xuh(:),xhm(:)] - ones(length(xsh),1)*[xsh(1),xfd(1),xuh(1),xhm(1)],'-*'); grid on; hold on;
    for k=1:length(plt)
        plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
    end
    legend(['Subharmonic, \Delta = ',num2str(xsh(end)-xsh(1),'%2.2f'),' dB'],...
        ['Fundamental, \Delta = ',num2str(xfd(end)-xfd(1),'%2.2f'),' dB'],...
        ['Ultraharmonic, \Delta = ',num2str(xuh(end)-xuh(1),'%2.2f'),' dB'],...
        ['2^{nd}-Harmonic, \Delta = ',num2str(xhm(end)-xhm(1),'%2.2f'),' dB'],'Location','southwest');
    xlabel('Ambient Overpressure (kPa)'); ylabel('Amplitude (dB)');
    title({'Harmonics Variation vs. Ambient Overpressure',...
        ['Microbubble Population, Model = ', bubble_model],...
        ['Driving Frequency = ',num2str(frqs/1e6),' MHz, ', 'Pulse Magnitude = ',num2str(pacs/1e3),' kPa']});
    

    fig=figure(4001); fig.Position = [680 610 1120 480];
    subplot(1,2,1);
    plt = plot(povs/1e3,[xsh(:),xfd(:),xuh(:),xhm(:)],'-*'); grid on; hold on;
    for k=1:length(plt)
        plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
    end
    legend(['Subharmonic  ',num2str((xsh(end)-xsh(1)),'%2.2f'),' dB'],...
        ['Fundamental  ',num2str(xfd(end)-xfd(1),'%2.2f'),' dB'],...
        ['Ultraharmonic  ',num2str(xuh(end)-xuh(1),'%2.2f'),' dB'],...
        ['2^{nd}-Harmonic  ',num2str(xhm(end)-xhm(1),'%2.2f'),' dB'],'Location','southwest');
    xlabel('Ambient Overpressure(kPa)'); ylabel('Amplitude(dB)');
    title({'Harmonics Amplitude vs. Ambient Overpressure',...
        ['Microbubble Population, Model = ', bubble_model],...
        ['Driving Frequency = ',num2str(frqs/1e6),' MHz, ', 'Pulse Magnitude = ',num2str(pacs/1e3),'kPa']});

    subplot(1,2,2);
    plt = plot(povs/1e3,[xsh(:),xfd(:),xuh(:),xhm(:)] - ones(length(xsh),1)*[xsh(1),xfd(1),xuh(1),xhm(1)],'-*'); grid on; hold on;
    for k=1:length(plt)
        plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
    end
    legend(['Subharmonic  ',num2str(xsh(end)-xsh(1),'%2.2f'),' dB'],...
        ['Fundamental  ',num2str(xfd(end)-xfd(1),'%2.2f'),' dB'],...
        ['Ultraharmonic  ',num2str(xuh(end)-xuh(1),'%2.2f'),' dB'],...
        ['2^{nd}-Harmonic  ',num2str(xhm(end)-xhm(1),'%2.2f'),' dB'],'Location','southwest');
    xlabel('Ambient Overpressure(kPa)'); ylabel('Amplitude(dB)');
    title({'Harmonics Variation vs. Ambient Overpressure',...
        ['Microbubble Population, Model=', bubble_model],...
        ['Driving Frequency=',num2str(frqs/1e6),'MHz, ', 'Pulse Magnitude=',num2str(pacs/1e3),'kPa']});

    if exist('PreceivedPI','var') && exist('PreceivedCPS','var')
        
        fig=figure(4008); fig.Position = [100 100 1000 800];
        subplot(2,2,1);
        plt = plot(povs/1e3,[xsh(:),xfd(:),xuh(:),xhm(:)],'-*'); grid on; hold on;
        for k=1:length(plt)
            plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
        end
        legend(['Subharmonic, \Delta = ',num2str((xsh(end)-xsh(1)),'%2.2f'),' dB'],...
            ['Fundamental, \Delta = ',num2str(xfd(end)-xfd(1),'%2.2f'),' dB'],...
            ['Ultraharmonic, \Delta = ',num2str(xuh(end)-xuh(1),'%2.2f'),' dB'],...
            ['2^{nd}-Harmonic, \Delta = ',num2str(xhm(end)-xhm(1),'%2.2f'),' dB'],'Location','southwest');
        xlabel('Ambient Overpressure(kPa)'); ylabel('Amplitude(dB)');
        title({'Harmonics Amplitude vs. Ambient Overpressure',...
            ['Linear Imaging, Model = ', bubble_model],...
            ['Driving Frequency = ',num2str(frqs/1e6),' MHz, ', 'Pulse Magnitude = ',num2str(pacs/1e3),' kPa']});

        subplot(2,2,2);
        plt = plot(povs/1e3,[xsh_pi(:),xfd_pi(:),xuh_pi(:),xhm_pi(:)],'-*'); grid on; hold on;
        for k=1:length(plt)
            plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
        end
        legend(['Subharmonic, \Delta = ',num2str(xsh_pi(end)-xsh_pi(1),'%2.2f'),' dB'],...
            ['Fundamental, \Delta = ',num2str(xfd_pi(end)-xfd_pi(1),'%2.2f'),' dB'],...
            ['Ultraharmonic, \Delta = ',num2str(xuh_pi(end)-xuh_pi(1),'%2.2f'),' dB'],...
            ['2^{nd}-Harmonic, \Delta = ',num2str(xhm_pi(end)-xhm_pi(1),'%2.2f'),' dB'],'Location','southwest');
        xlabel('Ambient Overpressure(kPa)'); ylabel('Amplitude(dB)');
        title({'Harmonics Amplitude vs. Ambient Overpressure',...
            ['PI-mode Imaging, Model = ', bubble_model],...
            ['Driving Frequency = ',num2str(frqs/1e6),' MHz, ', 'Pulse Magnitude = ',num2str(pacs/1e3),' kPa']});

        subplot(2,2,3);
        plt = plot(povs/1e3,[xsh_am(:),xfd_am(:),xuh_am(:),xhm_am(:)],'-*'); grid on; hold on;
        for k=1:length(plt)
            plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
        end
        legend(['Subharmonic, \Delta = ',num2str(xsh_am(end)-xsh_am(1),'%2.2f'),' dB'],...
            ['Fundamental, \Delta = ',num2str(xfd_am(end)-xfd_am(1),'%2.2f'),' dB'],...
            ['Ultraharmonic, \Delta = ',num2str(xuh_am(end)-xuh_am(1),'%2.2f'),' dB'],...
            ['2^{nd}-Harmonic, \Delta = ',num2str(xhm_am(end)-xhm_am(1),'%2.2f'),' dB'],'Location','southwest');
        xlabel('Ambient Overpressure(kPa)'); ylabel('Amplitude(dB)');
        title({'Harmonics Amplitude vs. Ambient Overpressure',...
            ['AM-mode Imaging, Model = ', bubble_model],...
            ['Driving Frequency = ',num2str(frqs/1e6),' MHz, ', 'Pulse Magnitude = ',num2str(pacs/1e3),' kPa']});

        subplot(2,2,4);
        plt = plot(povs/1e3,[xsh_cps(:),xfd_cps(:),xuh_cps(:),xhm_cps(:)],'-*'); grid on; hold on;
        for k=1:length(plt)
            plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
        end
        legend(['Subharmonic, \Delta = ',num2str(xsh_cps(end)-xsh_cps(1),'%2.2f'),' dB'],...
            ['Fundamental, \Delta = ',num2str(xfd_cps(end)-xfd_cps(1),'%2.2f'),' dB'],...
            ['Ultraharmonic, \Delta = ',num2str(xuh_cps(end)-xuh_cps(1),'%2.2f'),' dB'],...
            ['2^{nd}-Harmonic, \Delta = ',num2str(xhm_cps(end)-xhm_cps(1),'%2.2f'),' dB'],'Location','southwest');
        xlabel('Ambient Overpressure(kPa)'); ylabel('Amplitude(dB)');
        title({'Harmonics Amplitude vs. Ambient Overpressure',...
            ['CPS-mode Imaging, Model = ', bubble_model],...
            ['Driving Frequency = ',num2str(frqs/1e6),' MHz, ', 'Pulse Magnitude = ',num2str(pacs/1e3),' kPa']});
        
        
        fig=figure(4009); fig.Position = [100 200 1000 800];
        subplot(2,2,1);
        plt = plot(povs/1e3,[xsh(:),xfd(:),xuh(:),xhm(:)] - ones(length(xsh),1)*[xsh(1),xfd(1),xuh(1),xhm(1)],'-*'); grid on; hold on;
        for k=1:length(plt)
            plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
        end
        legend(['Subharmonic, \Delta = ',num2str(xsh(end)-xsh(1),'%2.2f'),' dB'],...
            ['Fundamental, \Delta = ',num2str(xfd(end)-xfd(1),'%2.2f'),' dB'],...
            ['Ultraharmonic, \Delta = ',num2str(xuh(end)-xuh(1),'%2.2f'),' dB'],...
            ['2^{nd}-Harmonic, \Delta = ',num2str(xhm(end)-xhm(1),'%2.2f'),' dB'],'Location','southwest');
        xlabel('Ambient Overpressure(kPa)'); ylabel('Amplitude(dB)');
        title({'Harmonics Variation vs. Ambient Overpressure',...
            ['Linear Imaging, Model = ', bubble_model],...
            ['Driving Frequency = ',num2str(frqs/1e6),' MHz, ', 'Pulse Magnitude = ',num2str(pacs/1e3),' kPa']});

        subplot(2,2,2);
        plt = plot(povs/1e3,[xsh_pi(:),xfd_pi(:),xuh_pi(:),xhm_pi(:)] - ones(length(xsh_pi),1)*[xsh_pi(1),xfd_pi(1),xuh_pi(1),xhm_pi(1)],'-*'); grid on; hold on;
        for k=1:length(plt)
            plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
        end
        legend(['Subharmonic, \Delta = ',num2str(xsh_pi(end)-xsh_pi(1),'%2.2f'),' dB'],...
            ['Fundamental, \Delta = ',num2str(xfd_pi(end)-xfd_pi(1),'%2.2f'),' dB'],...
            ['Ultraharmonic, \Delta = ',num2str(xuh_pi(end)-xuh_pi(1),'%2.2f'),' dB'],...
            ['2^{nd}-Harmonic, \Delta = ',num2str(xhm_pi(end)-xhm_pi(1),'%2.2f'),' dB'],'Location','southwest');
        xlabel('Ambient Overpressure(kPa)'); ylabel('Amplitude(dB)');
        title({'Harmonics Variation vs. Ambient Overpressure',...
            ['PI-mode Imaging, Model = ', bubble_model],...
            ['Driving Frequency = ',num2str(frqs/1e6),' MHz, ', 'Pulse Magnitude = ',num2str(pacs/1e3),' kPa']});

        subplot(2,2,3);
        plt = plot(povs/1e3,[xsh_am(:),xfd_am(:),xuh_am(:),xhm_am(:)] - ones(length(xsh_am),1)*[xsh_am(1),xfd_am(1),xuh_am(1),xhm_am(1)],'-*'); grid on; hold on;
        for k=1:length(plt)
            plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
        end
        legend(['Subharmonic, \Delta = ',num2str(xsh_am(end)-xsh_am(1),'%2.2f'),' dB'],...
            ['Fundamental, \Delta = ',num2str(xfd_am(end)-xfd_am(1),'%2.2f'),' dB'],...
            ['Ultraharmonic, \Delta = ',num2str(xuh_am(end)-xuh_am(1),'%2.2f'),' dB'],...
            ['2^{nd}-Harmonic, \Delta = ',num2str(xhm_am(end)-xhm_am(1),'%2.2f'),' dB'],'Location','southwest');
        xlabel('Ambient Overpressure(kPa)'); ylabel('Amplitude(dB)');
        title({'Harmonics Variation vs. Ambient Overpressure',...
            ['AM-mode Imaging, Model = ', bubble_model],...
            ['Driving Frequency = ',num2str(frqs/1e6),' MHz, ', 'Pulse Magnitude = ',num2str(pacs/1e3),' kPa']});

        subplot(2,2,4);
        plt = plot(povs/1e3,[xsh_cps(:),xfd_cps(:),xuh_cps(:),xhm_cps(:)] - ones(length(xsh_cps),1)*[xsh_cps(1),xfd_cps(1),xuh_cps(1),xhm_cps(1)],'-*'); grid on; hold on;
        for k=1:length(plt)
            plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
        end
        legend(['Subharmonic, \Delta = ',num2str(xsh_cps(end)-xsh_cps(1),'%2.2f'),' dB'],...
            ['Fundamental, \Delta = ',num2str(xfd_cps(end)-xfd_cps(1),'%2.2f'),' dB'],...
            ['Ultraharmonic, \Delta = ',num2str(xuh_cps(end)-xuh_cps(1),'%2.2f'),' dB'],...
            ['2^{nd}-Harmonic, \Delta = ',num2str(xhm_cps(end)-xhm_cps(1),'%2.2f'),' dB'],'Location','southwest');
        xlabel('Ambient Overpressure(kPa)'); ylabel('Amplitude(dB)');
        title({'Harmonics Variation vs. Ambient Overpressure',...
            ['CPS Imaging, Model = ', bubble_model],...
            ['Driving Frequency = ',num2str(frqs/1e6),' MHz, ', 'Pulse Magnitude = ',num2str(pacs/1e3),' kPa']});
        
    end
    
    fig=figure(40011); fig.Position = [180,610,560,480];
    ax = axes();
    if length(povs)==6
        list = [1,3,6];
    elseif length(povs)==3
        list = [1,2,3];
    else
        list = [1];
    end        
    
    srecv = cell(length(povs),1);
    for pov_ndx=1:length(povs)
        tr = Treceived{pov_ndx};
        sglist = tr>=sg_window_beg & tr<sg_window_end;
        pr = Preceived{pov_ndx,1,1};
        pr_roi = pr(sglist);
        [pr_roi_fft,pr_roi_frq] = freq_spectrum(pr_roi/pacs(1),probe_fs,'abs');
        srecv{pov_ndx} = [pr_roi_frq(:),pr_roi_fft(:)];
    end
    ndx = 0;
    for pov_ndx = list
        ndx = ndx + 1;
        plt=plot(srecv{pov_ndx}(:,1)/1e6,20*log10(srecv{pov_ndx}(:,2)));
        plt.Color=colors{ndx}; plt.LineStyle=styles{ndx}; plt.Marker='none'; plt.LineWidth=widths{ndx};
        grid on; hold on;
    end
    legend(compose('Pov = %2.0f kPa', povs(list)/1e3));
    xlabel('Frequency (MHz)');
    ylabel('Magnitude (dB)')
    title({'Magnitude Spectrum of Received Acoustic Signal',...
        ['Microbubble Population, Model = ', bubble_model],...
        ['Driving Frequency = ',num2str(frqs/1e6),' MHz, ', 'Pulse Magnitude = ',num2str(pacs/1e3),' kPa']});
    xlim([0,10]);
    ax.XTick = 0:(frqs/1e6/2):10;
%     ylim([-80,-30]);
    
    
    if exist('Param','var') && isfield(Param,'use_energy_spectrum') && Param.use_energy_spectrum==1
        fig=figure(40012); fig.Position = [680 610 560 480];
        ax = axes();
        ndx = 0;
        for pov_ndx = list
            ndx = ndx + 1;
            plt=plot(yrecv{pov_ndx,1,1}(:,1)/1e6,20*log10(yrecv{pov_ndx,1,1}(:,2)));
            plt.Color=colors{ndx}; plt.LineStyle=styles{ndx}; plt.Marker='none'; plt.LineWidth=widths{ndx};
            grid on; hold on;
        end
        legend(compose('Pov = %2.0f kPa', povs(list)/1e3));
        xlabel('Frequency (MHz)');
        ylabel('Magnitude (dB)')
        title({'Magnitude Spectrum of Received Acoustic Signal',...
            ['Microbubble Population, Model = ', bubble_model],...
            ['Driving Frequency = ',num2str(frqs/1e6),' MHz, ', 'Pulse Magnitude = ',num2str(pacs/1e3),' kPa']});
        
        xlim([0,7.5]);
        ylim([-85,-60]);
        ax.XTick = 0:frqs/1e6/2:7.5;
    end
end

%% plot population response --- Amplitude vs. Driving Frequency
if length(rads)>1 && length(povs)==1 && length(frqs)>1 && length(pacs)==1
    if exist('Param','var') && isfield(Param,'use_energy_spectrum') && Param.use_energy_spectrum==1
        xsh = reshape(ysh,length(frqs),1);
        xfd = reshape(yfd,length(frqs),1);
        xuh = reshape(yuh,length(frqs),1);
        xhm = reshape(yhm,length(frqs),1);
    else
        xsh = reshape(xsh,length(frqs),1);
        xfd = reshape(xfd,length(frqs),1);
        xuh = reshape(xuh,length(frqs),1);
        xhm = reshape(xhm,length(frqs),1);
    end
    
    fig=figure(4002); fig.Position = [680 610 560 480];
    plt = plot(frqs/1e6,[xsh(:),xfd(:),xuh(:),xhm(:)],'-*'); grid on; hold on;
    for k=1:length(plt)
        plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
    end
    legend(['Subharmonic ',num2str(xsh(end)-xsh(1),'%2.2f'),' dB'],...
        ['Fundamental ',num2str(xfd(end)-xfd(1),'%2.2f'),' dB'],...
        ['Ultraharmonic ',num2str(xuh(end)-xuh(1),'%2.2f'),' dB'],...
        ['2^{nd}-Harmonic ',num2str(xhm(end)-xhm(1),'%2.2f'),' dB']);
    xlabel('Driving Frequency(MHz)'); ylabel('Amplitude(dB)');
    title({'MB-Population Harmonics vs. Driving Frequency',...
        ['Model=', bubble_model, ', Microbubble Population'],...
        ['OverPressure=',num2str(povs/1e3),'kPa, ', 'Pulse Magnitude=',num2str(pacs/1e3),'kPa']});
end

%% plot population response --- Amplitude vs. Pulse Magnitude
if length(rads)>1 && length(povs)==1 && length(frqs)==1 && length(pacs)>1
    if exist('Param','var') && isfield(Param,'use_energy_spectrum') && Param.use_energy_spectrum==1
        c_xsh = ysh; c_xfd = yfd; c_xuh = yuh; c_xhm = yhm;
    else
        c_xsh = xsh; c_xfd = xfd; c_xuh = xuh; c_xhm = xhm;
    end
    xsh = reshape(c_xsh,length(pacs),1);
    xfd = reshape(c_xfd,length(pacs),1);
    xuh = reshape(c_xuh,length(pacs),1);
    xhm = reshape(c_xhm,length(pacs),1);
    fig=figure(4003); fig.Position = [680 610 560 480];
    plt = plot(pacs/1e3,[xsh(:),xfd(:),xuh(:),xhm(:)],'-*'); grid on; hold on;
    for k=1:length(plt)
        plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
    end
    legend(['Subharmonic ',num2str(c_xsh(end)-c_xsh(1),'%2.2f'),' dB'],...
        ['Fundamental ',num2str(c_xfd(end)-c_xfd(1),'%2.2f'),' dB'],...
        ['Ultraharmonic ',num2str(c_xuh(end)-c_xuh(1),'%2.2f'),' dB'],...
        ['2^{nd}-Harmonic ',num2str(c_xhm(end)-c_xhm(1),'%2.2f'),' dB']);
    xlabel('Pulse Magnitude(kPa)'); ylabel('Amplitude(dB)');
    title({'MB-Population Harmonics vs. Driving Frequency',...
        ['Model=', bubble_model, ', Microbubble Population'],...
        ['Driving Frequency=',num2str(frqs/1e6),'kPa, ', 'OverPressure=',num2str(povs/1e3),'kPa']});
end

%% plot population response --- Amplitude va. (Pulse Magnitude & Driving Frequency)
if length(rads)>1 && length(povs)==1 && length(frqs)>1 && length(pacs)>1
    c_frqs = frqs; c_pacs = pacs;
    if exist('Param','var') && isfield(Param,'use_energy_spectrum') && Param.use_energy_spectrum==1
        c_xsh = ysh(1,:,:); c_xfd = yfd(1,:,:); c_xuh = yuh(1,:,:); c_xhm = yhm(1,:,:);
    else
        c_xsh = xsh(1,:,:); c_xfd = xfd(1,:,:); c_xuh = xuh(1,:,:); c_xhm = xhm(1,:,:);
    end
    c_xsh = reshape(c_xsh,length(c_frqs),length(c_pacs));
    c_xfd = reshape(c_xfd,length(c_frqs),length(c_pacs));
    c_xuh = reshape(c_xuh,length(c_frqs),length(c_pacs));
    c_xhm = reshape(c_xhm,length(c_frqs),length(c_pacs));
    
    fig=figure(5001); fig.Position = [580 50 560 480];
    plt = plot(c_pacs/1e3,c_xsh','Marker','o'); grid on;
    for k=1:length(plt)
        plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
    end
    legend(compose("f=%2.2f MHz",c_frqs(:)/1e6),'Location','southeast');
    xlabel('Pulse Magnitude(kPa)'); ylabel('Amplitude(dB)');
    title({'MB-Population Subharmonic vs. Pulse Magnitude',...
        ['Model=', bubble_model,', Microbubble Population'],...
        ['OverPressure=',num2str(povs/1e3),'kPa']});
    
    fig=figure(5002); fig.Position = [80 50 560 480];
    plt = plot(c_frqs/1e6,c_xsh','Marker','o'); grid on;
    for k=1:length(plt)
        plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
    end
    legend(compose("Pac=%2.2f kPa",c_pacs(:)/1e3),'Location','southeast');
    xlabel('Driving Frequency(MHz)'); ylabel('Amplitude(dB)');
    title({'MB-Population Subharmonic vs. Driving Frequency',...
        ['Model=', bubble_model,', Microbubble Population'],...
        ['OverPressure=',num2str(povs/1e3),'kPa']});
    
    fig=figure(50020); fig.Position = [80 50 560 480];
    lgd_strs = [];
    X = frqs;
    for rad_ndx = 1:length(rads)
        Y = reshape(sh(rad_ndx,:,:,:),length(c_frqs),length(c_pacs));
        plot(X(:)/1e6,Y); grid on; hold on;
        lgd_strs = cat(1, lgd_strs, compose("R0 = %2.2f \x03BCm, Pac = %2.f kPa",rads(rad_ndx)*1e6, pacs(:)/1e3));
    end
    plt = fig.Children.Children;
    for k=1:length(plt)
        plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).LineWidth=widths{k};% plt(k).Marker=mrkers{k}; plt(k).MarkerSize=4;
    end
    legend(lgd_strs,'Location','southeast');
    xlabel('Driving Frequency (MHz)'); ylabel('Subharmonic Amplitude (dB)');
    title({'Bubble Subharmonic vs. Driving Frequency',...
        ['Model = ', bubble_model,', OverPressure = ',num2str(povs/1e3),' kPa']});
    
end

%% plot population response --- Amplitude va. (Overpressure & Driving Frequency)
if length(rads)>1 && length(povs)>1 && length(frqs)>1 && length(pacs)==1
    c_frqs = frqs; c_povs = povs;
    if exist('Param','var') && isfield(Param,'use_energy_spectrum') && Param.use_energy_spectrum==1
        c_xsh = ysh(:,:,1); % c_xfd = yfd(:,:,1); c_xuh = yuh(:,:,1); c_xhm = yhm(:,:,1);
    else
        c_xsh = xsh(:,:,1); % c_xfd = xfd(:,:,1); c_xuh = xuh(:,:,1); c_xhm = xhm(:,:,1);
    end
    
    c_xsh = reshape(c_xsh,length(c_povs),length(c_frqs));
    %c_xfd = reshape(c_xfd,length(c_povs),length(c_frqs));
    %c_xuh = reshape(c_xuh,length(c_povs),length(c_frqs));
    %c_xhm = reshape(c_xhm,length(c_povs),length(c_frqs));
    
    fig=figure(5003); fig.Position =  [100 100 1000 400];
    subplot(1,2,1);
    plt = plot(c_frqs(:)/1e6,c_xsh'); grid on;
    for k=1:length(plt)
        plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
    end
    legend(compose("Pov = %2.0f kPa",c_povs(:)/1e3),'Location','southwest');
    xlabel('Driving Frequency (MHz)'); ylabel('Subharmonic Amplitude (dB)');
    title({'MB-Population Subharmonic Amplitude vs. Driving Frequency',...
        ['Microbubble Population, Model = ', bubble_model,],...
        ['Pulse Magnitude = ',num2str(pacs/1e3),' kPa']});
    

    subplot(1,2,2);
    frq_demo = [2.2,2.4,2.6,3.2,3.4,3.6,4.0,4.4]*1e6;
    [~,ia,ib] = intersect(c_frqs,frq_demo);
    plt = plot(c_povs/1e3,c_xsh(:,ia) - ones(size(c_xsh,1),1) * c_xsh(1,ia),'Marker','o'); grid on;
    for k=1:length(plt)
        plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
    end
    legend(compose("f = %2.2f MHz",c_frqs(ia)/1e6),'Location','southwest');
    xlabel('Ambient Overpressure (kPa)'); ylabel('Subharmonic Variation (dB)');
    title({'MB-Population Subharmonic Variation vs. Overpressure',...
        ['Model = ', bubble_model,', Microbubble Population'],...
        ['Pulse Magnitude = ',num2str(pacs/1e3),' kPa']});
    
end

%% plot population response --- Amplitude va. (Overpressure & Pulse Magnitude)
if length(rads)>1 && length(povs)>1 && length(frqs)==1 && length(pacs)>1
    c_pacs = pacs; c_povs = povs;
    if exist('Param','var') && isfield(Param,'use_energy_spectrum') && Param.use_energy_spectrum==1
        c_xsh = ysh(:,1,:); %c_xfd = yfd(:,1,:); c_xuh = yuh(:,1,:); c_xhm = yhm(:,1,:);
    else
        c_xsh = xsh(:,1,:); %c_xfd = xfd(:,1,:); c_xuh = xuh(:,1,:); c_xhm = xhm(:,1,:);
    end
    
    c_xsh = reshape(c_xsh,length(c_povs),length(c_pacs));
%     c_xfd = reshape(c_xfd,length(c_povs),length(c_pacs));
%     c_xuh = reshape(c_xuh,length(c_povs),length(c_pacs));
%     c_xhm = reshape(c_xhm,length(c_povs),length(c_pacs));

    fig=figure(5006); fig.Position =  [100 100 1000 400];
    subplot(1,2,1);
    plt = plot(c_pacs/1e3,c_xsh','Marker','o'); grid on;
    for k=1:length(plt)
        plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
    end
    legend(compose("Pov = %2.2f kPa",c_povs(:)/1e3),'Location','southeast');
    xlabel('Pulse Magnitude (kPa)'); ylabel('Subharmonic Amplitude (dB)');
    title({'MB-Population Subharmonic Amplitude vs. Pulse Magnitude',...
        ['Microbubble Population, Model = ', bubble_model],...
        ['Driving Frequency = ',num2str(frqs/1e6),' MHz']});
    
    subplot(1,2,2);
    pac_demo = [100,300,350,400,600]*1e3;
    [~,ia,ib] = intersect(c_pacs,pac_demo);
    
    plt = plot(c_povs/1e3,c_xsh(:,ia) - ones(size(c_xsh,1),1) * c_xsh(1,ia),'Marker','o'); grid on;
    for k=1:length(plt)
        plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
    end
    legend(compose("Pac=%2.2f kPa",c_pacs(ia)/1e3),'Location','southwest');
    xlabel('Overpressure(kPa)'); ylabel('Amplitude(dB)');
    title({'MB-Population Subharmonic Variation vs. Overpressure',...
        ['Model=', bubble_model,', Microbubble Population'],...
        ['Driving Frequency=',num2str(frqs/1e6),'MHz']});
    
%     fig=figure(50050); fig.Position = [580 50 560 480];
%     X = c_povs/1e3;
%     Y = c_xsh - c_xfd; Y = Y - ones(size(Y,1),1) * Y(1,:);
%     plt = plot(X,Y,'Marker','o'); grid on;
%     for k=1:length(plt)
%         plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
%     end
%     legend(compose("Pac = %2.0f kPa",c_pacs(:)/1e3),'Location','southeast');
%     xlabel('Overpressure (kPa)'); ylabel('Amplitude (dB)');
%     title({'MB-Population (SH-FD) variation vs. Overpressure',...
%         ['Microbubble Population',', Model = ', bubble_model],...
%         ['Driving Frequency = ',num2str(frqs/1e6),' MHz']});
    
end

