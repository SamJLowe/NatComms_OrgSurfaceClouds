%/----------------------------------X------------------------------------/%
%Sam Lowe, ACES (SU) 2019: Plotting Figure 4, Lowe et al. 2019, Nat. Comms
%/----------------------------------X------------------------------------/%

close all; clear all;
%Concentration axes
nRes   = 20;
N2_vec = logspace(1, 3, nRes); N1_vec = logspace(0, 4, nRes);

%Ouput strings
OUTPUT_DATA_STR{1} = 'Marine_OPTICS_N1N2_DELTACDNC.mat';
OUTPUT_DATA_STR{2} = 'REV_HYY_Boreal_OPTICS_N1N2_DELTACDNC.mat';
OUTPUT_DATA_STR{3} = 'NUMevent_OPTICS_N1N2_DELTACDNC.mat';
ALBEDO5_DELTACDNC  = [22 28 39];

%For plotting-------------------------------------------------------------
ncont        = 5;
COL_SCALE{1} = 'dense';   COLSENV{1} = 'c';
COL_SCALE{2} = 'tempo';   COLSENV{2} = 'g';
COL_SCALE{3} = 'amp';     COLSENV{3} = 'r';
COLBAR_Y0    = 0.015;    COLBAR_LEN = 0.22; %Color bar dimensions
XTick        = [1 10 100 1000 10000]; XTicklabs  = {'1'; '10'; '100';'1,000';'10,000'};
YTick        = [10 100 1000];         YTicklabs  = {'10'; '100'; '1,000'};
MK_EDGE_WDTH = 1;
CHAR_MKSZ    = 40;
NI_MKSZ      = 8;
BOX_LIN_WDTH = 1.5;
CONT_WIDTH   = 1.2;
CLIM         = [0 35; 0 80];      %Delta CDNC
CLIM_TAU     = [0 10; 0 22];
CLIM_ALPHA   = [0 8; 0 15];
COLORS       = [0.00, 0.45, 0.74;...
                0.47, 0.67, 0.19;...
                0.85, 0.33, 0.1];
ALPHABET     = char(97:122);
ENV{1} = 'MAV'; ENV{2} = 'HYY';
%--------------------------------------------------------------------------
%Get aerosol concnetration and Delta cdnc for MH time series
[MH_MODE_CONC, DELTA_CDNC_MH, DELTA_CDNC_MH_THRES] = GETN1N2_MH;
save('Figure4a_scatter.mat','MH_MODE_CONC')

return
ii = 1
figure(1)
%set(findall(gcf,'-property','FontSize'),'FontSize',14)
for i = 1:2
    load(OUTPUT_DATA_STR{i});
    DELTA_CDNC = (CDNC(:,:,2)./CDNC(:,:,1) - 1) *100; DIFF_CDNC  = CDNC(:,:,2) - CDNC(:,:,1);
    DELTA_SMAX = (SMAX(:,:,2)./SMAX(:,:,1) - 1) *100; DIFF_SMAX  = abs(SMAX(:,:,2) - SMAX(:,:,1));
    DELTA_RACT = squeeze((RACT_OUT(:,:,2,2)./RACT_OUT(:,:,1,2) - 1)) *100;
    DIFF_RACT  = -1 * squeeze(RACT_OUT(:,:,2,2) - RACT_OUT(:,:,1,2));
    DELTA_ALPHA = (ALPHA(:,:,2)./ALPHA(:,:,1) - 1) *100; DIFF_ALPHA = ALPHA(:,:,2) - ALPHA(:,:,1);
    DELTA_TAU = (TAU(:,:,2)./TAU(:,:,1) - 1) *100; DIFF_TAU = TAU(:,:,2) - TAU(:,:,1);
    DELTA_LWP = (LW_PATH(:,:,2)./LW_PATH(:,:,1) - 1) *100; DIFF_TAU = LW_PATH(:,:,2) - LW_PATH(:,:,1);
    
    
    figure(1)
    ax1 = subplot(2, 3, 3*i-2)
    [C, h] = contourf(ax1, N1_vec, N2_vec, DELTA_CDNC', ncont); hold on
    
    ax2 = axes('Position',get(ax1, 'Position'));
    if i == 1
        scatter(ax2,MH_MODE_CONC(:,1), MH_MODE_CONC(:,2), NI_MKSZ, DELTA_CDNC_MH,'LineWidth',MK_EDGE_WDTH); hold on
        scat_NUM = scatter(ax2, 2000, 30, CHAR_MKSZ,'filled','LineWidth',MK_EDGE_WDTH,'MarkerEdgeColor',[1 1 1], 'MarkerFaceColor',COLORS(3,:)); hold on
        scat_MAV = scatter(ax2,226,135,CHAR_MKSZ, 'filled','LineWidth',MK_EDGE_WDTH,'MarkerEdgeColor',[1 1 1], 'MarkerFaceColor',COLORS(1,:)); hold on
        rec_marine=patch(ax2,[150 600 600 150],[60, 60, 250, 250],...
            'b','FaceAlpha',0,'Edgecolor','c','Linewidth',BOX_LIN_WDTH); hold on
        rec_arctic=patch([40 180 180 40],[20, 20, 150, 150],'b','FaceAlpha',0,'Edgecolor',[0 0 1],'Linewidth',BOX_LIN_WDTH)
        colormap(ax1, flipud(cmocean(COL_SCALE{i}, ncont)));
        colormap(ax2, flipud(gray(ncont)));
        title(ax1, strcat(ALPHABET(ii),')')); ii = ii + 1;
    elseif i == 2
        scat_HYY = scatter(ax2,1110,540,CHAR_MKSZ, 'filled','LineWidth',MK_EDGE_WDTH,'MarkerEdgeColor',[1 1 1], 'MarkerFaceColor',COLORS(2,:)); hold on
        rec_boreal=patch(ax2,[440 1780 1780 440],[160, 160, 920, 920],...
            'b','FaceAlpha',0,'Edgecolor','g','Linewidth',BOX_LIN_WDTH); hold on
        colormap(ax1, flipud(cmocean(COL_SCALE{i}, ncont)));
        title(ax1, strcat(ALPHABET(ii),')')); ii = ii + 1;
    elseif i == 3
        colormap(ax1, flipud(cmocean(COL_SCALE{i}, ncont)));
    end
    [C_10,h_10] = contour(ax2, N1_vec, N2_vec, DELTA_CDNC', [10 10])
    set(h_10, 'Linewidth', CONT_WIDTH, 'Linestyle', ':', 'Color', COLSENV{i})
    [C_A5,h_A5] = contour(ax2, N1_vec, N2_vec, DELTA_ALPHA', [5 5])
    set(h_A5, 'Linewidth', CONT_WIDTH, 'Linestyle', '--', 'Color', COLSENV{i})
        
    %Save out contours for fitting
    save(strcat(ENV{i},'_contours.mat'),'C_10','C_A5')
    
    %link axes and make top axio disappeario
    set([ax1 ax2],'Yscale','log','Xscale','log','XLim',[1 10000],'YLim',[10 1000]);
%    linkaxes([ax1 ax2]);
    ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = [];
    clear COLBAR_XY
    AX_POS  = get(ax1, 'Position');
    COLBAR_XY = AX_POS.*[1 1 0 0] + [0.005 COLBAR_Y0 0.013 COLBAR_LEN];
    COLBAR = colorbar(ax1,'Position',COLBAR_XY); ylabel(COLBAR, '\Delta_{CDNC} [%]')
     linkaxes([ax1 ax2]);
   
    if i == 1
        COLBAR2 = colorbar(ax2,'Position',COLBAR_XY+[0.018 0 0 0]);
        set(COLBAR2,'Color','White','AxisLocation','in','FontWeight','bold','TickDir','out')
        caxis(ax2, CLIM(i,:))
    end
    caxis(ax1, CLIM(i,:));
    
    set(COLBAR,'Color','White','AxisLocation','in','FontWeight','bold')
    YTICKCB = CLIM(i,1) : diff(CLIM(i,:))/ncont : CLIM(i,2);
    set(COLBAR, 'YTick', YTICKCB)
    set(ax1,'XTick',XTick,'XTickLabel',XTicklabs,'YTick',YTick,'YTickLabel',YTicklabs,'TickDir','out');
    set(h, 'linestyle', 'none')
    xlabel(ax1, 'Aitken mode concentration [cm^{-3}]'); ylabel(ax1, 'Accumulation mode concentration [cm^{-3}]');
    
    
    if i == 1
        set(COLBAR, 'YTick',[] , 'YTickLabel',[]); ylabel(COLBAR, ''); ylabel(COLBAR2, '\Delta_{CDNC} [%]')
        YTICKCB = CLIM(i,1) : diff(CLIM(i,:))/ncont : CLIM(i,2)
        set(COLBAR2, 'YTick', YTICKCB)
    end
        
    subplot(2, 3, 3*i - 1)
    [C, h] = contourf(N1_vec, N2_vec, DELTA_TAU', ncont); hold on
    set(gca, 'Xscale','log', 'Yscale','log','XLim',[1 10000],'YLim',[10 1000]); set(h, 'linestyle', 'none')
    xlabel('Aitken mode concentration [cm^{-3}]');
    ylabel('Accumulation mode concentration [cm^{-3}]');
    title(strcat(ALPHABET(ii),')')); ii = ii + 1;
    [C_10_COA,h_10] = contour( N1_vec, N2_vec, DELTA_CDNC', [10 10])
    set(h_10, 'Linewidth', CONT_WIDTH, 'Linestyle', ':', 'Color', COLSENV{i})
    [C_A5_COA,h_A5] = contour(N1_vec, N2_vec, DELTA_ALPHA', [5 5])
    set(h_A5, 'Linewidth', CONT_WIDTH, 'Linestyle', '--', 'Color', COLSENV{i})
    
    if i == 1
        colormap(gca, flipud(cmocean(COL_SCALE{i}, ncont)));
    elseif i == 2
        colormap(gca, flipud(cmocean(COL_SCALE{i}, ncont)));
    end
    caxis(gca,CLIM_TAU(i,:))
    
    AX_POS  = get(gca, 'Position');
    COLBAR_XY = AX_POS.*[1 1 0 0] + [0.005 COLBAR_Y0 0.013 COLBAR_LEN];
    COLBAR = colorbar(gca,'Position',COLBAR_XY); ylabel(COLBAR, '\Delta_{\tau} [%]')
    set(COLBAR,'Color','White','AxisLocation','in','FontWeight','bold')
    YTICKCB = CLIM_TAU(i,1) : diff(CLIM_TAU(i,:))/ncont : CLIM_TAU(i,2)
    set(COLBAR, 'YTick', YTICKCB,'TickDir','out')
    set(gca,'XTick',XTick,'XTickLabel',XTicklabs,'YTick',YTick,'YTickLabel',YTicklabs,'TickDir','out')
    
    
    if i == 1
        snTAU   = scatter(gca, 2000, 30, CHAR_MKSZ,'filled','LineWidth',...
                          MK_EDGE_WDTH,'MarkerEdgeColor',[1 1 1],...
                          'MarkerFaceColor',COLORS(3,:)); hold on
        smTAU   = scatter(gca, 226,135,CHAR_MKSZ, 'filled','LineWidth',...
                          MK_EDGE_WDTH,'MarkerEdgeColor',[1 1 1],...
                          'MarkerFaceColor',COLORS(1,:)); 
        rmTAU   = patch(gca, [150 600 600 150],[60, 60, 250, 250],...
                        'b','FaceAlpha',0,'Edgecolor','c',...
                        'Linewidth',BOX_LIN_WDTH);
        raTAU   = patch(gca,[40 180 180 40],[20, 20, 150, 150],...
                        'b','FaceAlpha',0,'Edgecolor',[0 0 1],...
                        'Linewidth',BOX_LIN_WDTH);
        LMAR      = legend([rmTAU, smTAU, raTAU, snTAU],'Marine','MA','Arctic','NE', 'NumColumns', 2, 'Location','Northwest')
%LMAR      = legend([rmTAU, raTAU], 'marine', 'Arctic marine', 'NumColumns', 2, 'Location','Northwest');

        set(LMAR,'Color','None','TextColor','White','FontWeight','bold'...
            ,'EdgeColor','None','Fontsize',12)
    else
        sbTAU   = scatter(gca, 1110, 540, CHAR_MKSZ, 'filled','LineWidth',...
                           MK_EDGE_WDTH,'MarkerEdgeColor',[1 1 1], ...
                           'MarkerFaceColor',COLORS(2,:)); hold on
        rbTAU = patch(gca, [440 1780 1780 440],[160, 160, 920, 920],...
                           'b','FaceAlpha',0,'Edgecolor','g','Linewidth',...
                           BOX_LIN_WDTH);
        LHYY      = legend([rbTAU, sbTAU], 'Boreal', 'HYY', 'Location','Northwest');
        
        set(LHYY,'Color','None','TextColor','White','FontWeight','bold'...
        ,'EdgeColor','None','Fontsize',12)
    end
    
    subplot(2, 3, 3*i)
    [C, h] = contourf(N1_vec, N2_vec, DELTA_ALPHA', ncont); hold on
    set(gca, 'Xscale','log', 'Yscale','log','XLim',[1 10000],'YLim',[10 1000]); set(h, 'linestyle', 'none')
    xlabel('Aitken mode concentration [cm^{-3}]');
    ylabel('Accumulation mode concentration [cm^{-3}]');
    title(strcat(ALPHABET(ii),')')); ii = ii + 1;
    [C_10_COA,h_10] = contour( N1_vec, N2_vec, DELTA_CDNC', [10 10])
    set(h_10, 'Linewidth', CONT_WIDTH, 'Linestyle', ':', 'Color', COLSENV{i})
    [C_A5_COA_SS,h_A5_SS] = contour(N1_vec, N2_vec,DELTA_CDNC', [ALBEDO5_DELTACDNC(i) ALBEDO5_DELTACDNC(i)])
    set(h_A5_SS, 'Linewidth', CONT_WIDTH, 'Linestyle', '-', 'Color', COLSENV{i})
    [C_A5_COA,h_A5] = contour(N1_vec, N2_vec, DELTA_ALPHA', [5 5])
    set(h_A5, 'Linewidth', CONT_WIDTH, 'Linestyle', '--', 'Color', COLSENV{i})
    
    if i == 1
        snALPHA   = scatter(gca, 2000, 30, CHAR_MKSZ,'filled','LineWidth',...
                           MK_EDGE_WDTH,'MarkerEdgeColor',[1 1 1],...
                           'MarkerFaceColor',COLORS(3,:));
        smALPHA   = scatter(gca, 226,135,CHAR_MKSZ, 'filled','LineWidth',...
                           MK_EDGE_WDTH,'MarkerEdgeColor',[1 1 1],...
                           'MarkerFaceColor',COLORS(1,:)); 
        rmALPHA = patch(gca, [150 600 600 150],[60, 60, 250, 250],...
                           'b','FaceAlpha',0,'Edgecolor','c',...
                           'Linewidth',BOX_LIN_WDTH);
        raALPHA   = patch(gca,[40 180 180 40],[20, 20, 150, 150],...
                        'b','FaceAlpha',0,'Edgecolor',[0 0 1],...
                        'Linewidth',BOX_LIN_WDTH);   
        LMAR = legend([h_10, h_A5, h_A5_SS],'\Delta_{CDNC} = 10 %' ,'\Delta_{\alpha} = 5 %',...
                '\Delta_{\alpha} = 5 % (SS approx.)','Location','Northwest',...
                'Color','none','TextColor','White','FontWeight','bold',...
                'EdgeColor','None','Fontsize',12)           
    else
        shALPHA   = scatter(gca, 1110, 540, CHAR_MKSZ, 'filled','LineWidth',...
                           MK_EDGE_WDTH,'MarkerEdgeColor',[1 1 1], ...
                           'MarkerFaceColor',COLORS(2,:));
        rbALPHA = patch(gca, [440 1780 1780 440],[160, 160, 920, 920],...
                           'b','FaceAlpha',0,'Edgecolor','g','Linewidth',...
                           BOX_LIN_WDTH);
        LHYY      = legend([h_10, h_A5, h_A5_SS],'\Delta_{CDNC} = 10 %' ,'\Delta_{\alpha} = 5 %',...
            '\Delta_{\alpha} = 5 % (SS approx.)','Location','Northwest',...
            'Color','none','TextColor','White','FontWeight','bold',...
            'EdgeColor','None','Fontsize',12) 
    end
    
    if i == 1
        colormap(gca, flipud(cmocean(COL_SCALE{i}, ncont)));
    elseif i == 2
        colormap(gca, flipud(cmocean(COL_SCALE{i}, ncont)));
    end
    caxis(gca,CLIM_ALPHA(i,:))
    AX_POS  = get(gca, 'Position');
    COLBAR_XY = AX_POS.*[1 1 0 0] + [0.005 COLBAR_Y0 0.013 COLBAR_LEN];
    COLBAR = colorbar(gca,'Position',COLBAR_XY); ylabel(COLBAR, '\Delta_{\alpha} [%]')
    set(COLBAR,'Color','White','AxisLocation','in','FontWeight','bold')
    YTICKCB = CLIM_ALPHA(i,1) : diff(CLIM_ALPHA(i,:))/ncont : CLIM_ALPHA(i,2)
    set(COLBAR, 'YTick', YTICKCB,'TickDir','out')
    set(gca,'XTick',XTick,'XTickLabel',XTicklabs,'YTick',YTick,'YTickLabel',YTicklabs,'TickDir','out')
    
end


fig = gcf;
set(findall(fig,'-property','FontSize'),'FontSize',7)
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8 4];
print('Figure4','-dpng','-r400')
