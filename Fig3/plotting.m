clear all
close all

ENV_STRINGS{1} = 'Marine';  ENV_STRINGS{2} = 'HYY_REV_Boreal'; ENV_STRINGS{3} = 'NUMevent';
CONT_COL{1}    = 'c'; CONT_COL{2}    = 'g'; CONT_COL{3}    = 'r';
ALBEDO5_DELTACDNC = [25 28 39];

ncont   = 10;mkwdth  = 0.5; box_col = [0 0 0];

box = [0.05 0.1 (0.4 - 0.05) 2.35;...
    0.45 0.1 (0.75 - 0.45) 2.35];

forg = [0.2 0.60];
COLORS       = [0.00, 0.45, 0.74;...
                0.47, 0.67, 0.19;...
                0.85, 0.33, 0.1];
            
            TITLE = {'a) \rm Accumulation mode conc. N_2 = 30 cm^{-3}',...
                     'b) \rm Accumulation mode conc. N_2 = 135 cm^{-3}';...
                     'c) \rm Accumulation mode conc. N_2 = 160 cm^{-3}',...
                     'd) \rm Accumulation mode conc. N_2 = 540 cm^{-3}'}
j = 1;
for ienv = 1:2
    load(strcat(ENV_STRINGS{ienv}, '_FWN2_DELTACDNC_NEW.mat'))
    UPDRAFT = UPDRAFT_vec; FORG = FORG_vec;
    
    for in = 1:2
        CDNC_CF = squeeze(CDNC(in, :, :, 2)); CDNC_BK = squeeze(CDNC(in, :, :, 1));
        SMAX_CF = squeeze(SMAX(in, :, :, 2)); SMAX_BK = squeeze(SMAX(in, :, :, 1));
        
        diffSMAX(in, :, :) = SMAX_BK - SMAX_CF;
        
        DELTA_CDNC = (CDNC_CF./CDNC_BK - 1) * 100;
        
        RACT_BK = squeeze(RACT_OUT(in, :, :, 1, 1)); RACT_CF = squeeze(RACT_OUT(in, :, :, 2, 1));
        
        dcdnc_max(ienv,in) = max(max(DELTA_CDNC))
        
        %         figure(5+in)
        %         subplot(2,2,1)
        %             contourf( FORG, UPDRAFT, CDNC_BK, ncont,'Linestyle','none'); hold on
        %             set(gca,'Yscale','log'); colormap(gca,cmocean('matter',10)); colorbar
        %         subplot(2,2,2)
        %             contourf( FORG, UPDRAFT, CDNC_CF, ncont,'Linestyle','none'); hold on
        %             set(gca,'Yscale','log'); colormap(gca,cmocean('matter',10)); colorbar
        %         subplot(2,2,3)
        %             contourf( FORG, UPDRAFT, RACT_BK, ncont,'Linestyle','none'); hold on
        %             set(gca,'Yscale','log'); colormap(gca,cmocean('matter',10)); colorbar
        %         subplot(2,2,4)
        %             contourf( FORG, UPDRAFT, RACT_CF, ncont,'Linestyle','none'); hold on
        %             set(gca,'Yscale','log'); colormap(gca,cmocean('matter',10)); colorbar
        
        
        figure(1)
        ax1 = subplot(2,2,in + 2*(ienv - 1))
        [C h] = contourf(ax1, FORG, UPDRAFT, DELTA_CDNC, ncont); hold on
        [C10, h10] = contour(ax1, FORG, UPDRAFT, DELTA_CDNC, [10 10]);
        [CA5, hA5] = contour(ax1, FORG, UPDRAFT, DELTA_CDNC, ...
            [ALBEDO5_DELTACDNC(ienv) ALBEDO5_DELTACDNC(ienv)]);
        set(hA5,'Linewidth',2.5,'Linestyle','--','Color',CONT_COL{ienv})
        set(h10,'Linewidth',2.5,'Linestyle', ':','Color',CONT_COL{ienv})
        title(TITLE{ienv,in})
        set(gca,'YTickLabel',[0.1, 1, 10])
        if in == 2
            if ienv < 3
                rectangle('Position',box(ienv,:), 'EdgeColor',box_col,'Linewidth',2,'linestyle',':');
                scatter(ax1,forg(ienv), 0.3192,75,COLORS(ienv,:),'o','filled','LineWidth',mkwdth,'MarkerEdgeColor',[1 1 1])
            end
        end
        
        set(ax1, 'Yscale','log'); COLBAR = colorbar;
        set(h, 'Linestyle','none')
        ylabel(COLBAR,'\Delta_{CDNC} [%]')
        if ienv == 1
            colormap(gca , cmocean('dense',ncont))
            %set(ax1,'colorscale','log')
            caxis(ax1, [0 90])
                 YTICKCB = 0 : 2 * 90/ncont : 90
    set(COLBAR, 'YTick', YTICKCB)
        elseif ienv == 2
            colormap(gca ,cmocean('tempo',ncont))
            %set(ax1,'colorscale','log')
            caxis(ax1, [0 65])
                 YTICKCB = 0 :2 *  65/ncont : 65
    set(COLBAR, 'YTick', YTICKCB)
        elseif ienv == 3
            colormap(ax1 ,cmocean('amp',ncont))
            %caxis(ax1, [-80 50])
        end

    
        xlabel('Organic mass fraction'); ylabel('Updraft [ms^{-1}]')
        
        
    end
end

fig = gcf;
set(findall(fig,'-property','FontSize'),'FontSize',9)
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8 5];
print('Figure3','-dpng','-r400')


return



ncont = 30
for ienv = 1:3
    load(strcat(ENV_STRINGS{ienv}, '_FWN2_DELTACDNC.mat'))
    UPDRAFT = UPDRAFT_vec; FORG = FORG_vec;
    
    for in = 1:2
        
        RACT_BK  = squeeze(RACT_OUT(in, :, :, 1, 1));
        RACT_CF  = squeeze(RACT_OUT(in, :, :, 2, 1));
        diffRACT = RACT_BK - RACT_CF;      ;
        %diffRACT(diffRACT< - 20) = NaN;
        
        
        figure(2)
        ax1 = subplot(3,2,in + 2*(ienv - 1))
        
        [C h] = contourf(ax1, FORG, UPDRAFT, diffRACT, ncont); hold on
        
        %         [C10, h10] = contour(ax1, FORG, UPDRAFT, DELTA_CDNC, [10 10])
        %         set(h10, 'Linewidth', 2.5, 'Linestyle', ':', 'Color', CONT_COL{ienv})
        %         [CA5, hA5] = contour(ax1, FORG, UPDRAFT, DELTA_CDNC, [ALBEDO5_DELTACDNC(ienv) ALBEDO5_DELTACDNC(ienv)])
        %         set(hA5,'Linewidth',2.5,'Linestyle','--','Color', CONT_COL{ienv})
        
        set(ax1, 'Yscale','log'); COLBAR = colorbar; %caxis([0 inf])
        set(h, 'Linestyle','none')
        
        clim =  [min(min(diffRACT)) max(max(diffRACT))]
        if ienv == 1
            colormap(ax1 ,cmocean('ice',ncont))
            %caxis(ax1, clim)
        elseif ienv == 2
            colormap(ax1 ,flipud(cmocean('algae',ncont)))
            set(ax1,'colorscale','log')
            %caxis(ax1, clim)
        elseif ienv == 3
            colormap(ax1 ,flipud(cmocean('amp',ncont)))
            set(ax1,'colorscale','log')
            %caxis(ax1, clim)
        end
        
        xlabel('Organic mass fraction'); ylabel('Updraft [ms^{-1}]')
        
        %    figure(j+1)
        %    s = surf(FORG, UPDRAFT, DELTA_CDNC); set(gca, 'Yscale','log');
        %    set(s,'Linestyle','none')
        %         if ienv == 1
        %             colormap(gca ,cmocean('ice',ncont))
        %         elseif ienv == 2
        %             colormap(gca ,cmocean('algae',ncont))
        %         end
        %      j = j +1
    end
end

