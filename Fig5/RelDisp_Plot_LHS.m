%/----------------------------------X------------------------------------/%
%Sam Lowe, ACES (SU) 2019: Plotting Figure 5, Lowe et al. 2019, Nat. Comms
%Latine hypercube model dispersion
%/----------------------------------X------------------------------------/%
close all; clear all
%function [] = RelDisp_Plot_LHS(pindex, titstr)

PAR.MIN(1, :)    = [150,  15.5, 1.4,  60,  70,  1.4,  0.05, 0.05, 0.05, 0.16, 30];
PAR.MAX(1, :)    = [600,  23.5, 1.8,  250, 100, 1.8,  0.40, 0.40, 2.40, 0.30, 50];
PAR.MIN(2, :)    = [450,  16,   1.55, 215, 75,  1.4,  0.60, 0.60, 0.05, 0.16, 30];
PAR.MAX(2, :)    = [1600, 32,   1.9,  690, 105, 1.75, 0.71, 0.71, 2.40, 0.30, 50];
NPARS            = size(PAR.MAX,2);
PAR_HS_RANGE     = zeros(2,NPARS); PAR_RANGE        = zeros(2,NPARS);
CDNC_THRESHOLDS  = [25, 28, 39];   %CDNC % change for albedo ehancement 5%
NSAMP            = 500;
OUTPUT_STRING{1} = '5000_Marine_OPTICS_LHSOUTPUT.mat';
OUTPUT_STRING{2} = 'REV_HYY_5000_Boreal_OPTICS_LHSOUTPUT.mat';
titles{1}        = '\rm a) Marine';  titles{2}        = '\rm b) Boreal';
UNIT_CONV        = [1e6, 1e-3, 1, 1e6, 1e-3, 1, 1, 1, 1, 1e-9, 1e-3 ];
UNIT_CONV        = repmat(UNIT_CONV, NSAMP, 1);
BARWIDTH         = 1.5;
LINW             = 2.5;
NSAMP            = 100;
HIST_NBIN        = 15;
STYLE            = 'stairs';
cpsz             = [20 20];
ALPHABET         = char(97:122); MKER = 'os';
bwidth           = 2;
alpha            = 0.55;        %Back ground transparency
EPS_SURF         = [0.7, 0.3]; %Global fraction surface coverage (marine, continental)
EPS_CLOUD        = [0.7, 0.5]; %Cloud fraction over surface (marine, continental)
XTickLabsFull    = {'     Maximum * supersaturation';...
                    'Smallest activated *         dry size';...
                    '       Cloud droplet * number concentration';...
                    'Liquid water *       path';...
                    'Cloud optical *    thickness';
                    'Cloud-top *   albedo'; '         Shortwave cloud * radiative effect differential'};
PRCS             = 95;         %Percentile for spread metric (and 100 - PRCS)
%Function for replacing "*" with newline for axis labels (XTickLabsFull)
labels           = cellfun(@(x) strrep(x,'*','\newline'), XTickLabsFull,'UniformOutput',false);
COLORS           = [0 0 1; 0 0.67 0.10]; 
COLORS           = [0 0 1; 0 0.77 0.15]; 
YLIMS            = [0.15 2.5]; 
XLIMS            = [0 35];


figure(1); %subplot(5,1,pindex)
patch([-20 XLIMS(2) XLIMS(2) -20], [YLIMS(1) YLIMS(1) YLIMS(2) YLIMS(2)],...
      [-20 XLIMS(2) XLIMS(2) -20], 'FaceColor','interp','EdgeColor','none',...
      'FaceAlpha',alpha); hold on; xlim(XLIMS); ylim(YLIMS)
      colormap(gca,cmocean('balance','Pivot',18))
l1 = line(XLIMS, [1 1]); set(l1,'Linewidth',1.2,'Color',[1 1 1]*0.5, 'Linestyle','-')

t = text(0.8,1.705,'Microphysical','FontSize',14,'FontWeight','bold');
[Xf, Yf] = ds2nfu(gca, [9 19], [1.7125 1.7125]); 
annotation('arrow', Xf,Yf)

t = text(19.5,1.705,'Cloud','FontSize',14,'FontWeight','bold');
[Xf, Yf] = ds2nfu(gca, [23 30.5], [1.7125 1.7125]); 
annotation('arrow', Xf,Yf)

t = text(31,1.705,'Global','FontSize',14,'FontWeight','bold');

t = text(15.5,1.7,['Columnar' newline 'integration'],'FontSize',12,...
        'FontWeight','bold','HorizontalAlignment','Center','Color',[0 0.45 0.74]);
t = text(25.75,1.7,['Droplet' newline 'scattering'],'FontSize',12,...
        'FontWeight','bold','HorizontalAlignment','Center','Color',[0.64 0.08 0.18]);

for ienv = 1:2
clear OUTPUT
    %Load LHS output-------------------------------------------------------
    load(OUTPUT_STRING{ienv});
%     PARAMS              = PAR.MIN(ienv, :) + (PAR.MAX(ienv, :) - PAR.MIN(ienv, :)) ...
%                           .* QUOT_ARRAY;
%     PARAMS              = PARAMS ./ UNIT_CONV; %convert units
%     
    %Concatenate outputs, get stats----------------------------------------
    OUTPUT              = [SMAX; RACT_OUT(:,:,1); CDNC; LIQWATPATH; TAU; ALPHA];
    DELTA_OUTPUT        = 100 * (OUTPUT(2:2:end, :) ./ OUTPUT(1:2:end, :) - 1);
    DIFF_OUTPUT         = OUTPUT(2:2:end, :) - OUTPUT(1:2:end, :);
    
    %Determine forcing differential from alpha differential
    DIFF_FORCING        = 1370 / 4 * DIFF_OUTPUT(6,:) * EPS_SURF (ienv) * EPS_CLOUD (ienv);
    OUTPUT              = [OUTPUT; DIFF_FORCING];
    
    %Specify statistical defintions (average, spread)
    OUTPUT_AVG          = mean(OUTPUT, 2);
    OUTPUT_STD          = var(OUTPUT, [], 2).^ (1/2);
    OUTPUT_SPAN         = prctile(OUTPUT,95,2) - prctile(OUTPUT,5,2);
    OUTPUT_DISP         = abs(OUTPUT_SPAN ./ OUTPUT_AVG);
    
    %Plotting-------------------------------------------------------------
    YPLOT         = OUTPUT(1:13,:)' ./ OUTPUT_AVG(1:13)';
    ERRORB(:,1:2) = [prctile(OUTPUT,PRCS,2), prctile(OUTPUT,100-PRCS,2)];
    x0            = 1;  
    x0store(1)    = 1;
    ymv           = ERRORB(:,1)./OUTPUT_AVG(1:13);
    ypv           = ERRORB(:,2)./OUTPUT_AVG(1:13);
    
    if ienv == 2
        bwidth = bwidth/2; x0store = x0store + bwidth/2;
    end
   
    DISPERSION(:,:,ienv) = [ymv, ypv];

    
    for ii = 1:length(ERRORB)
        ym         = ymv(ii);
        yp         = ypv(ii);
        if ienv == 1
            if mod(ii,2) == 0
                xm = x0; xp = x0 + bwidth; x0 = xp + 1;
                linsty = ':';
            else
                xm = x0; xp = x0 + bwidth; x0 = xp;
                linsty = '-';
            end
            x0store(ii+1) = x0;
        elseif ienv == 2
            xm = x0store(ii); xp = x0store(ii) + bwidth;
            if mod(ii,2) == 0
                linsty = ':';
            else
                linsty = '-';
            end
        end
        
        if ii == 13
            linsty = '--';
        end
      
        
        p(ienv, ii) = patch([xm xp xp xm],[ym ym yp yp],'r'); hold on; box on
        set(p(ienv, ii), 'FaceColor', 'none', 'EdgeColor', COLORS(ienv,:),...
            'Linestyle', linsty, 'Linewidth',2);
        set(gca, 'YGrid', 'off', 'XGrid', 'on')
        %title('Total parametric uncertainties')
    end
    

%ABSOLUTES------------------------------------------------------------
%      figure(2)
%      for ip = 1:6
%         subplot(2, 3, ip)
%         YPLOT = OUTPUT_AVG(2*ip-1:2*ip,:); ERRORB = OUTPUT_SPAN(2*ip-1:2*ip,:)/2;
%         XLABS{1} = XTickLabs{2*ip-1}; XLABS{2} = XTickLabs{2*ip};
%
%         XPLOT = 1:2;
%         if ienv == 2
%             XPLOT = 3:4;
%         end
%
%         for j = 1:2
%         eb_h = errorbar(XPLOT(j),YPLOT(j),ERRORB(j),'Marker',MKER(j),'MarkerSize',10,...
%             'MarkerEdgeColor','k','MarkerFaceColor',COLORS(ienv,:),...
%             'Color',COLORS(ienv,:),'linewidth',1.5,'capsize',cpsz(ienv)); hold on
%         end
%         xlim([0 5]); set(gca,'XTick',[]); ylabel(BOX_YLABS{ip}) %,'XTickLabel',[XLABS XLABS])
%         title(strcat('\rm ',ALPHABET(ip+1),')'))
%      end
%----------------------------------------------------------------------
    
%DISCRETE PROB DENSITY FUNCTIONS---------------------------------------
%     figure(2 + ienv)
%     for io = 1:12
%         subplot(4,6,io)
%         h(io)    = histogram(OUTPUT(io,:), HIST_NBIN,'Normalization', 'pdf',...
%                       'DisplayStyle',STYLE, 'Linewidth', LINW,'EdgeColor','k'); hold on
%
%         h(io+24) = histogram(OUTPUT_HS(io,:), HIST_NBIN,'Normalization', 'pdf',...
%                       'DisplayStyle',STYLE, 'Linewidth', LINW/2,'EdgeColor','r');
%         xlabel(XTickLabs{io})
%     end
%     title(titles{ienv})
%---------------------------------------------------------------------
    
end
save('Figure5_data.mat','DISPERSION')


xtix      = x0store(2:2:end) - bwidth/2; 
xtix(end) = xtix(end) - bwidth;

figure(1)
%set(gca,'XTick',xtix,'XTickLabel', xtix, 'XTickLabelRotation', 45,'Xgrid','on')
set(gca,'XTick',xtix,'XTickLabel', labels, 'XTickLabelRotation', 0,'Xgrid','on',...
    'YTick',0.25:0.25:2.5,'YLim',[0.15 2.5])

ylabel('Mean-normalised parameteric uncertainties');

figure(1)
h = zeros(5, 1);
h(1) = plot(-5,-5,'k','Linewidth',2);
h(2) = plot(-5,-5,':k','Linewidth',2);
h(3) = plot(-5,-5,'--k','Linewidth',2);
h(4) = plot(-5,-5,'color',COLORS(1,:),'Linewidth',2);
h(5) = plot(-5,-5,'color',COLORS(2,:),'Linewidth',2);

legend(h, 'Bulk Köhler (BK)','Approx. compressed film (CF)','CF-BK differential',...
       'Marine','Boreal','Location','South','NumColumns',2)
       legend boxoff;

%end
       
%--------------------------------------------------------------------------      
%CODE GRAVEYARD------------------------------------------------------------

%title(strcat('\rm ',ALPHABET(1),')'));
%leg = legend('Marine BK', 'Marine CF', 'Boreal BK', 'Boreal CF',...
%'NumColumns', 2,'Location','Northeast');
%legend boxoff;

%'\Delta_{CDNC}';    '\Delta_{S_{max}}'; '\Delta_{r^*}';...
%'\delta_{CDNC}';    '\delta_{S_{max}}'; '\delta_{r^*}';...
%'\Delta_{\alpha}';  '\Delta_{\tau}';    '\Delta_{LWP}';...
%'\delta_{\alpha}';  '\delta_{\tau}';    '\delta_{LWP}'};

%XPLOT = [1:2:13]; XPLOT = sort([XPLOT, XPLOT+0.5]);

% patch([0 15 15 0], [0.5 0.5 1.5 1.5],[0 15 15 0],'FaceColor','interp','EdgeColor','none','FaceAlpha',alpha); hold on
% colormap(gca,cmocean('ice')); ax1 = gca
% ax2 = axes('Position',get(ax1,'Position'))
% patch(ax2,[15 30 30 15], [0.5 0.5 1.5 1.5],[15 30 30 15],'FaceColor','interp','EdgeColor','none','FaceAlpha',alpha)
% colormap(ax2,cmocean('amp'))
% linkaxes([ax1, ax2]); ax2.Visible = 'off';

%COLORS       = [0.00, 0.45, 0.74;...
%                0.47, 0.67, 0.19;...
%                0.85, 0.33, 0.1];
%--------------------------------------------------------------------------      
%--------------------------------------------------------------------------      

