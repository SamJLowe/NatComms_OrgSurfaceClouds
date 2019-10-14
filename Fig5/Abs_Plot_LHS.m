%/----------------------------------X------------------------------------/%
%Sam Lowe, ACES (SU) 2019: Plotting Figure 5, Lowe et al. 2019, Nat. Comms
%Latine hypercube model dispersion
%/----------------------------------X------------------------------------/%
close all; clear all

PAR.MIN(1, :)    = [150,  15.5, 1.4,  60,  70,  1.4,  0.05, 0.05, 0.05, 0.16, 30];
PAR.MAX(1, :)    = [600,  23.5, 1.8,  250, 100, 1.8,  0.40, 0.40, 2.40, 0.30, 50];
PAR.MIN(2, :)    = [450,  16,   1.55, 215, 75,  1.4,  0.60, 0.60, 0.05, 0.16, 30];
PAR.MAX(2, :)    = [1600, 32,   1.9,  690, 105, 1.75, 0.71, 0.71, 2.40, 0.30, 50];
NPARS            = size(PAR.MAX,2);
PAR_HS_RANGE     = zeros(2,NPARS); PAR_RANGE        = zeros(2,NPARS);
CDNC_THRESHOLDS  = [25, 28, 39];   %CDNC % change for albedo ehancement 5%
NSAMP            = 5000;
OUTPUT_STRING{1} = '5000_Marine_OPTICS_LHSOUTPUT.mat';
OUTPUT_STRING{2} = 'REV_HYY_5000_Boreal_OPTICS_LHSOUTPUT.mat';
titles{1}        = 'Marine MA';  titles{2}        = 'Boreal HYY';
UNIT_CONV        = [1e6, 1e-3, 1, 1e6, 1e-3, 1, 1, 1, 1, 1e-9, 1e-3 ];
UNIT_CONV        = repmat(UNIT_CONV, NSAMP, 1);
BARWIDTH         = 1.5;
LINW             = 2.5;
HIST_NBIN        = 10;
STYLE            = 'stairs';
cpsz             = [20 20];
ALPHABET         = char(97:122); MKER = 'os';
bwidth           = 2;
alpha            = 0.55;       %Back ground transparency
EPS_SURF         = [0.7, 0.3]; %Global fraction surface coverage (marine, continental)
EPS_CLOUD        = [0.7, 0.5]; %Cloud fraction over surface (marine, continental)
PRCS             = 95;         %Percentile for spread metric (and 100 - PRCS)
LINSTY{1}        = 'solid';
LINSTY{2}        = 'dotted';
YLabelFull       = {['Maximum supersaturation' newline 'S_{max} [%]'],...
    ['Smallest activated' newline 'dry size r^* [nm]'],...
    ['Cloud droplet number' newline 'concentration CDNC [cm^{-3}]'],...
    ['Liquid water path' newline 'LWP [gm^{-2}]'],...
    ['Cloud optical' newline 'thickness \tau'],...
    ['Cloud Albedo \alpha'],...
    ['Shortwave cloud radiative' newline 'effect differential [Wm^{-2}]']};

%Function for replacing "*" with newline for axis labels (XTickLabsFull)
labels           = cellfun(@(x) strrep(x,'*','\newline'), YLabelFull,'UniformOutput',false);
COLORS           = [0 0 1; 0 0.67 0.10];
%COLORS           = [0 0 1; 0 0.77 0.15];
YLIMS            = [0.2 2];
XLIMS            = [0 35];
bcol = 'cg'

for ienv = 1:2
    %Load LHS output-------------------------------------------------------
    load(OUTPUT_STRING{ienv});
    PARAMS              = PAR.MIN(ienv, :) + (PAR.MAX(ienv, :) - PAR.MIN(ienv, :)) ...
        .* QUOT_ARRAY;
    PARAMS              = PARAMS ./ UNIT_CONV; %convert units
    
    %Concatenate outputs, get stats----------------------------------------
    OUTPUT              = [SMAX; RACT_OUT(:,:,1); CDNC; LIQWATPATH; TAU; ALPHA];
    DELTA_OUTPUT        = 100 * (OUTPUT(2:2:end, :) ./ OUTPUT(1:2:end, :) - 1);
    DIFF_OUTPUT         = OUTPUT(2:2:end, :) - OUTPUT(1:2:end, :);
    
    %Determine forcing differential from alpha differential
    DIFF_FORCING        = - 1370 / 4 * DIFF_OUTPUT(6,:) * EPS_SURF (ienv) * EPS_CLOUD (ienv);
    OUTPUT              = [OUTPUT; DIFF_FORCING];
    
    %Specify statistical defintions (average, spread)
    OUTPUT_AVG          = mean(OUTPUT, 2);
    OUTPUT_STD          = var(OUTPUT, [], 2).^ (1/2);
    OUTPUT_SPAN         = prctile(OUTPUT,95,2) - prctile(OUTPUT,5,2);
   % EBP                 = prctile(OUTPUT,95,2) - OUTPUT_AVG;
   % EBM                 = OUTPUT_AVG - prctile(OUTPUT,5,2);
    OUTPUT_DISP         = abs(OUTPUT_SPAN ./ OUTPUT_AVG);
    
    %ABSOLUTES------------------------------------------------------------
    figure(2)
    sigfig = [2, 2, 2, 3, 2 ,3, 2]
    for ip = 1:6
        
        subplot(3, 3, ip)
        YPLOT  = OUTPUT_AVG(2*ip-1:2*ip);
        ERRORB = OUTPUT_SPAN(2*ip-1:2*ip)/2;
        XPLOT  = [1:2] + 2*(ienv-1);
        
        for j = 1:2 %Kohler scheme loop
%             eb_h = errorbar(XPLOT(j),YPLOT(j),ERRORB(j),'o','MarkerSize',10,...
%                 'MarkerEdgeColor','k','MarkerFaceColor',COLORS(ienv,:),...
%                 'Color',COLORS(ienv,:),'linewidth',2,'capsize',cpsz(ienv)); hold on
        EBM    = OUTPUT_AVG(2*ip+j-2) - prctile(OUTPUT(2*ip+j-2,:),5);
        EBP    = prctile(OUTPUT(2*ip+j-2,:),95) - OUTPUT_AVG(2*ip+j-2);
                eb_h = errorbar(XPLOT(j),YPLOT(j),EBM,EBP,'o','MarkerSize',10,...
                'MarkerEdgeColor','k','MarkerFaceColor',COLORS(ienv,:),...
                'Color',COLORS(ienv,:),'linewidth',2,'capsize',cpsz(ienv)); hold on
            eb_h.Bar.LineStyle = LINSTY{j};
        end
        %set(gca,'Ylim',round(get(gca,'YLim'), 3, 'significant'))
        yl = get(gca,'YLim');
        ytix = round(yl(1) : diff(yl) / 5 : yl(2), sigfig(ip), 'significant');
        set(gca,'YTick',ytix,'YTickLabel',ytix);
        
        ylabel(YLabelFull{ip}); xlim([0 5]);
        title(strcat('\rm ',ALPHABET(ip),')'));
        
    end
    sigfig = 2;
    subplot(3, 3, ip+1)
    XPLOT = [1.5 2];
    EBM    = OUTPUT_AVG(end) - prctile(OUTPUT(end,:),5);
    EBP    = prctile(OUTPUT(end,:),95) - OUTPUT_AVG(end);
%     eb_h = errorbar(XPLOT(ienv),OUTPUT_AVG(end),OUTPUT_SPAN(end)/2,'o','MarkerSize',10,...
%         'MarkerEdgeColor','k','MarkerFaceColor',COLORS(ienv,:),...
%         'Color',COLORS(ienv,:),'linewidth',2,'capsize',cpsz(ienv)); hold on
    eb_h = errorbar(XPLOT(ienv),OUTPUT_AVG(end),EBM,EBP,'o','MarkerSize',10,...
        'MarkerEdgeColor','k','MarkerFaceColor',COLORS(ienv,:),...
        'Color',COLORS(ienv,:),'linewidth',2,'capsize',cpsz(ienv)); hold on

eb_h.Bar.LineStyle = 'dashed';
    set(gca,'Ylim', [-4 0],'Xlim',[1.0 2.5]);
    title(strcat('\rm ',ALPHABET(ip+2),')'));
    ylabel(YLabelFull{7})
    %set(gca,'Ylim',round(get(gca,'YLim'), 3, 'significant'))
    
    yl = get(gca,'YLim');
    ytix = round(yl(1) : diff(yl) / 5 : yl(2), sigfig, 'significant');
    set(gca,'YTick',ytix,'YTickLabel',ytix);
    
    cloudfrac = linspace(0.1,0.7,100)
    avg = OUTPUT_AVG(end)*cloudfrac/EPS_CLOUD(ienv)
    ulim = prctile(OUTPUT(end,:),95)*cloudfrac/EPS_CLOUD(ienv)
    dlim = prctile(OUTPUT(end,:),5)*cloudfrac/EPS_CLOUD(ienv)

    figure(5)
    subplot(1,2,ienv)
    p1=plot(cloudfrac,avg,'color',COLORS(ienv,:),'Linewidth',2);hold on
    p2=plot(cloudfrac,ulim,':','color',COLORS(ienv,:),'Linewidth',1);
    p3=plot(cloudfrac,dlim,':','color',COLORS(ienv,:),'Linewidth',1);
    s1=scatter(EPS_CLOUD(ienv),OUTPUT_AVG(end),65,'MarkerFaceColor',COLORS(ienv,:),'MarkerEdgeColor','k')
    xlabel('Cloud fraction'); ylabel(YLabelFull{7})
    xlim([0.1 0.7])
    legend([p1, p2, p3, s1],'Mean SW CRE','95^{th} percentile SW CRE','5^{th} percentile SW CRE',...
    'Assumed cloud fraction','Location','SouthWest')
     set(findall(gcf,'-property','FontSize'),'FontSize',14)

    %----------------------------------------------------------------------
    
    %DISCRETE PROB DENSITY FUNCTIONS---------------------------------------
    MODEL_OUTPUTS(ienv,:,:) = OUTPUT;
    
    figure(3)
    j=1
    for io = 1:2:14
        subplot(2,4,j)
        h(io)    = histogram(OUTPUT(io,:), HIST_NBIN,'Normalization', 'pdf',...
            'DisplayStyle',STYLE, 'Linewidth', LINW,'EdgeColor',COLORS(ienv,:)); hold on
        if io == 13
            set(h(io), 'Linestyle','--')
            %h7 = h(io)
        elseif j < 7
            h(io+1)    = histogram(OUTPUT(io+1,:), HIST_NBIN,'Normalization', 'pdf',...
                'DisplayStyle',STYLE, 'Linewidth', LINW,'LineStyle',':','EdgeColor',COLORS(ienv,:)); hold on
        end
        %         h(io+24) = histogram(OUTPUT_HS(io,:), HIST_NBIN,'Normalization', 'pdf',...
        %                       'DisplayStyle',STYLE, 'Linewidth', LINW/2,'EdgeColor','r');
        xlabel(YLabelFull{j})
        title(strcat('\rm ',ALPHABET(j),')'));
        j=j+1;
       
    end
    %---------------------------------------------------------------------
    
    
    %figure(4)
    %for ifrac = 1:length(CloudFrac)
        
    %end
    
end

save('FigureS11_Data.mat','MODEL_OUTPUTS')
save('FigureS12_Data.mat','')

hand = zeros(6, 1);
hand(1) = plot(-5,-5,'k','Linewidth',2);
hand(2) = plot(-5,-5,'k:','Linewidth',2);
hand(3) = plot(-5,-5,'color',COLORS(1,:),'Linewidth',2);
hand(4) = plot(-5,-5,'color',COLORS(2,:),'Linewidth',2);
hand(5) = plot(-5,-5,'--','color',COLORS(1,:),'Linewidth',2);
hand(6) = plot(-5,-5,'--','color',COLORS(2,:),'Linewidth',2);
xlim([-6 0]);ylim([0 inf])

subplot(2,4,1)
    legend([hand(1),hand(2),hand(3), hand(4),hand(5),hand(6)],...
        'Bulk Köhler BK','Approx. compressed film CF',...
        'Marine MA','Boreal HYY','MA shortwave cloud \newline radiative effect differential',...
        'HYY shortwave cloud \newline radiative effect differential')
legend boxoff

figure(2)
subplot(3,3,7)
h = zeros(5, 1);
h(1) = plot(-5,-5,'k','Linewidth',2);
h(2) = plot(-5,-5,':k','Linewidth',2);
h(3) = plot(-5,-5,'--k','Linewidth',2);
h(4) = plot(-5,-5,'color',COLORS(1,:),'Linewidth',2);
h(5) = plot(-5,-5,'color',COLORS(2,:),'Linewidth',2);

legend(h, 'Bulk Köhler','Compressed film','Shortwave cloud radiative \newline effect differential',...
    'Marine','Boreal','Location','Southwest','NumColumns',1)
legend boxoff;
set(findall(gcf,'type','axes'), 'XTick',[],'Ygrid','on','FontSize',12)
%sgtitle('All parametric uncertainties','FontSize',14)
return

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

