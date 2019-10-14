close all; clear all

LINE_WDTH   = 1.6;
colors      = [0.00, 0.45, 0.74; 
               0.47, 0.67, 0.19; 
               0.85, 0.33, 0.10; 
               0.86, 0.63, 0.53];
rdry = zeros(3,3);
%Size distribution, panel a)-----------------------------------------------           
XTicklabs   = {'1'; '10'; '100';'1,000'};
YTicklabs   = {'1'; '10'; '100';'1,000'};
XTick       = [1, 10, 100, 1000];
YTick       = [1, 10, 100, 1000];

MAV_AERODRY = load('MAV/Input_Aerosol.txt');
HYY_AERODRY = load('HYY/Input_Aerosol.txt');
NUM_AERODRY = load('NUM/Input_Aerosol.txt');

DRY_RADIUS                = MAV_AERODRY(1, 3:end)*1e9;
[DUMMY, KOHLER_BIN_INDEX] = min(abs(DRY_RADIUS - 50)); %Index of model size class from plotted kohler curve

DNDLOGR       = diff(log(DRY_RADIUS));
DRY_RADIUS(1) = [];

MAV_AERO_CONC = sum(MAV_AERODRY(2:3, 4:end),1)./DNDLOGR;
HYY_AERO_CONC = sum(HYY_AERODRY(2:3, 4:end),1)./DNDLOGR;
NUM_AERO_CONC = sum(NUM_AERODRY(2:3, 4:end),1)./DNDLOGR;

FIG1_DATA.MAV_AERODRY = MAV_AERODRY;
FIG1_DATA.HYY_AERODRY = MAV_AERODRY;
FIG1_DATA.NUM_AERODRY = MAV_AERODRY;
FIG1_DATA.DRY_RADIUS  = DRY_RADIUS;

FIG1 = figure(1); 
subplot(2,2,1)
    semilogx(DRY_RADIUS,MAV_AERO_CONC,'color',colors(1,:),'Linewidth',LINE_WDTH); hold on
    semilogx(DRY_RADIUS,HYY_AERO_CONC,'color',colors(2,:),'Linewidth',LINE_WDTH);
    semilogx(DRY_RADIUS,NUM_AERO_CONC,'color',colors(3,:),'Linewidth',LINE_WDTH);
    set(gca, 'Ylim', [1 2500], 'YTick', YTick, 'YTickLabel', YTicklabs, 'Yscale','log')
    set(gca, 'Xlim', [5 1000], 'XTick', XTick, 'XTickLabel', XTicklabs); 
    xlabel('Dry aerosol radius [nm]'); 
    ylabel({'Aerosol number concentration';'size distribution [cm^{-3}]'})
    dist_leg = legend('Marine (MA)','Boreal (HYY)','NUM event (NE)',...
                      'Location','South'); legend boxoff
title('a)')
drawnow

%COMPOSITION BAR CHART ----------------------------------------------------
NSPECIES = 6; %Number of chemical compounds
colors   = flipud(cmocean('phase', NSPECIES + 1)); 
x        = 1:NSPECIES;    

                     %[NUM1 NUM2 HYY1 HYY2 MA1 MA2]
SOA_DIMER_SURROGATE = [52   0    66.8   0    0   0 ];
SOA_SURROGATE       = [0    30   0      66.8 0   0 ];
PALMITIC            = [0    0    0      0    20  20 ];
AMMONIUM_SULF       = [48   70   16.6   16.6 80  0 ];
AMMONIUM_NITR       = [0    0    16.6   16.6 0   0 ];
SODIUM_CHLORIDE     = [0    0    0      0    0   80];

MASS_FRACS = [SOA_DIMER_SURROGATE; SOA_SURROGATE; PALMITIC;...
              AMMONIUM_SULF; AMMONIUM_NITR; SODIUM_CHLORIDE];
          
MASS_FRACS = horzcat(MASS_FRACS(:,5:6),MASS_FRACS(:,3:4),MASS_FRACS(:,1:2));

FIG1_DATA.MASS_FRACS  = MASS_FRACS;


XBAR = [1,2,4,5,7,8];
subplot(2,2,2)
    bb=bar(XBAR, MASS_FRACS','stacked')
    set(bb,{'FaceColor'},{colors(1,:);colors(2,:);colors(3,:);colors(4,:);...
        colors(5,:);colors(6,:)}, 'Edgecolor','none');
    set(gca, 'XTick', XBAR,...
             'XTickLabel',{'MA AIT', 'MA ACC', 'HYY AIT', 'HYY ACC', 'NE NUF', 'NE ACC'},...
             'Xlim',[0 16], 'Ylim', [0 100],...
             'XTickLabelRotation',45);
    ylabel('Mass fraction [%]')
    MF_leg = legend(['SOA' newline 'surrogate 1'],['SOA' newline 'surrogate 2'],...
                    ['Palmitic' newline 'acid'],['Ammonium' newline 'sulfate'],...
                    ['Ammonium' newline 'nitrate'],['Sodium' newline 'chloride'],...
                     'Location','East'); legend boxoff
    title('b)')
    drawnow

%KOHLER CURVE PLOT, panel c -----------------------------------------------
colors = [0.00, 0.45, 0.74; 
          0.47, 0.67, 0.19; 
          0.85, 0.33, 0.10; 
          0.86, 0.63, 0.53];

XTICK     = [100, 300, 500, 700, 1000, 3000];
XTICK_LAB = {'100'; '300'; '500';'700';'1,000';'3,000'};
      
load('MAV/OUTPUT/kohler_curves.txt'); load('MAV/OUTPUT/kohler_curves_CF.txt')

RADIUS_BK     = kohler_curves(1:3:end-3,:) * 1e9;      %[nm]
EQ_SUSAT_BK   = (kohler_curves(2:3:end-2,:)-1) * 100;  %[%]
SURF_TENS_BK  = kohler_curves(3:3:end-1,:)*1000;       %[mNm-1]
rdry(1,:)  = kohler_curves(end,:);       %[mNm-1]

RADIUS_CF     = kohler_curves_CF(1:3:end-2,:)*1e9;
EQ_SUSAT_CF   = 100*(kohler_curves_CF(2:3:end-1,:)-1);
SURF_TENS_CF  = kohler_curves_CF(3:3:end,:)*1000;
sigw          = 72.8; 
rsigw         = 160:1:3000;

FILM_DISSOLUTION_INDX = max(find(SURF_TENS_CF(:,2) == 40));
R_FILM_DISS(1)        = RADIUS_CF(FILM_DISSOLUTION_INDX);

subplot(2,2,3)
    plot(RADIUS_BK(:,2), EQ_SUSAT_BK(:,2), 'color',colors(1,:),'Linewidth',LINE_WDTH); hold on
    plot(RADIUS_CF(:,2), EQ_SUSAT_CF(:,2), '--','color',colors(1,:),'Linewidth',LINE_WDTH);
    ylim([-0.2 0.35]); xlim([150 3000]); 
    set(gca, 'Xscale', 'log', 'XTick', XTICK, 'XTickLabel', XTICK_LAB)
   xlabel('Wet radius [nm]')
    ylabel({'Equilibrium';'supersaturation [%]'})
    title('c)')

FIG1_DATA.MAV_RADIUS_BK = RADIUS_BK(:,2);
FIG1_DATA.MAV_RADIUS_CF = RADIUS_CF(:,2);
FIG1_DATA.MAV_EQ_SUSAT_BK = EQ_SUSAT_BK(:,2);
FIG1_DATA.MAV_EQ_SUSAT_CF = EQ_SUSAT_CF(:,2);
FIG1_DATA.MAV_ST_CF = SURF_TENS_CF(:,2);

    
subplot(2,2,4)
    plot(rsigw, ones(1,length(rsigw))*sigw,'color',[0.5 0.5 0.5],'Linewidth',1.6); hold on
    plot(RADIUS_CF(:,2),SURF_TENS_CF(:,2),'--','color',colors(1,:),'Linewidth',1.6);
    ylim([10 85]); xlim([200 3000])
    set(gca,'Xscale','log')%,'XTick',[100, 300, 500, 700, 1000, 3000],...
        %'XTickLabeL',{'100'; '300'; '500';'700';'1,000';'3,000'})
    xlabel('Wet radius [nm]'); 
    ylabel('Surface tension [mN m^{-1}]')

    title('d)')

ftsz = 12
ht = text(205, 80 ,'Surface tension of water = 72.8 mNm^{-1}');
ht.Color=[0.5 0.5 0.5];
ht.FontWeight='bold';ht.FontSize=ftsz;
%------
    
load('HYY/OUTPUT/kohler_curves.txt'); load('HYY/OUTPUT/kohler_curves_CF.txt');

RADIUS_BK     = kohler_curves(1:3:end-3,:) * 1e9;      %[nm]
EQ_SUSAT_BK   = (kohler_curves(2:3:end-2,:)-1) * 100;  %[%]
SURF_TENS_BK  = kohler_curves(3:3:end-1,:)*1000;       %[mNm-1]
rdry(2,:)  = kohler_curves(end,:);       %[mNm-1]
R0 = min(RADIUS_BK(:,2))

RADIUS_CF     = kohler_curves_CF(1:3:end-2,:)*1e9;
EQ_SUSAT_CF   = 100*(kohler_curves_CF(2:3:end-1,:)-1);
SURF_TENS_CF  = kohler_curves_CF(3:3:end,:)*1000;

FILM_DISSOLUTION_INDX = max(find(SURF_TENS_CF(:,2) == 40));
R_FILM_DISS(2)        = RADIUS_CF(FILM_DISSOLUTION_INDX);

subplot(2,2,3)
    plot(RADIUS_BK(:,2),EQ_SUSAT_BK(:,2),'color',colors(2,:),'Linewidth',LINE_WDTH); hold on
    plot(RADIUS_CF(:,2),EQ_SUSAT_CF(:,2),'--','color',colors(2,:),'Linewidth',LINE_WDTH);
    ylim([-0.1  0.35])
    
FIG1_DATA.HYY_RADIUS_BK = RADIUS_BK(:,2);
FIG1_DATA.HYY_RADIUS_CF = RADIUS_CF(:,2);
FIG1_DATA.HYY_EQ_SUSAT_BK = EQ_SUSAT_BK(:,2);
FIG1_DATA.HYY_EQ_SUSAT_CF = EQ_SUSAT_CF(:,2);
FIG1_DATA.HYY_ST_CF = SURF_TENS_CF(:,2);

    
subplot(2,2,4)
    plot(RADIUS_CF(:,2),SURF_TENS_CF(:,2),'--','color',colors(2,:),'Linewidth',1.6); hold on
    ylim([10 85]); xlim([100 3000])
    set(gca,'Xscale','log')%,'XTick',[100, 300, 500, 700, 1000, 3000],...
        %'XTickLabeL',{'100'; '300'; '500';'700';'1,000';'3,000'})
    xlabel('Wet radius [nm]'); 
    ylabel('Surface tension [mN m^{-1}]')
    title('d)')
    
%------------
load('NUM/OUTPUT/kohler_curves.txt'); load('NUM/OUTPUT/kohler_curves_CF.txt');

RADIUS_BK     = kohler_curves(1:3:end-3,:) * 1e9;      %[nm]
EQ_SUSAT_BK   = (kohler_curves(2:3:end-2,:)-1) * 100;  %[%]
SURF_TENS_BK  = kohler_curves(3:3:end-1,:)*1000;       %[mNm-1]
rdry(3,:)  = kohler_curves(end,:);       %[mNm-1]


RADIUS_CF     = kohler_curves_CF(1:3:end-2,:)*1e9;
EQ_SUSAT_CF   = 100*(kohler_curves_CF(2:3:end-1,:)-1);
SURF_TENS_CF  = kohler_curves_CF(3:3:end,:)*1000;

FILM_DISSOLUTION_INDX = max(find(SURF_TENS_CF(:,2) == 40));
R_FILM_DISS(3)        = RADIUS_CF(FILM_DISSOLUTION_INDX); 

subplot(2,2,3)
    plot(RADIUS_BK(:,2),EQ_SUSAT_BK(:,2),'color',colors(3,:),'Linewidth',LINE_WDTH); hold on
    plot(RADIUS_CF(:,2),EQ_SUSAT_CF(:,2),'--','color',colors(3,:),'Linewidth',LINE_WDTH);
     xlim([R0 3000])
     
FIG1_DATA.NUM_RADIUS_BK = RADIUS_BK(:,2);
FIG1_DATA.NUM_RADIUS_CF = RADIUS_CF(:,2);
FIG1_DATA.NUM_EQ_SUSAT_BK = EQ_SUSAT_BK(:,2);
FIG1_DATA.NUM_EQ_SUSAT_CF = EQ_SUSAT_CF(:,2);
FIG1_DATA.NUM_ST_CF = SURF_TENS_CF(:,2);
     
h2 = plot(-5,-5,'k--')
h1 = plot(-5,-5,'k')
legend([h1, h2],'Bulk köhler BK','Approx. compressed film CF'); legend boxoff

subplot(2,2,4)
    plot(RADIUS_CF(:,2),SURF_TENS_CF(:,2),'--','color',colors(3,:),'Linewidth',1.6); hold on
    ylim([10 85]); xlim([R0 3000])
    set(gca,'Xscale','log','XTick',[100, 300, 500, 700, 1000, 3000],...
        'XTickLabeL',{'100'; '300'; '500';'700';'1,000';'3,000'})
    xlabel('Wet radius [nm]');
    ylabel('Surface tension [mN m^{-1}]')
    title('d)')    
    
set(findall(FIG1,'-property','FontSize'),'FontSize',8)
set(dist_leg, 'FontSize',7); 
set(MF_leg, 'FontSize',7)

save('Figure1_Data.mat','FIG1_DATA')

return
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8 5];
print('Figure1','-dpng','-r400')

% FIG1.PaperUnits = 'inches';
% FIG1.PaperPosition = [0 0 6.5 5.5];
% print('FIG1','-dpng','-r300')

return

%Rdry = 20.5 nm for review round 2---------------------
figure(2)
subplot(2,1,1)
plot(RADIUS_BK(:,1)*2e-3, EQ_SUSAT_BK(:,1), 'color',colors(3,:),'Linewidth',LINE_WDTH); hold on
plot(RADIUS_CF(:,1)*2e-3, EQ_SUSAT_CF(:,1), '--','color',colors(3,:),'Linewidth',LINE_WDTH);
ylim([-0.2 0.9]); xlim([0.025 0.575]);
%set(gca, 'XTick', XTICK*2e-3, 'XTickLabel', XTICK*2e-3)
xlabel('Wet diameter [\mum]')
ylabel({'Equilibrium';'supersaturation [%]'})
title('a)')    

subplot(2,1,2)
plot(RADIUS_CF(:,1)*2e-3,SURF_TENS_CF(:,1),'--','color',colors(3,:),'Linewidth',1.6); hold on
ylim([10 85]); xlim([0.025 0.575]);
%set(gca, 'XTick', XTICK*2e-3, 'XTickLabel', XTICK*2e-3)
xlabel('Wet diameter [\mum]');
ylabel('Surface tension [mN m^{-1}]')
title('b)')    
