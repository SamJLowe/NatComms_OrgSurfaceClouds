close all, clear all

ENV_STRING{1} = 'Marine';
ENV_STRING{2} = 'Boreal'
ENV_STRING{3} = 'NUM event'
NAT           = [2 2 2];

SURFACE_FLAG  = [0 1];
UPDRAFT       = 0.32;
%UPDRAFT       = linspace(0.05, 2.4, 10);

%Plotting vars
colors      = [0.00, 0.45, 0.74; 
               0.47, 0.67, 0.19; 
               0.85, 0.33, 0.10; 
               0.86, 0.63, 0.53];
           
lsty{1}    = '-'; lsty{2} = '--'; linw = 1.5
YTicklabs = {'0';'1,000'; '2,000'; '3,000'; '4,000'; '5,000'; '6,000'; '7,000'};
CLOUD_TOP = 200;

Labelstr{1} = 'Marine';  
Labelstr{2} = 'Boreal continental';     
Labelstr{3} = 'NUM-event';

labposx     = [130 380 45];  
labposy     = [97 97 97];
for iw = 1:length(UPDRAFT)
    for i = 1:3
        for j = 1:2
            iw
            i
            j
            GenInputFiles(SURFACE_FLAG(j), UPDRAFT(iw), ENV_STRING{i}, CLOUD_TOP);
            [PROFILES, CLOUDTOP, RACT, OPTICS, liqwp] = RUN_CLOUD(ENV_STRING{i}, NAT(i));
            CDNC(i,j,iw)     = CLOUDTOP.CDNC_rc
            SMAX(i,j,iw)     = CLOUDTOP.Smax;
            tau(i,j,iw)      = OPTICS.COD;
            alpha(i,j,iw)    = OPTICS.ALBEDO;
            R_ACT(i,j,iw,:)  = RACT;
            LWP(i,j,iw)      = liqwp;
            ZPROFILES(i,j,:,:)        = PROFILES;

            figure(1)
            subplot(2,2,1)
               plot(PROFILES(:,5), PROFILES(:,1),'color',colors(i,:),...
                   'Linewidth',linw, 'Linestyle',lsty{j}); hold on;
               xlabel('Supersaturation [%]'); ylabel('Parcel displacement [m]'); 
               ylim([0 100]); xlim([0 0.5]); title('a)')

            subplot(2,2,2)
            plot(PROFILES(:,4), PROFILES(:,1),'color',colors(i,:),...
               'Linewidth',linw, 'Linestyle',lsty{j}); hold on
            xlabel('Cloud droplet concentration [cm^{-3}]'); 
            ylabel('Parcel displacement [m]'); 
            set(gca,'XTick',0:100:600)
            ylim([0 100]);
            title('b)')

            if j==1 
               ht = text(labposx(i), labposy(i),Labelstr{i});
               ht.Color=colors(i,:)
               ht_all(i)=ht
               set(ht,'Rotation',270)           

            end

            rdrops=CLOUDTOP.SPEC(2,CLOUDTOP.SPEC(2,:)>1e-6)*1e6;
            Ndrops=CLOUDTOP.SPEC(1,CLOUDTOP.SPEC(2,:)>1e-6);
            
            DROPLET_RADII(i,j,:) = rdrops;
            DROPLET_CONC(i,j,:)  = Ndrops;

            subplot(2,2,3:4)
            semilogx(rdrops,Ndrops,'color',colors(i,:), 'Linestyle',lsty{j},'Linewidth',2.2); hold on
            xlabel('Hydrometeor radius [\mum]'); 
            ylabel({'Hydrometeor number concentration';'size distribution [cm^{-3}]'})

            xlim([4 12]);
            set(gca,'YTick',[0:7]*1000,'YtickLabel',YTicklabs);
            title('c)')
            drawnow
        end
    end
end

 save(strcat('Figure2_Data',num2str(CLOUD_TOP),'.mat'),...
     'CDNC','SMAX','R_ACT','tau','alpha','LWP', 'DROPLET_RADII',...
     'DROPLET_CONC','ZPROFILES')

%Response of CDNC and optical properties to surface phase------------------
%delta_CDNC  = 100 * squeeze(CDNC(:,2,:) ./ CDNC(:,1,:)-1);
%delta_tau   = 100 * squeeze(tau(:,2,:) ./ tau(:,1,:)-1);
%delta_alpha = 100 * squeeze(alpha(:,2,:) ./ alpha(:,1,:)-1);
%diff_albedo = squeeze(alpha(:,2,:) - alpha(:,1,:));
%--------------------------------------------------------------------------

%Following supplementary of Topping et al. Nat. GeoSci (2013)--------------
%Assumes 0.7 global ocean and cloud coverage
%MA_F0  = 1370/4 * diff_albedo(1) * 0.7 * 0.7;
%HYY_F0 = 1370/4 * diff_albedo(2) * 0.3 * 0.5;
%FORCING = [MA_F0, HYY_F0];
%--------------------------------------------------------------------------
%return
%fixed albedo calcs
% fixed_delta_alpha = [1, 5, 10]; % [%]
% for i = 1:3
%     fixed_delta_CDNC(i,:) = 100 * ( exp( alpha(i,1) .*...
%         fixed_delta_alpha / 100 / 0.075 ) - 1);
% end
%%%%
return
% Finishing touches to plot-------------------------------------------------
 figure(1);subplot(2,2,1); hold on
 h = zeros(2, 1);
 h(1) = plot(-5,-5,'k');
 h(2) = plot(-5,-5,'--k');
 legend(h, 'Bulk Köhler','Compressed film','Location','Northeast');

 legend boxoff

%--------------------------------------------------------------------------



%save('HYY_CDNC_W_PDF.mat','CDNC','SMAX','R_ACT','tau','alpha','LWP')
%save('HYY_CDNC_W_PDFALL.mat')

fig = gcf;
set(findall(fig,'-property','FontSize'),'FontSize',9)
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8 5];
print('Figure2','-dpng','-r400')
