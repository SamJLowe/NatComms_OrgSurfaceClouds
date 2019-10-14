close all; clear all
%presribed synthetic sensitivity surface span
nRes         = 20;
FORG_vec     = linspace(0.05, 0.95, nRes);
UPDRAFT_vec  = linspace(0.10, 10.0, nRes);
N2_ARRAY     = [30 134.600; 160 540; 30 134.600]

ENV_STRING{1} = 'Marine';
ENV_STRING{2} = 'Boreal'
ENV_STRING{3} = 'NUM event'
SURFACE_FLAG  = [0 1];

for ienv = 2:2
    
    for in2 = 1:2
        N2 = N2_ARRAY(ienv,in2);
        
        for iFORG = 1:nRes
            FORG = FORG_vec(iFORG);
            
            for iw = 1:nRes
                UPDRAFT = UPDRAFT_vec(iw);
                ienv
                in2
                iFORG
                iw
                
                for isurf=1:2
                    
                    GenInputFiles(SURFACE_FLAG(isurf), UPDRAFT, ENV_STRING{ienv}, N2, FORG);
                    [CLOUDTOP, RACT] = RUN_CLOUD(ENV_STRING{ienv}, 2);
                    
                    CLOUDTOP.CDNC_rc
                    CDNC(in2, iw, iFORG, isurf) = CLOUDTOP.CDNC_rc;
                    CDNC2(in2, iw, iFORG, isurf) = CLOUDTOP.CDNC_1micron;
                    SMAX(in2, iw, iFORG, isurf) = CLOUDTOP.Smax;
                    RACT_OUT(in2, iw, iFORG, isurf, :) = RACT;
                    
                    %figure(3)
                    %subplot(1,2,1)
                    %plot(PROFILES(:,5), PROFILES(:,1),'r','Linewidth',linw, 'Linestyle',lsty{isurf}); hold on;
                    %xlabel('Supersaturation [%]'); ylabel('Parcel displacement [m]');

                    %subplot(1,2,2)
                    %plot(PROFILES(:,4), PROFILES(:,1),'k','Linewidth',linw, 'Linestyle',lsty{isurf}); hold on
                    %xlabel('Cloud droplet concentration [cm^{-3}]'); ylabel('Parcel displacement [m]');

                    drawnow
                    
                end
            end
        end
    end
    save(strcat('HYY_REV_',ENV_STRING{ienv},'_FWN2_DELTACDNC_NEW.mat'),...
        'CDNC','SMAX', 'FORG_vec', 'UPDRAFT_vec', 'N2_ARRAY','CLOUDTOP',...
        'RACT_OUT','CDNC2');
end

