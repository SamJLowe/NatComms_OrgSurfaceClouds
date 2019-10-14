close all
clear all

ENV_STRING{1} = 'Marine';
ENV_STRING{2} = 'Boreal';
SURFACE_FLAG  = [0 1];
CLOUD_DEPTH = 200;
if CLOUD_DEPTH >= 200
    doOPTICS = true
end 

%Parameter ranges to sample from
%                N1, R1, GSD1, N2, R2, GSD2, FORG1, FORG2, w, delta_film, ST org
PAR.MIN(1, :) = [150,  15.5, 1.4,  60,  70,  1.4,  0.05, 0.05, 0.05, 0.16, 30];
PAR.MAX(1, :) = [600,  23.5, 1.8,  250, 100, 1.8,  0.40, 0.40, 2.40, 0.30, 50];
% PAR.MIN(2, :) = [450,  16,   1.55, 215, 75,  1.4,  0.45, 0.45, 0.05, 0.16, 30];
% PAR.MAX(2, :) = [1600, 32,   1.9,  690, 105, 1.75, 0.75, 0.75, 2.40, 0.30, 50];
PAR.MIN(2, :) = [440,  14,   1.56,  160,  66,  1.45,  0.45, 0.45, 0.05, 0.16, 30];
PAR.MAX(2, :) = [1780, 31,   1.95,  920,  98.5, 1.80, 0.75, 0.75, 2.40, 0.30, 50];
%Conversion to model-friendly 98.5
UNIT_CONV = [1e6, 1e-3, 1, 1e6, 1e-3, 1, 1, 1, 1, 1e-9, 1e-3 ];
UNIT_CONV = repmat(UNIT_CONV, 2, 1);
PAR.MIN   = PAR.MIN .* UNIT_CONV; 
PAR.MAX   = PAR.MAX .* UNIT_CONV; 
%Set up LHS sampling
NPAR  = length(PAR.MAX);
NSAMP = 5000
%Sample quotients [0,1]
QUOT_ARRAY = lhsdesign(NSAMP, NPAR);
%Save LHS quotients for reproducibility
save(strcat(num2str(NSAMP),'_dz',num2str(CLOUD_DEPTH),'m_LHSsample.mat'),...
    'QUOT_ARRAY');
tmp = [2, 1]
%Call cloud model over environs and sampled parameters
for itmp = 1:1
   ienv = tmp(itmp)
     for ip = 1: NSAMP
        %LHS quotients to input parameters
        PAR_IN   = PAR.MIN(ienv, :) + (PAR.MAX(ienv, :) - PAR.MIN(ienv, :)) .* QUOT_ARRAY(ip,:);

        for isurf = 1:2
            ienv
            ENV_STRING{ienv}
            ip
            
            GenInputFiles(SURFACE_FLAG(isurf), ENV_STRING{ienv}, PAR_IN, CLOUD_DEPTH);
            [CLOUDTOP, OPTICS] = RUN_CLOUD(ENV_STRING{ienv}, doOPTICS);
            ALPHA(isurf, ip) = OPTICS.ALBEDO;
            TAU(isurf, ip)   = OPTICS.COD;
            LIQWATPATH(isurf, ip) = OPTICS.LWP;
            CLOUDTOP.CDNC_rc
            CDNC(isurf, ip)        = CLOUDTOP.CDNC_rc;
            SMAX(isurf, ip)        = CLOUDTOP.Smax;
            RACT_OUT(isurf, ip, :) = CLOUDTOP.RACT;
        end       
    end
    if mod(ip,1000) == 0
       save(strcat('REV_HYY_',num2str(NSAMP),'_',ENV_STRING{ienv},'_OPTICS_LHSOUTPUT.mat'),...
            'CDNC','SMAX','RACT_OUT','QUOT_ARRAY', 'LIQWATPATH', 'ALPHA', 'TAU')
    end
end

