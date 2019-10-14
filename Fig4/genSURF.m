close all; clear all
%presribed synthetic sensitivity surface span
nRes   = 20;
N2_vec = logspace(1, 3, nRes); 
N1_vec = logspace(0, 4, nRes);
linw=1.5; lsty{1}='-'; lsty{2}=':';

ENV_STRING{1} = 'Marine';
ENV_STRING{2} = 'Boreal'
ENV_STRING{3} = 'NUM event'
NAT           = [2 2 2];
SURFACE_FLAG  = [0 1];
UPDRAFT       = 0.32;

for ienv = 2:2
    
	for in2 = 1:nRes
    	N2 = N2_vec(in2);
    
    	for in1 = 1:nRes
        	N1 = N1_vec(in1);
       		in1
        	in2
        	ienv
        	for isurf=1:2
            
        		N1N2_GenInputFiles(SURFACE_FLAG(isurf), 0.32, ENV_STRING{ienv}, N1, N2);
        		[CLOUDTOP, OPTICS] = N1N2_RUN_CLOUD(ENV_STRING{ienv}, NAT(ienv));
        		ALPHA(in1, in2, isurf)   = OPTICS.ALBEDO;
        		TAU(in1, in2, isurf)     = OPTICS.COD;
        		LW_PATH(in1, in2, isurf) = OPTICS.LWP;
        		CLOUDTOP.CDNC_rc 
        		OPTICS.ALBEDO
        		CDNC(in1, in2, isurf)        = CLOUDTOP.CDNC_rc; 
        		SMAX(in1, in2, isurf)        = CLOUDTOP.Smax;
        		RACT_OUT(in1, in2, isurf, :) = CLOUDTOP.RACT;      
        end
    end
    end
               
save(strcat('REV_HYY_', ENV_STRING{ienv}, '_OPTICS_N1N2_DELTACDNC.mat'),...
 'CDNC','SMAX','RACT_OUT','ALPHA','TAU','LW_PATH');
end

