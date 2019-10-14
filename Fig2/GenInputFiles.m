%/----------------------------------X------------------------------------/%
%Sam Lowe, ACES (SU) 2019: Function for writing input parameter and control
%flags files for the ICPM-Compressed film model (Lowe et al. Nature Comms)
%/----------------------------------X------------------------------------/%

function GenInputFiles(iCOMPFILM, UPDRAFT, ENV, CLOUD_TOP)

TEMP0       = 280;         %Initial temperature       [K] 
PRESS0      = 98000;       %Initial pressure          [Pa]
CLOUD_DEPTH = CLOUD_TOP;   %Simulation depth          [m]    
z0          = 150;         %Initial altitude          [m]
MASS_ACCOM  = 1.0;         %Water mass accommodation coefficient 
SURF_TENS_W = 72.8;        %Surface tension water     [mNm-1]
SIM_TIME    = 1e6;         %Total simulation time     [s]
%UPDRAFT    = 0.32;        %Updraft velocity          [ms-1]
RH0         = 99.9;        %Initial relative humidity [%] 

%Model compoennt densities
DENSITY                     = [1.841, 1.78, 1.77, 1.5, 1.5, 2., 2.65, 2.165]; %[gcm-3] 

%Set up log-normal aersol number concentration size distribution and
%mass fraction composition
if strcmp(ENV, 'Marine')

    cd MAV
    MODAL_PARS.CONC         = [223  137];   %[number/cm3]
    MODAL_PARS.GSD          = [1.68 1.68];     
    MODAL_PARS.GEOMEAN_DIAM = [0.0390  0.139];     %[um]
    DENSITY(4:5)            = 0.852;               %Set organic density. Palmitic acid 256.4 [gmol-1]
    MASS_FRAC               = [0 0.8 0 0.2 0 0 0 0.0;...
                               0 0.0 0 0.2 0 0 0 0.8];
    NAT                     = 2;
    NMODE                   = [1 1];

    %MASS_FRAC               = [0 0 0 0.2 0 0 0 0.8];
    %!NAT                     = 1;
    %NMODE                   = 2;
    
elseif strcmp(ENV, 'Marine Heintzenberg')

    MODAL_PARS.CONC         = [265   165];
    MODAL_PARS.GSD          = [1.45  1.50];
    MODAL_PARS.GEOMEAN_DIAM = [0.042 0.165];    
    DENSITY(4)              = [];                  %Set organic density
    MASS_FRAC               = [0 0 0 0.2 0 0 0 0.8];
    NAT                     = 1;

elseif strcmp(ENV, 'Boreal')

    cd HYY
%     MODAL_PARS.CONC         = [1010   451];
%     MODAL_PARS.GSD          = [1.71   1.58];
%     MODAL_PARS.GEOMEAN_DIAM = [0.0474 0.1796];
    MODAL_PARS.CONC         = [1110   540];
    MODAL_PARS.GSD          = [1.75   1.62];
    MODAL_PARS.GEOMEAN_DIAM = [0.0453 0.1644];    
    DENSITY(4:5)            = [1.2 1.4];                  %Set organic density
    DENSITY(1)              = 1.72;
    INORG_MASS_RATIO        = 0.1515/0.1559;  %ammonium sulfate:ammonium nitrate  mass
    FORG                    = 0.60;
    MASS_FRAC               = [(1-FORG)/(1+INORG_MASS_RATIO) 0 ...
                              INORG_MASS_RATIO*(1-FORG)/(1+INORG_MASS_RATIO) FORG...
                              0 0 0 0;...
                              (1-FORG)/(1+INORG_MASS_RATIO) 0 ...
                              INORG_MASS_RATIO*(1-FORG)/(1+INORG_MASS_RATIO) 0 ...
                              FORG 0 0 0];
    NAT                     = 2;
    NMODE                   = [1 1];

elseif strcmp(ENV, 'NUM event')

    cd NUM
    MODAL_PARS.CONC         = [2000  30];
    MODAL_PARS.GSD          = [1.71  1.703];
    MODAL_PARS.GEOMEAN_DIAM = [0.023 0.200];
    DENSITY(4:5)            = [1.2 1.24];                  %Set organic density
    MASS_FRAC               = [0 0 0.48 0.52 0 0 0 0;...
                               0 0 0.70 0.0 0.3 0 0 0];
    NAT                     = 2;
    NMODE                   = [1 1];

elseif strcmp(ENV, 'Coastal Lowe')                 %Mace Head Summer (May - Aug, 2012)

    MODAL_PARS.CONC         = [340  140];
    MODAL_PARS.GSD          = [1.54  1.62];
    MODAL_PARS.GEOMEAN_DIAM = [0.0474 0.1796];
    DENSITY(4)              = [];                  %Set organic density

end

SURF_TENS_ORG      = 40e-3; 
MIN_FILM_THICKNESS = 0.2e-9;



for i = 1:2

    AERO_DIST_PARAMS(i,:) = [MODAL_PARS.CONC(i)*1e6,...
        MODAL_PARS.GEOMEAN_DIAM(i)/2, MODAL_PARS.GSD(i)];

end




%Write to par.in-----------------------------------------------------------
fid_params = fopen('ini/par.in','w');
fprintf(fid_params, '%10.4e \n',     TEMP0, PRESS0, CLOUD_DEPTH, z0, MASS_ACCOM,...
                                     SURF_TENS_W, SIM_TIME, UPDRAFT, RH0);
fprintf(fid_params,'%d \n',          NAT);
    
jj = 0;
for i = 1:NAT
    %convert mass to volume fractions
    VOL_FRAC(i,:) = MASS_FRAC(i,:) ./ DENSITY ...
    / sum(MASS_FRAC(i,:)./DENSITY);

    fprintf(fid_params,'%d \n',      NMODE(i));
    for j = 1:NMODE(i)
        fprintf(fid_params,'%10.4e \t',  AERO_DIST_PARAMS(j+jj,:));
        fprintf(fid_params,'\n');
    end
    fprintf(fid_params,'%10.4e \t ', VOL_FRAC(i,:));
    fprintf(fid_params,'\n');
    jj = jj + 1;
end

fprintf(fid_params,'%10.4e \n',      SURF_TENS_ORG);      %Organic component surface tension
fprintf(fid_params,'%10.4e \n',      MIN_FILM_THICKNESS); %Minimum film thickness 
fclose(fid_params); 
%--------------------------------------------------------------------------

%Write to control.dat-----------------------------------------------------------
fid_flags  = fopen('ini/control.dat','w');
fprintf(fid_flags,'   ***** 1=constant updraft velocity and up/down 0=explicit buoyancy calculated');
fprintf(fid_flags,'\n');
fprintf(fid_flags, '%d \n', 1);
fprintf(fid_flags,'   ***** neutral buoyancy after cloud top is reached (only for explicit buoyancy; 1: yes)');
fprintf(fid_flags,'\n');
fprintf(fid_flags, '%d \n', 0);
fprintf(fid_flags,'"   ***** entrainment parameter (0=no entrainment, max. 0.6)"');
fprintf(fid_flags,'\n');
fprintf(fid_flags, '%d \n', 0);
fprintf(fid_flags,'   ***** initflag=1/2/3...      sets which aerosol input is desired');
fprintf(fid_flags,'\n');
fprintf(fid_flags, '%d \n', 5);
fprintf(fid_flags,'   ***** kprecip=0/1/2/3/... : switch for precip formation through coalescence');
fprintf(fid_flags,'\n');
fprintf(fid_flags, '%d \n', 0);
fprintf(fid_flags,'   ***** irainfall: 0/1 gravitational settling of precip sized drops');
fprintf(fid_flags,'\n');
fprintf(fid_flags, '%d \n', 0);
fprintf(fid_flags,'   ***** iupadiabat: 0/1 choose between adiabatic (1) or pseudo-adiabatic ascent (0)');
fprintf(fid_flags,'\n');
fprintf(fid_flags, '%d \n', 1);
fprintf(fid_flags,'   ***** iCOMPFILM: 1 = compressed film surface phase, 0 = standard kohler ');
fprintf(fid_flags,'\n');
fprintf(fid_flags, '%d \n', iCOMPFILM);
fclose(fid_flags); 
%--------------------------------------------------------------------------
cd ..

end

