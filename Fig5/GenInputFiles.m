%/----------------------------------X------------------------------------/%
%Sam Lowe, ACES (SU) 2019: Function for writing input parameter and control
%flags files for the ICPM-Compressed film model (Lowe et al. Nature Comms)
%/----------------------------------X------------------------------------/%

function N1N2_GenInputFiles(iCOMPFILM, ENV, PARS, CLOUD_DEPTH)

TEMP0       = 280;         %Initial temperature       [K] 
PRESS0      = 98000;       %Initial pressure          [Pa]
CLOUD_DEPTH = CLOUD_DEPTH; %Simulation depth          [m]    
z0          = 150;         %Initial altitude          [m]
MASS_ACCOM  = 1.0;         %Water mass accommodation coefficient 
SURF_TENS_W = 72.8;        %Surface tension water     [mNm-1]
SIM_TIME    = 1e6;         %Total simulation time     [s]
UPDRAFT     = PARS(9);        %Updraft velocity          [ms-1]
RH0         = 99.9;        %Initial relative humidity [%] 

%Model compoennt densities
DENSITY                     = [1.841, 1.78, 1.77, 1.5, 1.5, 2., 2.65, 2.165]; %[gcm-3] 

%Set up log-normal aersol number concentration size distribution and
%mass fraction composition
if strcmp(ENV, 'Marine')
    
    cd MAV
    DENSITY(4:5)            = 0.852;               %Set organic density. Palmitic acid 256.4 [gmol-1]
    MASS_FRAC               = [0 1 - PARS(7) 0 PARS(7) 0 0 0 0.0;...
                               0 0.0 0 PARS(8) 0 0 0 1 - PARS(8)];
    

elseif strcmp(ENV, 'Boreal')

    cd HYY
    DENSITY(4:5)            = [1.2 1.4];                  %Set organic density
    DENSITY(1)              = 1.72;
    INORG_MASS_RATIO        = 0.1515/0.1559;  %ammonium sulfate:ammonium nitrate  mass
    MASS_FRAC               = [(1-PARS(7))/(1+INORG_MASS_RATIO) 0 ...
                              INORG_MASS_RATIO*(1-PARS(7))/(1+INORG_MASS_RATIO) PARS(7)...
                              0 0 0 0;...
                              (1-PARS(8))/(1+INORG_MASS_RATIO) 0 ...
                              INORG_MASS_RATIO*(1-PARS(8))/(1+INORG_MASS_RATIO) 0 ...
                              PARS(8) 0 0 0];

end
MIN_FILM_THICKNESS = PARS(10);
SURF_TENS_ORG      = PARS(11); 

for i = 1:2
    %convert mass to volume fractions
    VOL_FRAC(i,:) = MASS_FRAC(i,:) ./ DENSITY ...
    / sum(MASS_FRAC(i,:) ./ DENSITY);
end

%Write to par.in-----------------------------------------------------------
fid_params = fopen('ini/par.in','w');
fprintf(fid_params, '%10.4e \n',     TEMP0, PRESS0, CLOUD_DEPTH, z0, MASS_ACCOM,...
                                     SURF_TENS_W, SIM_TIME, UPDRAFT, RH0);
    fprintf(fid_params,'%d \n',      2);
    fprintf(fid_params,'%d \n',      1);
    fprintf(fid_params,'%10.4e \t',  PARS(1:3));
    fprintf(fid_params,'\n');
    fprintf(fid_params,'%10.4e \t ', VOL_FRAC(1,:));
    fprintf(fid_params,'\n');
    fprintf(fid_params,'%d \n',      1);
    fprintf(fid_params,'%10.4e \t',  PARS(4:6));
    fprintf(fid_params,'\n');
    fprintf(fid_params,'%10.4e \t ', VOL_FRAC(2,:));
    fprintf(fid_params,'\n');
    fprintf(fid_params,'%10.4e \n',  SURF_TENS_ORG);      %Organic component surface tension
    fprintf(fid_params,'%10.4e \n',  MIN_FILM_THICKNESS); %Minimum film thickness  
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

