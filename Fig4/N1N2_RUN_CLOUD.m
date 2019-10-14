function [CLOUDTOP, OPTICS] = N1N2_RUN_CLOUD(ENV, NAT)

if strcmp(ENV, 'Marine')
    cd MAV
elseif strcmp(ENV, 'Boreal')
    cd HYY
elseif strcmp(ENV, 'NUM event')
    cd NUM
end

%remove previous model output----------------------------------------------
unix('rm -f uit/zplot_micro');
unix('rm -f uit/zplot_dyn');
unix('rm -f bulk_cdnc.txt');
%unix('rm -f MC_output.txt');
%-------------------------------------------------------------------------

%Run  model executable----------------------------------------------------
unix('./bubble.x > /dev/null');
%unix('./bubble.x');
%-------------------------------------------------------------------------
bulk_cdnc             = load('bulk_cdnc.txt');
%Loading model output from bulk_cdnc.txt-----------------------------------
CLOUDTOP.CDNC_1micron = bulk_cdnc(1);                                                %1 micron integration limit    [#/cm-3]
CLOUDTOP.CDNC_sc      = bulk_cdnc(2);                                                %CCN at S. Sc<Smax             [#/cm-3]
CLOUDTOP.CDNC_rc      = bulk_cdnc(3);                                                %Kohler rc<rwet                [#/cm-3]
CLOUDTOP.Smax         = bulk_cdnc(4);                                                %Maximum parcel supersat.      [%]
CLOUDTOP.RACT         = [bulk_cdnc(7) bulk_cdnc(8)]*1e9;                                                 %Samllest activated size class [nm]

%Read microphysics and dynamics profiles-----------------------------------
micro_id = fopen('uit/zplot_micro');
nCols=15;
if NAT == 2
    nCols = 20;
end
format = repmat('%f', [1 nCols]);                                % Number of columns and read format for doubles in each column
micro_cell = textscan(micro_id, format, 'HeaderLines',2);                  % read data in, skipping column headings
fclose(micro_id);                                                          % close file
micro_array=cell2mat(micro_cell);                                          % convert data cell to array
micro_array(end-1:end,:)=[];                                               % remove end rows (nothing useful there)

dyn_id=fopen('uit/zplot_dyn');                                             % repeat for zplot_dyn
nCols=11;
format = repmat('%f', [1 nCols]);
dyn_cell = textscan(dyn_id, format, 'HeaderLines',2);
fclose(dyn_id);
dyn_array=cell2mat(dyn_cell);
dyn_array(end-1:end,:)=[];

%remove miscellaneous values..... (To check: why do these -99.0 vlues get
%saved)
nan_micro.index=micro_array(:,1)==-99.0;
micro_array(nan_micro.index==1,:)=[];
nan_dyn.index=dyn_array(:,1)==-99.0;
dyn_array(nan_dyn.index==1,:)=[];

if NAT == 1

    T       = micro_array(:,11);                                               %Temperature [K]
    P       = micro_array(:,12);                                               %Pressure    [hPa]
    t       = micro_array(:,1);                                                %time        [s]
    Nd      = micro_array(:,7);                                                %droplet concentration
    Nd      = (Nd.*P)./(T*287)*1e-3;                                           %Number concentration dropets, defined r>1um       [#/cm-3]
    Nd_sc   = micro_array(:,13); Nd_sc = (Nd_sc.*P)./(T*287)*1e-3;             %Number concentration dropets, defined Smax>sc     [#/cm-3]
    Nd_rc   = micro_array(:,14); Nd_rc = (Nd_rc.*P)./(T*287)*1e-3;             %Number concentration dropets, defined r>rc & s>sc [#/cm-3]
    reff    = micro_array(:,15);                                               %Effective radius [m]

    dz      = 1000*(dyn_array(:,2)-dyn_array(1,2));                            %Parcel displacement  [m]
    u       = dyn_array(:,3);                                                  %updraft              [ms-1]
    ss      = dyn_array(:,4);                                                  %supersaturation      [%]
    zrad    = dyn_array(:,10);                                                 %parcel radius        [m]
    rho_air = dyn_array(:,9);                                                  %density of air       [kg m-3]
    lwc     = dyn_array(:,6);                                                  %Liquid water content [g_water kg_air-1]
    z       = 1000*(dyn_array(:,2)-dyn_array(1,2));                                             %Altitude             [m]


    %droplet spectra ----------------------------------------------------------
%     MC_out=load('MC_output.txt');
%     N=MC_out(5,3:end); rd=MC_out(2,3:end);                                      %Concentration [#/cm-3], radius [m]
%     dlogdr=diff(log10(rd)); dN=N(2:end)./dlogdr;                                %Normalise against bin-width
%     CLOUDTOP.SPEC = [dN; rd(2:end)];
    %--------------------------------------------------------------------------

elseif NAT == 2
    T     = micro_array(:,17);
    P     = micro_array(:,18);
    t     = micro_array(:,1);
    lwc   = micro_array(:,7);
    Nd    = micro_array(:,14);
    Nd    = (Nd.*P)./(T*287)*1e-3;
    Nd_sc = micro_array(:,15); Nd_sc = (Nd_sc.*P)./(T*287)*1e-3;
    Nd_rc = micro_array(:,16); Nd_rc = (Nd_rc.*P)./(T*287)*1e-3;
    reff  = micro_array(:,19);

    dz      = 1000*(dyn_array(:,2)-dyn_array(1,2));
    u       = dyn_array(:,3);                                                          %updraft (m/s) (should/fa be constant for now)
    ss      = dyn_array(:,4);                                                         %supersaturation (%)
    zrad    = dyn_array(:,10);
    rho_air = dyn_array(:,6);                                                     %rho air
    z       = (micro_array(:,2) - micro_array(1,2))*1000;                                                        %altitude (km)

    %------------------------------------------------------------------------
    %droplet spectra
%     MC_out=load('MC_output.txt');
%     N=MC_out(5:4:end,3:end); rd=MC_out(2:4:end,3:end);
%     rcrit = MC_out(3:4:end,3:end);
%     figure(100)
%     plot(rd','.')
% 
%     figure(101)
%     plot(rd(2,:),'b'); hold on
%     plot(rcrit(2,:),'r');
% 
%     dlogdr(1,:)=diff(log10(rd(1,:)));
%     dlogdr(2,:)=diff(log10(rd(2,:)));
%     dNdlogr(1,:)=N(1,2:end)./dlogdr(1,:);
%     dNdlogr(2,:)=N(2,2:end)./dlogdr(2,:);
%     %interp to same grid
%     %rmin=min(min(rd(rd>1e-6))); rmax=max(max(rd(rd>1e-6)));
%     rmin1=min(rd(1,rd(1,:)>1e-6)); rmax1=max(rd(1,rd(1,:)>1e-6));
%     rmin2=min(rd(2,rd(2,:)>1e-6)); rmax2=max(rd(2,rd(2,:)>1e-6));
%     rinterp=logspace(log10(rmin1), log10(rmax2), 200);
%     ninterp(1,:)=interp1(rd(1,2:end),dNdlogr(1,:),rinterp);
%     ninterp(2,:)=interp1(rd(2,2:end),dNdlogr(2,:),rinterp);
%     ninterp(isnan(ninterp)) =0;
%     CLOUDTOP.SPEC = [sum(ninterp,1); rinterp];
% 
%     figure(5)
%     semilogx(rd(1,2:end)*1e6,(N(1,2:end)./diff(log10(rd(1,:)))),'r'); hold on
%     semilogx(rd(2,2:end)*1e6,(N(2,2:end)./diff(log10(rd(1,:)))),'b')
%     semilogx(rinterp*1e6,sum(ninterp,1),'k')
%     xlim([2 10])
    %------------------------------------------------------------------------

    %RACT                  = [bulk_cdnc(7) bulk_cdnc(8)]*1e9;                                                 %Samllest activated size class [nm]

end
    %Concatenate profiles for output
    %PROFILES    = [z, Nd, Nd_sc, Nd_rc, ss, lwc, T, P, u, reff];



%Calculate optical properties----------------------------------------------
rhow      = 1e6;                                                            %Density water [g/m3]
g         = 0.85;                                                           %Kravitz 2014, asymmetry parameter [dimensionless]
LWP_trapz = trapz(z, rho_air .* lwc);                                       %Numerical integration to get liquid water path [g m-2]
OPTICS.COD       = 3 / (2 * rhow) * LWP_trapz /reff(end);                   %Cloud optical thickness Stephens 1978 [dimensionless]
OPTICS.ALBEDO    = (1 - g) * OPTICS.COD   / (2 + (1 - g) * OPTICS.COD  );   %Cloud albedo Bohren 1987 [dimensionless]
OPTICS.LWP       = LWP_trapz;
cd ..
end
