function [Nmode_HS, delt_HS, delt_d_threshold] = GETN1N2_MH

%--------------------------------------------------------------------------
%Load model output for scatter runs----------------------------------------
%load('0117_MH_TS_DELTA.mat','CDNC'); delt_d_fit=100*(CDNC(2,:)./CDNC(1,:)-1);
%load('../MACEHEAD_TS/MACEHEAD_TS_CDNC_RACT_SMAX.mat','CDNC'); 
load('../MACEHEAD_TS/200_MACEHEAD_TS_CDNC_RACT_SMAX.mat','CDNC'); 

delt_d_fit = 100 * (CDNC(2,:) ./ CDNC(1,:) - 1);
%--------------------------------------------------------------------------
%RAW DATA------------------------------------------------------------------
load('../MACEHEAD_TS/MH_PSD_MJJA.mat');
Dp          = MHPSDMJJA(1, 3:3:(end - 1)) * 1e-9; 
dnlogdd_raw = MHPSDMJJA(2:end, 3:3:(end - 1));
%--------------------------------------------------------------------------

%Fill in missing time-stamps of PSD with NaN values
%Define required time stamps based proper timestamp differential
UTC_req   = [(MHPSDMJJA(2,1)+MHPSDMJJA(2,2))/2:0.041666499999991:...
    (MHPSDMJJA(end,1)+MHPSDMJJA(end,2))/2];
%Initialise PSD matrix to be "inflated" with missing timesttamps
HS_index_req = NaN.*ones(length(UTC_req));
delt_d_req   = NaN.*ones(length(UTC_req));
Nmode_req    = NaN.*ones(length(UTC_req),2);
IDX_req      = NaN.*ones(length(UTC_req));

%Existting UTC times - average start and end
UTC_num   = (MHPSDMJJA(2:end,1)+MHPSDMJJA(2:end,2))/2;

%load bi-modal fitted logmormal parameters----------------------
load('../MACEHEAD_TS/MJJA_MH_2MODE_AEROFIT','-ascii');
Num_modes = MJJA_MH_2MODE_AEROFIT(:,19);    %number of modes fitted
GSD       = MJJA_MH_2MODE_AEROFIT(:,7:8);   %Geometric standard deviation
Dmode     = MJJA_MH_2MODE_AEROFIT(:,2:3);   %Mode diameter [m]
Nmode     = MJJA_MH_2MODE_AEROFIT(:,12:13); %Mode concentration [cm-3]
LLS       = MJJA_MH_2MODE_AEROFIT(:,17);    %linear least squares of fits
LOGLS     = MJJA_MH_2MODE_AEROFIT(:,18);    %logarithmic LLS of fits

for i=1:length(UTC_num)
    [dt(i), dt_ind(i)]   = min(abs(UTC_num(i)-UTC_req));
    delt_d_req(dt_ind(i))  = delt_d_fit(i);
    Nmode_req(dt_ind(i),:) = Nmode(i,:);
end

%Calculate average modal sizes and GSDs to run synthetic N space on
%Remove missing time stamped and unimodal fits
%Nmode = Nmode(~isnan(Dmode(:,1)) &~ isnan(Dmode(:,2)),:);
%GSD   = GSD(~isnan(Dmode(:,1))   &~ isnan(Dmode(:,2)),:);
%Dmode = Dmode(~isnan(Dmode(:,1)) &~ isnan(Dmode(:,2)),:);
%Remove parameters maybe associated with european influnce N2>Ncut [cm^-3]

% Ncut=300;
% GSD_cut=GSD(Nmode(:,2)<Ncut,:); 
% Dmode_cut=Dmode(Nmode(:,2)<Ncut,:); 
% Nmode_cut=Nmode(Nmode(:,2)<Ncut,:);
% prcs_N_2     = prctile(Nmode_cut, [5 95]);
% prcs_GSD_2   = prctile(GSD_cut, [5 95]);
% prcs_Rmode_2 = prctile(Dmode_cut, [5 95])*0.5e9;
% med_N_2      = median(Nmode_cut)
% med_GSD_2    = median(GSD_cut)
% med_Rmode_2  = median(Dmode_cut)*0.5e9

nbin        = 400;
Nall        = zeros(size(Nmode,1),2); 
dnlogdd_fit = zeros(size(Nmode,1),nbin);
%Adjust fittedconcentration-------------
iadj = 1
if iadj == 1
    for i=1:size(Nmode,1)
        Dp_meas=logspace(log10(Dp(1)), log10(Dp(end)), nbin);
        for imode = 1:Num_modes(i)
            [dnlogdd_mode,dummy] = DMPS_PSD(Dp_meas, Dmode(i,imode),GSD(i,imode), Nmode(i,imode));
            Nall(i,imode)        = sum(dnlogdd_mode(2:end).*diff(log10(Dp_meas)));
            dnlogdd_fit(i,:)     = dnlogdd_mode+dnlogdd_fit(i,:);
        end
    end
else
    for i=1:size(Nmode,1)
        Dp_meas=logspace(log10(Dp(1)), log10(Dp(end)), nbin);
        for imode = 1:Num_modes(i)
            [dnlogdd_mode,dummy] = DMPS_PSD(Dp_meas, Dmode(i,imode),GSD(i,imode), Nmode(i,imode));
            dnlogdd_fit(i,:)=dnlogdd_mode+dnlogdd_fit(i,:);
        end
    end
    Nall = Nmode;
end

%Remove anthropegenics in HS set of PSDs
icut = 1
if icut == 1
    n2cut_index=Nall(:,2) < 1e9;  %300
    Nall=Nall(n2cut_index,:); delt_d_fit=delt_d_fit(n2cut_index)
    delt_d_req=delt_d_req(n2cut_index)
    GSD=GSD(n2cut_index,:); Dmode=Dmode(n2cut_index,:);
    dnlogdd_fit=dnlogdd_fit(n2cut_index,:);
    dnlogdd_raw=dnlogdd_raw(n2cut_index,:);
    
    n1cut_index=Nall(:,1)< 1e9 %1190
    Nall=Nall(n1cut_index,:); delt_d_fit=delt_d_fit(n1cut_index)
    delt_d_req=delt_d_req(n1cut_index)
    GSD=GSD(n1cut_index,:); Dmode=Dmode(n1cut_index,:);
    dnlogdd_fit=dnlogdd_fit(n1cut_index,:);
    dnlogdd_raw=dnlogdd_raw(n1cut_index,:);
end

%Define high sensitivity ((HS) criteria
delt_d_threshold = 0 %prctile(delt_d_fit,50);
%delt_d_threshold = 10; %percentage CDNC enhancement to give 5% albedo change
HS_index=delt_d_fit>=delt_d_threshold; HS_index2=delt_d_req>=delt_d_threshold;
HS_delt_d_fit=delt_d_fit(HS_index); LS_delt_d_fit=delt_d_fit(~HS_index);
PSDs_highsens_fit=dnlogdd_fit(HS_index,:);
PSDs_highsens_raw=dnlogdd_raw(HS_index,:);

%High sens parameters and sens
Nmode_HS = Nall(HS_index, :);    Nmode_LS = Nall(~HS_index, :);
Dmode_HS = Dmode(HS_index, :);    Dmode_LS = Dmode(~HS_index, :);
GSD_HS = GSD(HS_index, :);      GSD_LS = GSD(~HS_index, :);
delt_HS  = delt_d_fit(HS_index); delt_LS  = delt_d_fit(~HS_index)

%Sort ascedning delta_CDNC for plotting vissibbility
isort = 1
if isort ==1
    [delt_HS, IDX_delta_CDNC] = sort(delt_HS);
    Nmode_HS = Nmode_HS(IDX_delta_CDNC,:);
end

end
