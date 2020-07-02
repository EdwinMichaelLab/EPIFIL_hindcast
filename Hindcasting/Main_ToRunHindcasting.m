clear all; clc;

%% add folders containing model functions to path
addpath('../CommonFunctions/');
addpath('../BaselineFitting');
addpath('../Interventions/');

%% User-defined inputs

% Parameters related to sampling and resampling
SIR_samples = 500; % default: 500

% try to load user-editable variables from IO/setup_Vars.m
% if no file exists, spit error message and die
try
    load('../IO/IN/Baseline_IN.mat');
catch ME
    if (strcmp(ME.identifier, 'MATLAB:load:couldNotReadFile'))
        fprintf(1, '\nError: Could not read Baseline Input file from IO/.\nDid you run setup_Vars.m?\n');
        exit(2);  % Die with value of 2
    else
        rethrow(ME)
        exit(3); % Die with value of 3, unknown error
    end
end

% Import country-level age-demographic parameters
Import_CountryDemogParams;

%% Intervention specifications

% Load user defined intervention variables within IO directory
load('../IO/IN/Inter_IN.mat');

% Mass drug administration parameters
AgeLimits=[5 100]; % treatment for > 5 yrs old
RegimenEfficacy0 = [0.1 0.99 9; % IVM worm kill rate, mf kill rate, sterilization period (months)
    0.55 0.95 6; % DEC+ALB
    0.55 1 120; % IVM+DEC+ALB (permanent sterilization - 10 yrs)
    0.35 0.99 9; % IVM+ALB
    0.55 0.95 6; % ALB
    0.59 0.86 10;% DEC salt
    0.45 0.99 9]; %DEC+IVM

% Vector control parameters
AnnualDecrease = 0;
IRSParams = [87.5, 93, 277, 6, 80]; % Vector control by indoor residual spraying
ITNParams = [20, 90, 97, 26, 12*3, 80]; % Vector control by long lasting insecticide nets
VCparams = SSA_IRS_ITN_Parameters(IRSParams,ITNParams,AnnualDecrease);

% Set up parallel pool
setupParallelPool();

%% Run fitting procedure for each site

for iSites = 1:length(Sites)
    
    %% initial baseline simulation     
    load(sprintf('../IO/OUT/ParamVectors%s.mat',Sites{iSites}));

    %% simulate interventions
    
    % define intervention data
    MfData  = eval(sprintf('%sMf',Sites{iSites})); % longitudinal surveys
    CfaData  = eval(sprintf('%sCfa',Sites{iSites})); % longitudinal surveys
    L3Data  = eval(sprintf('%sL3',Sites{iSites})); % longitudinal surveys
    Vol     = eval(sprintf('%sVol',Sites{iSites}));
    MDAReg = eval(sprintf('%sReg',Sites{iSites}));
    MDAFreq = eval(sprintf('%sFreq',Sites{iSites}));
    MDACov  = eval(sprintf('%sMDACov',Sites{iSites}));
    VC  = eval(sprintf('%sVC',Sites{iSites}));
    IRSCov  = eval(sprintf('%sIRSCov',Sites{iSites}));
    ITNCov  = eval(sprintf('%sITNCov',Sites{iSites}));
    ITNYear  = eval(sprintf('%sITNYear',Sites{iSites}));
    NumYears  = eval(sprintf('%sNumYears',Sites{iSites}));
    SwitchYear  = eval(sprintf('%sSwitchYear',Sites{iSites}));
    APYear  = eval(sprintf('%sYear',Sites{iSites})); % first age profile
    APMf = eval(sprintf('%sMfAge',Sites{iSites})); % first age profile
    APCfa = eval(sprintf('%sCfaAge',Sites{iSites})); % first age profile
    
    % convert data into format compatible with model
    if isnan(APMf)
        ageMax = 69; % default max age in a population
    else
        ageMax = max(APMf(:,4));
    end
    
    % correct mf data for blood sample volume
    if ~isnan(MfData(1,1)) && (Vol == 20)
        MfData(:,3)=floor(min(MfData(:,2),MfData(:,3)*1.95));
        if ~isnan(APMf)
            APMf(:,3)=floor(min(APMf(:,2),APMf(:,3)*1.95));
        end
    elseif ~isnan(MfData(1,1)) && (Vol == 100 || Vol == 60)
        MfData(:,3)=floor(min(MfData(:,2),MfData(:,3)*1.15));
        if ~isnan(APMf)
            APMf(:,3)=floor(min(APMf(:,2),APMf(:,3)*1.15));
        end
    end
    
    % Month to switch drug treatment regimen
    SwitchMonth = zeros(1,length(SwitchYear));
    for h = 1:length(SwitchYear)
        if SwitchYear(h) == 0
            SwitchMonth(h) = 0;
        else
            SwitchMonth(h) = SwitchYear(h)*12;
        end
    end
    
    % MDA Frequency
    MDAInterval = zeros(1,NumYears*12);
    if length(MDAFreq)==1
        MDAInterval(1,1:end) = MDAFreq(1);
    else
        MDAInterval(1,1:SwitchMonth) = MDAFreq(1);
        MDAInterval(1,SwitchMonth+1:end) = MDAFreq(2);
    end
    
    % MDA Drug Regimen Efficacy
    RegimenEfficacy = zeros(length(MDAReg),length(RegimenEfficacy0(1,:)));
    for j = 1:length(MDAReg)
        RegimenEfficacy(j,:) = RegimenEfficacy0(MDAReg(j),:);
    end
    
    % estimate worm kill parameter from data
    WormKill = [];
    for j = 1:length(MDAReg)
        if MDAReg(j) == 1
            WormKill(:,j) = 0 + 0.2*rand(length(mfPrevArray(1,:)),1);
        elseif MDAReg(j) == 2
            WormKill(:,j) = 0.45 + 0.2*rand(length(mfPrevArray(1,:)),1);
        elseif MDAReg(j) == 4
            WormKill(:,j) =  0.25 + 0.2*rand(length(mfPrevArray(1,:)),1);
        end
    end

    % MDA Coverage
    MonthlyMDACov = zeros(NumYears*12,length(mfPrevArray(1,:)));
    for i = 1:length(MDACov)
        MonthlyMDACov((i-1)*MDAFreq(1)+1:i*MDAFreq(1),:) = MDACov(i)*(0.75 + 0.25*rand(MDAFreq(1),length(mfPrevArray(1,:))));
    end

    % VC coverages
    ITNCoverages=[zeros(1,(ITNYear-1)*12),ones(1,length((ITNYear-1)*12+1:NumYears*12)).*ITNCov];
    IRSCoverages=ones(1,NumYears*12)*IRSCov;
    
    MultiVecMBR = ABR/12;
    kId = 1:length(ABR);
    ageMthMax = ageMax*12;
    da = 1;
    
    [mfAgePrevArray,cfaAgePrevArray,mfPrevIntv,cfaPrevIntv,MBRIntv,L3IntvPrev] =...
        Modelling_MDAplusVC_L3prev(kId,...
        ParameterVectors,L3Values,demog,ageMthMax,da,bCulex,...
        AgeLimits,RegimenEfficacy,MDAInterval,NumYears,...
        MultiVecMBR,IRSCoverages,ITNCoverages,MonthlyMDACov,...
        SwitchMonth,VCparams,WormKill,APYear);
    
    %% fit to post-intervention data
    
    [mfPrevIntv,L3IntvPrev,MBRIntv,cfaPrevIntv,ParameterVectors,...
        mfPrevArray,cfaPrevArray,cfaAgePrevArray,mfAgePrevArray,APMf,APCfa,WormKill] = ...
        fitPostIntvData_group_ageprofiles(MfData,CfaData,APMf,mfAgePrevArray,demog,da,...
        ageMthMax,SIR_samples,mfPrevIntv,L3IntvPrev,MBRIntv,cfaPrevIntv,...
        ParameterVectors,mfPrevArray,cfaPrevArray,cfaAgePrevArray,APCfa,APYear,...
        WormKill);
     
      save(sprintf('../IO/OUT/Hindcast%s.mat',char(Sites{iSites})),...
          'L3IntvPrev','MBRIntv','mfPrevIntv','cfaPrevIntv',...
          'ParameterVectors','mfPrevArray', 'cfaPrevArray',...
          'mfAgePrevArray','cfaAgePrevArray','bCulex','demog',...
          'ageMthMax','APMf','APCfa','WormKill');
     
    
end


