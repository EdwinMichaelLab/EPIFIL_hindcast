clear all; clc;

%% add CommonFunctions folder containing basic model functions and Baseline folder containing baseline fits to path
addpath('../CommonFunctions/');
addpath('../BaselineFitting/Baseline_MATLAB/');
addpath('../Hindcasting/');

%% Load user defined variables within IO directory
load('../IO/IN/Inter_IN.mat');

% integration time step
da = 1; % 1 month

%% Intervention specifications

% Mass drug administration parameters
AgeLimits = [5 100]; % treatment for > 5 yrs old
RegimenEfficacy0 = [0.1 0.99 9; % IVM worm kill rate, mf kill rate, sterilization period (months)
    0.55 0.95 6; % DEC+ALB
    0.55 1 120; % IVM+DEC+ALB (permanent sterilization - 10 yrs)
    0.35 0.99 9; % IVM+ALB
    0.55 0.95 6; % ALB
    0.59 0.86 10]; % DEC salt

% Vector control parameters
AnnualDecrease = 0;
IRSParams = [87.5, 93, 277, 6, 80]; % Vector control by indoor residual spraying
ITNParams = [20, 90, 97, 26, 12*3, 80]; % Vector control by long lasting insecticide nets
VCparams = SSA_IRS_ITN_Parameters(IRSParams,ITNParams,AnnualDecrease);

setupParallelPool(); % Initialize the parallel workers

%% Run fitting procedure for each site

for iSites = 1:length(Sites)
    % set up inputs
    
    % define data
    try
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
        NumYears  = NumYears + 15;
        SwitchYear  = eval(sprintf('%sSwitchYear',Sites{iSites}));
        tEnd = [8 8 22 25 25 8 5 3]+10; % last year of MDA
        APYear  = tEnd(iSites)+1; %eval(sprintf('%sYear',Sites{iSites})); % first age profile
        APMf = eval(sprintf('%sMfAge',Sites{iSites})); % first age profile
        APCfa = eval(sprintf('%sCfaAge',Sites{iSites})); % first age profile
    catch ME
        % If the file is missing, complain but continue.
        if (strcmp(ME.identifier,'MATLAB:UndefinedFunction'))
            fprintf('\nWarning: Cannot find variables for %s\n', Sites{iSites});
            fprintf('Check PostIntv_data_SSA.m. Continuing...\n');
            continue;
            % Other error, throw it and die
        else
            rethrow(ME)
        end
    end
    
    % convert data into format compatible with model
    
    if isnan(APMf)
        ageMax = 69; % default max age in a population
    else
        ageMax = max(APMf(:,4));
    end
    % correct mf data for blood sample volume
    if ~isnan(MfData(1,1)) && (Vol == 20)
        MfData(:,3)=floor(min(MfData(:,2),MfData(:,3)*1.95));
    elseif ~isnan(MfData(1,1)) && (Vol == 100 || Vol == 60)
        MfData(:,3)=floor(min(MfData(:,2),MfData(:,3)*1.15));
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
    
    % VC coverages
    ITNCoverages=[zeros(1,(ITNYear-1)*12),ones(1,length((ITNYear-1)*12+1:NumYears*12)).*ITNCov];
    IRSCoverages=ones(1,NumYears*12)*IRSCov;
    
    % Elimination threshold
    EP_Th = 1;
    
    %% load baseline fits
    % ABR, L3Values, ParameterVectors, ageMthMax, bCulex, demog, mfPrevArray are loaded.
    % Load the baseline fits from IO. If file not found, try the next one.
    loadFile = sprintf('../IO/OUT/Hindcast%s.mat',Sites{iSites});
    try
        load(loadFile);
    catch ME
        % If file not found, try the next one
        if (strcmp(ME.identifier, 'MATLAB:load:couldNotReadFile'))
            fprintf('\nWarning: %s not found. Trying next baseline file...\n', loadFile);
            continue;
        else
            rethrow(ME) % If its not file not found, throw the error.
        end
    end
    
     % estimate worm kill parameter from data
    WormKill = [];
    for j = 1:length(MDAReg)
        if MDAReg(j) == 1
            WormKill(:,j) = 0 + 0.2*rand(length(mfPrevIntv(1,:)),1); %ones(length(mfPrevArray(1,:)),1)*0.1; % 0 + 0.2*rand(length(mfPrevArray(1,:)),1);
        elseif MDAReg(j) == 2
            WormKill(:,j) = 0.45 + 0.2*rand(length(mfPrevIntv(1,:)),1); % ones(length(mfPrevArray(1,:)),1)*0.55; % 0.45 + 0.2*rand(length(mfPrevArray(1,:)),1);
        elseif MDAReg(j) == 4
            WormKill(:,j) =  0.25 + 0.2*rand(length(mfPrevIntv(1,:)),1); % ones(length(mfPrevArray(1,:)),1)*0.35; % 0.25 + 0.2*rand(length(mfPrevArray(1,:)),1);
        end
    end
    
    % MDA Coverage
    MonthlyMDACov = zeros(NumYears*12,length(mfPrevIntv(1,:)));
    for i = 1:length(MDACov)
        MonthlyMDACov((i-1)*MDAFreq(1)+1:i*MDAFreq(1),:) = MDACov(i)*(0.75 + 0.25*rand(MDAFreq(1),length(mfPrevArray(1,:))));
    end
    
    MultiVecMBR = MBRIntv(1,:)'; %ABR/12;
    kId = 1:length(MultiVecMBR);
    ageMthMax = ageMax*12;
    da = 1;

    NParamVecs = length(kId);
    VoverH = MultiVecMBR./ParameterVectors(1,:)';
    toleranceX = 0.000001; % Used in endemic equilibrium of state variables
    [L3Values,~,~] = ...
        Calculate_EndemicEquil_L3_Mf_CFA(NParamVecs,...
        ageMthMax,bCulex,demog,ParameterVectors,VoverH,da,toleranceX);
    
    %% model interventions
    [mfAgePrevArray,cfaAgePrevArray,mfPrevIntv,cfaPrevIntv,MBRIntv,L3IntvPrev] =...
        Modelling_MDAplusVC_L3prev(kId,...
        ParameterVectors,L3Values,demog,ageMthMax,da,bCulex,...
        AgeLimits,RegimenEfficacy,MDAInterval,NumYears,...
        MultiVecMBR,IRSCoverages,ITNCoverages,MonthlyMDACov,...
        SwitchMonth,VCparams,WormKill,APYear);
    
    save(sprintf('../IO/OUT/Intv%s.mat',char(Sites{iSites})),...
        'mfPrevIntv','MBRIntv','cfaPrevIntv','L3IntvPrev','RegimenEfficacy',...
        'MonthlyMDACov','IRSCoverages',...
        'ITNCoverages','MfData','SwitchMonth','mfAgePrevArray','cfaAgePrevArray','APYear');
end

