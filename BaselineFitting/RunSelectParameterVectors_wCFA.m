% Function: RunSelectParameterVectors_wCFA
%
% Description: Runs through parameter vectors for baseline fitting and likelihood
%              calculation, adjusts for MF-data quality.
% Inputs:
%   Sites,NParamVecs,MfData,CfaData,ABR,...
%    SIR_samples,demoA,demoB,ageMax,bCulex,minPrev,maxPrev,ABRmax,ABRmin
%   Notes: Not all inputs will be used depending on quality of the MfData.
%
% Outputs:
%   Writes to file, 'ParamVectors%s.mat', Sites.
%   Writes: 'bCulex','demog','ABR','ParameterVectors','L3Values','ageMthMax',...
%    'mfPrevArray','cfaPrevArray','minPrev','maxPrev','ABRmax','ABRmin','MfData','CfaData'.
% ________________________________________

function [ABR,ParameterVectors,L3Values,mfPrevArray,cfaPrevArray,demog,ageMthMax] = ...
    RunSelectParameterVectors_wCFA(Sites,NParamVecs,MfData,CfaData,ABR,...
    SIR_samples,demoA,demoB,ageMax,bCulex,minPrev,maxPrev,ABRmax,ABRmin)

ageMthMax = ageMax*12; % convert max age in years to months
demog = [1, 1, demoA, demoB, 1]; % set up demographic parameters

da = 1; % integration step-size of 1 month
toleranceX = 0.000001; % Used in endemic equilibrium of state variables

% Set random number stream
RandStream.setGlobalStream...
    (RandStream('mt19937ar','seed',rand*sum(10000*clock)));  % mt19937ar:random generator

% Initialize Output Arrays
ParameterVectors = [];
L3Values1 = [];
mfPrevArray1 = [];
cfaPrevArray1 = [];
ABR1 = [];

% Loop to generate an initial sample of parameter vectors, run the model
% to equilibrium, and resample best fits until the desired number of best fits defined by SIR_samples is
% achieved. In absense of Mfdata, resample best fit until the desired number of plausible mf fits are selected
NumParam = 0;
numSamples = SIR_samples;

% Obtaining the boolean variables needed to process MfData.
% Each of the following variables represent the quality of the MfData
[AgeMf_asInput, OverallMf_asInput, OverallMfRange_asInput] = get_dataQualityVars(MfData);

if OverallMfRange_asInput % No Baseline Mf, different loop condition in this case.
    numSamples = NParamVecs; % Changing loop condition to match RSPV_OverallMfRange
end

while ( NumParam < numSamples ) % Either #SIR or NParamVecs, depending on quality of MfData
    
    % generate an initial sample of parameter vectors
    ParamVectors = ParameterVectors_LHSDesign_based(NParamVecs);
    
    % calculate VoverH (this is a vector for one species, matrix for more than one species)
    [~,n] = size(ABR); % n is the number of species
    VoverH = (ABR/12)./repmat(ParamVectors(1,:)',[1 n]); % VoverH = MBR/beta
    
    % Running to equilibrium
    [L3Values,mfPrevArray,cfaPrevArray] = ...
        Calculate_EndemicEquil_L3_Mf_CFA(NParamVecs,...
        ageMthMax,bCulex,demog,ParamVectors,VoverH,da,toleranceX);
    
    if AgeMf_asInput || OverallMf_asInput % if baseline data is available, fit to data
        if length(MfData(:,1)) > 1 && length(CfaData(:,1)) > 1
            kId = pass_fail_sampling_combined(MfData,mfPrevArray,CfaData,cfaPrevArray,SIR_samples);
        else
            if length(MfData(:,1)) == 1
                for icurve = 1:2 % loop over three types of age curves
                    
                    % calculate theoretical mf age profile data
                    MidAge = [5 14.5 24.5 34.5 44.5 54.5 64.5];
                    OverallMfPrev = MfData(:,3)/MfData(:,2);
                    APMf = getMfAgeProfile_fromLFCurves(MfData(:,2),icurve,OverallMfPrev,ageMthMax/12,demog,MidAge);
                    v=genvarname(sprintf('APMf%s',int2str(icurve)));
                    eval([v '= APMf;']);
                    
                end
                MfData = [APMf1;APMf2];
            else
                MfData = [MfData;MfData];
            end
            
            if length(CfaData(:,1)) == 1
                for icurve = 1:2 % loop over types of age curves
                    
                    % calculate theoretical cfa age profile data
                    MidAge = [5 14.5 24.5 34.5 44.5 54.5 64.5];
                    OverallCfaPrev = CfaData(:,3)/CfaData(:,2);
                    APCfa = getMfAgeProfile_fromLFCurves(CfaData(:,2),icurve,OverallCfaPrev,ageMthMax/12,demog,MidAge);
                    v=genvarname(sprintf('APCfa%s',int2str(icurve)));
                    eval([v '= APCfa;']);
                    
                end
                CfaData = [APCfa1;APCfa2];
            else
                CfaData = [CfaData;CfaData];
            end
            APMf1 = MfData(1:length(MfData(:,1))/2,:);
            APMf2 = MfData(length(MfData(:,1))/2+1:end,:);
            APCfa1 = CfaData(1:length(CfaData(:,1))/2,:);
            APCfa2 = CfaData(length(CfaData(:,1))/2+1:end,:);
            kId = pass_fail_sampling_combined_ageprofiles_separate(APMf1,APMf2,mfPrevArray,APCfa1,APCfa2,cfaPrevArray,SIR_samples);
        end
              
    elseif OverallMfRange_asInput % if no baseline data, use defined range for fitting
        for i = 1:length(mfPrevArray(1,:))  %% loop through 200,000 parameters
            OverallMf(i) = getOverallPrev_fromAgeCurves(100*mfPrevArray(:,i),...
                demog,da,ageMthMax);
        end
        kId0 = find(OverallMf > minPrev & OverallMf < maxPrev);
        
        % sample prev uniformly
        [N,edges] = histcounts(OverallMf(1,kId0)*100,int8(minPrev*100):5:int8(maxPrev*100)); % number in each percent bin
        
        kId=[];
        for i = 1:length(edges)-1
            Id = find(100*OverallMf(1,kId0)>edges(i) & 100*OverallMf(1,kId0) < edges(i+1));
            kId = [kId,kId0(Id(1:min(N)))];
        end
        
    end
    
    % Store parameters, equilibrium L3, and an equilibrium mf values
    % corresponding to the resampled parameter vectors
    ParameterVectors(:,NumParam+1:NumParam+length(kId)) = ...
        ParamVectors(:,kId);
    
    L3Values1(NumParam+1:NumParam+length(kId),:) = ...
        L3Values(kId,:);
    
    mfPrevArray1(:,NumParam+1:NumParam+length(kId)) = ...
        mfPrevArray(:,kId);
    
    cfaPrevArray1(:,NumParam+1:NumParam+length(kId)) = ...
        cfaPrevArray(:,kId);
    
    ABR1(NumParam+1:NumParam+length(kId),:) = ...
        12*VoverH(kId,:).*ParamVectors(1,kId)';
    
    
    % Update number of parameter vectors selected
    NumParam = NumParam + length(kId);
    fprintf(1,'Num selected kIds = %d, Total Num Params = %d\n',...
        length(kId),NumParam);
end

L3Values = L3Values1;
mfPrevArray = mfPrevArray1;
cfaPrevArray = cfaPrevArray1;
ABR = ABR1;

% save sampled best fitting outputs and basic descriptive parameters of
% the site
save(sprintf('../IO/OUT/ParamVectors%s.mat',Sites),'bCulex',...
    'demog','ABR','ParameterVectors','L3Values','ageMthMax',...
    'mfPrevArray','cfaPrevArray','minPrev','maxPrev','ABRmax','ABRmin','MfData','CfaData');

end
