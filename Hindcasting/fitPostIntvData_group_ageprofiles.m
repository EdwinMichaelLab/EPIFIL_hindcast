function [mfPrevIntv,L3Intv,MBRIntv,cfaPrevIntv,ParameterVectors,...
    mfPrevArray,cfaPrevArray,cfaAgePrevArray,mfAgePrevArray,APMf,APCfa,WormKill] = ...
    fitPostIntvData_group_ageprofiles(MfData,CfaData,APMf,mfAgePrevArray,demog,da,...
    ageMthMax,SIR_samples,mfPrevIntv,L3Intv,MBRIntv,cfaPrevIntv,...
    ParameterVectors,mfPrevArray,cfaPrevArray,cfaAgePrevArray,APCfa,APYear,...
    WormKill)

%% fit to first age profiles if available

% fit to MfData
if ~isnan(APMf(1,1)) && ~isnan(APCfa(1,1))
    
    kId = pass_fail_sampling_combined(APMf,mfAgePrevArray/100,APCfa,cfaAgePrevArray/100,SIR_samples);
    fprintf(1,'combined data: n = %d\n',length(kId));
    
    % save fits
    if ~isempty(kId)
        L3Intv = L3Intv(:,kId);
        MBRIntv = MBRIntv(:,kId);
        mfPrevIntv = mfPrevIntv(:,kId);
        cfaPrevIntv = cfaPrevIntv(:,kId);
        ParameterVectors = ParameterVectors(:,kId);
        mfPrevArray = mfPrevArray(:,kId);
        cfaPrevArray = cfaPrevArray(:,kId);
        mfAgePrevArray = mfAgePrevArray(:,kId);
        cfaAgePrevArray = cfaAgePrevArray(:,kId);
        WormKill = WormKill(kId,:);
    end
else
    
    idxMf = find(MfData(:,1) == (APYear-1)*12+1);
    idxCfa = find(CfaData(:,1) == (APYear-1)*12+1);
    for icurve = 1:2 % loop over three types of age curves
        
        % calculate theoretical mf age profile data
        MidAge = [5 14.5 24.5 34.5 44.5 54.5 64.5];
        OverallMfPrev = MfData(idxMf,3)/MfData(idxMf,2);
        APMf = getMfAgeProfile_fromLFCurves(MfData(idxMf,2),icurve,OverallMfPrev,ageMthMax/12,demog,MidAge);
        v=genvarname(sprintf('APMf%s',int2str(icurve)));
        eval([v '= APMf;']);
        
        % calculate theoretical cfa age profile data
        MidAge = [5 14.5 24.5 34.5 44.5 54.5 64.5];
        OverallCfaPrev = CfaData(idxCfa,3)/CfaData(idxCfa,2);
        APCfa = getMfAgeProfile_fromLFCurves(CfaData(idxCfa,2),icurve,OverallCfaPrev,ageMthMax/12,demog,MidAge);
        v=genvarname(sprintf('APCfa%s',int2str(icurve)));
        eval([v '= APCfa;']);
        
    end
    APMf = [APMf1;APMf2];
    APCfa = [APCfa1;APCfa2];
    kId = pass_fail_sampling_combined_ageprofiles_separate(APMf1,APMf2,mfAgePrevArray/100,APCfa1,APCfa2,cfaAgePrevArray/100,SIR_samples);
    
    fprintf(1,'combined data: n = %d\n',length(kId));
    
    % save fits
    if ~isempty(kId)
        L3Intv = L3Intv(:,kId);
        MBRIntv = MBRIntv(:,kId);
        mfPrevIntv = mfPrevIntv(:,kId);
        cfaPrevIntv = cfaPrevIntv(:,kId);
        ParameterVectors = ParameterVectors(:,kId);
        mfPrevArray = mfPrevArray(:,kId);
        cfaPrevArray = cfaPrevArray(:,kId);
        mfAgePrevArray = mfAgePrevArray(:,kId);
        cfaAgePrevArray = cfaAgePrevArray(:,kId);
        WormKill = WormKill(kId,:);
    end
end

end
