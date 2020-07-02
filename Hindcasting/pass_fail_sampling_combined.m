function kId = pass_fail_sampling_combined(MfData,mfPrevArray,CfaData,cfaPrevArray,SIR_samples)

% Initialize the index and likelihood arrays
LikArray1 = zeros(length(mfPrevArray(1,:)),1); % for likelihoods of pass/fail criteria fits
LikArray2 = zeros(length(cfaPrevArray(1,:)),1); % for likelihoods of pass/fail criteria fits

% calc likelihoods for mf
% calculate 95% CI upper and lower bounds of each mf data point
MfBounds = get_the95LU_bounds_agedata(MfData);

% calculate likelihood of each sampled parameter vector using pass/fail
for i=1:length(mfPrevArray(1,:))
    LikArray1(i)= likelihood_using_ranks_by_passFail_Mf_Only(mfPrevArray(:,i),MfBounds);
end

% calc likelihoods for cfa
% calculate 95% CI upper and lower bounds of each cfa data point
CfaBounds = get_the95LU_bounds_agedata(CfaData);
 
% calculate likelihood of each sampled parameter vector using pass/fail
for i=1:length(mfPrevArray(1,:))
     LikArray2(i) = likelihood_using_ranks_by_passFail_Mf_Only(cfaPrevArray(:,i),CfaBounds);
end


% resample by pass/fail criteria
if length(MfData(:,1)) == 1
    kId = find( LikArray > 0);
else
    %     kId = find( LikArray > length(MfData(:,1))-1);
   % kId = ReSampleMfPrevalence_passFail_Mf_only(LikArray,SIR_samples);

kId1 = find( LikArray1(:,1) >= length(MfData(:,1))-2);
kId2 = find( LikArray2(:,1) >= length(CfaData(:,1))-2);

kId =  intersect(kId1,kId2);

end

end