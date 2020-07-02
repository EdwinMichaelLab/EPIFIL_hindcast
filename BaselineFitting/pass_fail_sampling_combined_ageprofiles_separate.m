function kId = pass_fail_sampling_combined_ageprofiles_separate(MfData1,MfData2,mfPrevArray,CfaData1,CfaData2,cfaPrevArray,SIR_samples)

% Initialize the index and likelihood arrays
LikArray1 = zeros(length(mfPrevArray(1,:)),2); % for likelihoods of pass/fail criteria fits
LikArray2 = zeros(length(mfPrevArray(1,:)),2); 

% calc likelihoods for mf
% calculate 95% CI upper and lower bounds of each mf data point
% if ~isnan(MfData)
MfBounds1 = get_the95LU_bounds_agedata(MfData1);

% calculate likelihood of each sampled parameter vector using pass/fail
for i=1:length(mfPrevArray(1,:))
    LikArray1(i,1)= LikArray1(i,1) + likelihood_using_ranks_by_passFail_Mf_Only(mfPrevArray(:,i),MfBounds1);
end
% end
MfBounds2 = get_the95LU_bounds_agedata(MfData2);

% calculate likelihood of each sampled parameter vector using pass/fail
for i=1:length(mfPrevArray(1,:))
    LikArray1(i,2)= LikArray1(i,2) + likelihood_using_ranks_by_passFail_Mf_Only(mfPrevArray(:,i),MfBounds2);
end

% calc likelihoods for cfa
% calculate 95% CI upper and lower bounds of each cfa data point
% if ~isnan(CfaData)
CfaBounds1 = get_the95LU_bounds_agedata(CfaData1);

% calculate likelihood of each sampled parameter vector using pass/fail
for i=1:length(mfPrevArray(1,:))
    LikArray2(i,1) = LikArray2(i,1) + likelihood_using_ranks_by_passFail_Mf_Only(cfaPrevArray(:,i),CfaBounds1);
end
% end

CfaBounds2 = get_the95LU_bounds_agedata(CfaData2);

% calculate likelihood of each sampled parameter vector using pass/fail
for i=1:length(mfPrevArray(1,:))
    LikArray2(i,2) = LikArray2(i,2) + likelihood_using_ranks_by_passFail_Mf_Only(cfaPrevArray(:,i),CfaBounds2);
end

% resample by pass/fail criteria
if length(MfData1(:,1)) == 1
    kId = find( LikArray > 0);
else
    %     kId = find( LikArray > length(MfData(:,1))-1);
%     kId = ReSampleMfPrevalence_passFail_Mf_only(LikArray,SIR_samples);
kId1 = find( LikArray1(:,1) >= length(MfData1(:,1))-2);
kId2 = find( LikArray1(:,2) >= length(MfData1(:,1))-2);
kId3 = find( LikArray2(:,1) >= length(CfaData1(:,1))-2);
kId4 = find( LikArray2(:,2) >= length(CfaData1(:,1))-2);

kIda = intersect(kId1,kId2);
kIdb = intersect(kId3,kId4);
kId =  intersect(kIda,kIdb);

end

end
