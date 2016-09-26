function [ CMC, CMCExpectation, nAUCCMC, probePersonsIDs, matchingIDs ] = ...
    NM_CMC_ROC( personsIDs, similarityScore, varargin)
%NM_CMC Summary of this function goes here
%   Detailed explanation goes here

p = inputParser;
p.addOptional('isDist', false);
p.addOptional('statsFunHandle', @mean);
p.parse(varargin{:});
opts = p.Results;

% Fix old score size
if size(similarityScore,2) > 1
    similarityScore = similarityScore(:,2);
end

% Query persons IDs (using unique the IDs are in increasing order)
probePersonsIDs = unique(personsIDs(:,1), 'stable');
galleryPersonsIDs = unique(personsIDs(:,2), 'stable');
probeIDs = personsIDs(:,1);
galleryIDs = personsIDs(:,2);
    
% Matching IDs
matchingIDs = zeros(length(probePersonsIDs), length(galleryPersonsIDs));
sortedTestedPersonIDs = zeros(1, length(galleryPersonsIDs));

% Loop through all the query images%
match = zeros(1,length(galleryPersonsIDs));

%aa = num2cell(personsIDs,2);
%cellfun(@(x,score)(score(x), aa, repmat({similarityScore}, size(aa)), 'UniformOutput', false)

% Loop through all the probe persons
%par
for k=1:length(probePersonsIDs)

    % Take the similarities only of pairs where current probe ID is appearing
    idxk = probeIDs == probePersonsIDs(k);
    simsk  = similarityScore(idxk);
    
    % Take the corresponding gallery IDs (unique ones as well)
    idb = galleryIDs(idxk);
    galleryPersonsIDsIterK = unique(idb, 'stable');
    
    % Compute the score between the probe and every other gallery 
    avgScore = arrayfun(@(idbb)( opts.statsFunHandle(simsk(idb==idbb)) ), galleryPersonsIDsIterK);
    
    % Some negative samples might not been tested since they are
    % randomly chosen
    %if ~opts.excludeNaN
        if opts.isDist
            avgScore(isnan(avgScore)) = inf;
        	[~, sortedIdx] =  sort(avgScore, 'ascend');
        else
            avgScore(isnan(avgScore)) = -inf;
            [~, sortedIdx] =  sort(avgScore, 'descend');
        end
        sortedTestedPersonIDs = galleryPersonsIDsIterK(sortedIdx);
        posMatch = probePersonsIDs(k) == sortedTestedPersonIDs;
        match(posMatch) = match(posMatch) + 1;

        % Matching persons
        % We should use the 1:length trick because for some person we can
        % have that the match is computed with a subset of persons only..
        matchingIDs(k,1:length(sortedTestedPersonIDs)) = sortedTestedPersonIDs;
end

% CMC and normalized Area Under Curve (of CMC)
CMC = cumsum(match);
CMC = 100*(CMC / max(CMC));
nAUCCMC = sum(CMC)./(max(CMC)*length(match)); 

% Compute CMC expectation
CMCExpectation = NM_expectation_from_CMC(CMC);
end

