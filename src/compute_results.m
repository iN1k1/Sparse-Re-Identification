function [results] = compute_results( tests, pars, varargin )
% Author:    Niki Martinel
% Date:      2013/02/21 10:38:19
% Revision:  0.1
% Copyright: Niki Martinel, 2013


% Loop through all tests
for i=1:size(tests,1)
    for t=1:size(tests,2)
        [CMC{i,t}, CMCExpectation(i, :, t), nAUCCMC(i,t),...
            queryPersonsIDs{i,t}, matchingIDs{i,t}] ...
            = NM_CMC_ROC(tests(i,t).testSet.ID, tests(i,t).score(:,2), 'isDist', true);
    end
  % Average values
    
%     % Make CMC uniform across all the traials. This is needed as some
%     % persons might not appear in two cameras or and them are sometimes
%     % used for training or for classification
%     cmcLength = unique(cellfun(@length, CMC(i,:)));
%     if length(cmcLength) > 1
%         cmcLength = pars.settings.numPersons * pars.settings.learningSets(2);
%         for t=1:size(CMC,2)
%             cmci = CMC{i,t};
%             if length(cmci) < cmcLength
%                 CMC{i,t} = [cmci 100];
%             else
%                 CMC{i,t} = interp1(1:length(cmci), cmci, linspace(1, length(cmci), cmcLength), 'cubic');
%             end
%         end
%     end
    results(i).CMC = mean(vertcat(CMC{i,:}),1);
    results(i).CMCmed = median(vertcat(CMC{i,:}),1);
    
    results(i).CMCExpectation = mean(CMCExpectation(i,:,:), 3);
    results(i).nAUCCMC = mean(nAUCCMC(i,:)); 
    results(i).nAUCCMCmed = sum(results(i).CMCmed)/(100*length(results(i).CMCmed));

    % Best run
    [~, bestRunIdx] = max(nAUCCMC(i,:));
    results(i).CMCBest = CMC{i,bestRunIdx};
    results(i).nAUCCMCBest = nAUCCMC(i,bestRunIdx);
    results(i).CMCExpectationBest = CMCExpectation(i,:,bestRunIdx);
    
    % Matching results
    results(i).queryPersonsIDs = queryPersonsIDs(i,:);
    results(i).matchingIDs = matchingIDs(i,:);
 
end

end
