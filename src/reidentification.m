function [results] = reidentification(patches, features, testSplit, pars)

% Try to load predefined splits
fprintf('Performing re-identification...');
tt = tic;

% Try to load data
reidFile = fullfile(pars.settings.dataFolder, [pars.settings.outputPrefix '_reid.mat']);
if exist(reidFile, 'file')
    load(reidFile);
else
    % Prior
    prior = 3.5;
    
    % Threshold
    threshold = 0.9;
    
    % Loop over all camera pairs
    for c=1:size(pars.settings.testCams,1)
        
        % Loop over trials
        for t=1:size(testSplit,2)
            
            % Get split Idx
            ida = testSplit(c,t).index(:,1);
            idb = testSplit(c,t).index(:,2);
            
            % Init score 
            score = zeros(length(ida),2);
            
            % Get features
            featA = features(ida,:,:);
            featB = features(idb,:,:);
            
            % Get patches
            patchesA = patches(ida,:,:,:,:);
            patchesB = patches(idb,:,:,:,:);
            
            % Loop over all pairs
            parfor pair=1:length(ida)
                
                % Features for single pair
                fA = squeeze(featA(pair,:,:));
                fB = squeeze(featB(pair,:,:));
                
                % Patches for single pair
                pA = reshape(squeeze(patchesA(pair,:,:,:,:)), [], size(patchesA,5));
                pB = reshape(squeeze(patchesB(pair,:,:,:,:)), [], size(patchesB,5));
                
                % Select useful patches
                w = feature_selection_weights(pA',pB',threshold);
                
                % Compute feature dissimilarities
                d = NM_pdist(fA,fB, 'chisq');
                
                % Apply patch weights to feature dissimilarities
                w(isinf(w))=0;
                dw = d .* (w+prior*eye(size(d)));
                dd = mean(dw(dw~=0));
                if isnan(dd), dd = inf; end
                score(pair,2) = dd;
            end
            
            % Update scores in final structure
            results(c,t).score = score;
            results(c,t).testSet = testSplit(c,t);
            
        end
        
    end
    
    
    % Save data
    try
        save(reidFile, 'results');
    catch ME
        warning('reidentification:save', 'Unable to save reidentification data on file %s.', reidFile)
    end
end

fprintf('done in %.2f(s)\n', toc(tt));


end
