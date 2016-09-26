function [train, cv, test] = split_dataset( dataset, pars )

fprintf('Splitting dataset...');
tt = tic;

% Try to load data
splitsFile = fullfile(pars.settings.dataFolder, [pars.settings.outputPrefix '_split.mat']);
if exist(splitsFile, 'file')
    load(splitsFile);
else
    % Loop through all transformations
    for i=1:size(pars.settings.testCams,1)
        
        % Loop through all tests
        for t=1:pars.settings.numTests
            
            %% Split into train/test datasets
            allPeopleIDs = 1:dataset.peopleCount;
            
            % Randomly shuffled person IDs
            if ~isempty(pars.settings.testPeopleIDs)
                if size(pars.settings.testPeopleIDs, 2) > 1
                    allPeopleIDs = pars.settings.testPeopleIDs(:,t);
                else
                    allPeopleIDs = pars.settings.testPeopleIDs;
                end
            else
                allPeopleIDs = allPeopleIDs(randperm(length(allPeopleIDs)));
            end
            
            % Reduce dataset size if requried
            if ~isempty(pars.settings.numPersons)
                allPeopleIDs = allPeopleIDs(1:pars.settings.numPersons);
            end
            
            % Number of person per each set (we can have some overlapping)
            numPersonTrain = round(pars.settings.learningSets(1)*length(allPeopleIDs));
            numPersonTest = round(pars.settings.learningSets(2)*length(allPeopleIDs));
            if numPersonTest+numPersonTrain > length(allPeopleIDs)
                numPersonTest = length(allPeopleIDs)-numPersonTrain;
            end
            
            % Person IDs used for training and testing
            idPersonTrain = sort(allPeopleIDs(1:numPersonTrain), 'ascend');
            idPersonTest = sort(allPeopleIDs(numPersonTrain+1:numPersonTrain+numPersonTest), 'ascend');
            
            % Use the same persons for training and testing
            if pars.settings.trainAndTestWithSamePersons
                idPersonTrain = idPersonTest;
            end
            
            if ~isfield(pars.settings, 'useSameImageIndexForPositiveAndNegativeSamples')
                pars.settings.useSameImageIndexForPositiveAndNegativeSamples = false;
            end
            useSameImageIndexForPositiveAndNegativeSamples = pars.settings.useSameImageIndexForPositiveAndNegativeSamples;
            
            % Get training and test samples
            if ~pars.settings.trainAndTestWithSamePersons
                
                % Training set
                train(i, t) = getSamples(dataset, idPersonTrain, pars.settings.testCams(i,:), ...
                    pars.settings.numberSamplesPerPersonTraining(1), pars.settings.numberSamplesPerPersonTraining(2), ...
                    useSameImageIndexForPositiveAndNegativeSamples);
                
                % Cross validation set
                if isfield(pars.settings, 'kfold') && ~isempty(pars.settings.kfold) && pars.settings.kfold > 0
                    cvPartition = cvpartition(length(idPersonTrain), 'kfold', pars.settings.kfold);
                    for c=1:pars.settings.kfold
                        cv.train(i,t,c) = getSamples(dataset, idPersonTrain(cvPartition.training(c)), pars.settings.testCams(i,:), pars.settings.numberSamplesPerPersonTraining(1), pars.settings.numberSamplesPerPersonTraining(2), useSameImageIndexForPositiveAndNegativeSamples);
                        cv.test(i,t,c) = getSamples(dataset, idPersonTrain(cvPartition.test(c)), pars.settings.testCams(i,:), pars.settings.numberSamplesPerPersonTesting(1), pars.settings.numberSamplesPerPersonTesting(2), useSameImageIndexForPositiveAndNegativeSamples);
                    end
                else
                    cv = [];
                end
                
                % Test set
                test(i, t) = getSamples(dataset, idPersonTest, pars.settings.testCams(i,:), ...
                    pars.settings.numberSamplesPerPersonTesting(1), pars.settings.numberSamplesPerPersonTesting(2), ...
                    useSameImageIndexForPositiveAndNegativeSamples);
            else
                
                % Training and test set
                [train(i, t), test(i,t)] = getOverlappingSamples(dataset, idPersonTrain, idPersonTest, ...
                    pars.settings.testCams(i,:), ...
                    pars.settings.numberSamplesPerPersonTraining, pars.settings.numberSamplesPerPersonTesting, ...
                    useSameImageIndexForPositiveAndNegativeSamples);
                
                % Cross validation set
                if ~isempty(pars.settings.kfold) && pars.settings.kfold > 0
                    cvPartition = cvpartition(length(idPersonTrain), 'kfold', pars.classifier.kfold);
                    for c=1:pars.settings.kfold
                        [cv.train(i,t,c), cv.test(i,t,c)] = getOverlappingSamples(dataset, idPersonTrain(cvPartition.training(c)), idPersonTrain(cvPartition.test(c)), ...
                            pars.settings.testCams(i,:), ...
                            pars.settings.numberSamplesPerPersonTraining, pars.settings.numberSamplesPerPersonTraining, ...
                            useSameImageIndexForPositiveAndNegativeSamples);
                    end
                else
                    cv = [];
                end
            end
        end
    end
    
    % Save data
    try
        save(splitsFile, 'train', 'test', 'cv');
    catch ME
        warning('split_dataset:save', 'Unable to save splits data on file %s.', splitsFile)
    end
end

% Features extraction time
fprintf('done in %.2f(s)\n', toc(tt));


end

%% GET TRAINING SAMPLES SAMPLES
function [samples] = getSamples(dataset, personsIDs, cameraPair, numPositive, numNegative, useSameImageIndexForPositiveAndNegativeSamples)

samples.ID = [];
samples.index = [];
samples.label = [];
if (numPositive < 0 && numNegative < 0)
    numPositive = abs(numPositive);
    numNegative = abs(numNegative);
    
    if strcmpi(dataset.name, 'VIPeR')
        idx1 = dataset.imageIndex(dataset.cam==cameraPair(1));
        idx2 = dataset.imageIndex(dataset.cam==cameraPair(2));
        
        idx1 = idx1(personsIDs);
        randSortedPersonsIDs = personsIDs(randperm(length(personsIDs)));
        idx2 = idx2(randSortedPersonsIDs);
        samples.index = allcomb(idx1, idx2);
        samples.ID = allcomb(personsIDs, randSortedPersonsIDs);
    else
        
        idx1_pos = cell(length(personsIDs), 1);
        idx1_neg = cell(length(personsIDs), 1);
        idx2_pos = cell(length(personsIDs), 1);
        idx2_neg = cell(length(personsIDs), 1);
        %par
        for i=1:length(personsIDs)
            idx_tmp = dataset.imageIndex(dataset.cam==cameraPair(1) & dataset.personID==personsIDs(i));
            if ~isempty(idx_tmp)
                idx_tmp = idx_tmp(randperm(length(idx_tmp)))';
                if numPositive > length(idx_tmp)
                    idx1_pos{i} = idx_tmp;
                else
                    idx1_pos{i} = idx_tmp(1:numPositive);
                end
                if numNegative > length(idx_tmp)
                    idx1_neg{i} = idx_tmp;
                else
                    idx1_neg{i} = idx_tmp(1:numNegative);
                end
            else
                idx1_pos{i} = [];
                idx1_neg{i} = [];
            end
            
            idx_tmp = dataset.imageIndex(dataset.cam==cameraPair(2) & dataset.personID==personsIDs(i));
            if ~isempty(idx_tmp)
                idx_tmp = idx_tmp(randperm(length(idx_tmp)))';
                if numPositive > length(idx_tmp)
                    idx2_pos{i} = idx_tmp;
                else
                    idx2_pos{i} = idx_tmp(1:numPositive);
                end
                if numNegative > length(idx_tmp)
                    idx2_neg{i} = idx_tmp;
                else
                    idx2_neg{i} = idx_tmp(1:numNegative);
                end
            else
                idx2_pos{i} = [];
                idx2_neg{i} = [];
            end
        end
        
        id_pos = cell2mat(arrayfun(@(x)(allcomb(repmat(personsIDs(x), numel(idx1_pos{x}), 1), repmat(personsIDs(x), numel(idx2_pos{x}), 1))),1:length(personsIDs), 'UniformOutput', false)');
        idx_pos = cell2mat(arrayfun(@(x)(allcomb(idx1_pos{x}, idx2_pos{x})), 1:length(personsIDs), 'UniformOutput', false)');
        samples.ID = [samples.ID; id_pos];
        samples.index = [samples.index; idx_pos];
        
        id_1 = cell2mat(arrayfun(@(x)(repmat(personsIDs(x), numel(idx1_neg{x}), 1)), 1:length(personsIDs), 'UniformOutput', false)');
        id_2 = cell2mat(arrayfun(@(x)(repmat(personsIDs(x), numel(idx2_neg{x}), 1)), 1:length(personsIDs), 'UniformOutput', false)');
        id_neg = allcomb(id_1, id_2);
        idx_neg = allcomb(cell2mat(idx1_neg), cell2mat(idx2_neg));
        toRemove = id_neg(:,1) == id_neg(:,2);
        id_neg(toRemove,:) = [];
        idx_neg(toRemove,:) = [];
        samples.ID = [samples.ID; id_neg];
        samples.index = [samples.index; idx_neg];
        
    end
   
else
    for i=1:length(personsIDs)
        idx = dataset.imageIndex(dataset.cam==cameraPair(1) & dataset.personID==personsIDs(i));
        idxPos = dataset.imageIndex(dataset.cam==cameraPair(2) & dataset.personID==personsIDs(i));
        idxNeg = dataset.imageIndex(dataset.cam==cameraPair(2) & ismember(dataset.personID, personsIDs(personsIDs~=personsIDs(i))));
        
        % No images of this person in camera 1
        if isempty(idx)
            continue;
        end
        
        % Positive combinations
        if ~isempty(idxPos)
            samples.ID = [samples.ID; repmat([personsIDs(i) personsIDs(i)], numPositive, 1)];
            combinations = allcomb(idx, idxPos);
            combinations = combinations(randperm(size(combinations, 1)),:);
            if numPositive > size(combinations,1)
                numPositive = size(combinations,1);
            end
            samples.index = [samples.index; combinations(1:numPositive,:)];
            %samples.label = [samples.label; ones(numPositive,1)];
        end
        
        % Negative combinations
        if ~isempty(idxNeg)
            combinations = allcomb(idx, idxNeg);
            combinations = combinations(randperm(size(combinations, 1)),:);
            if numNegative>length(combinations)
                numNegative = length(combinations);
            end
            samples.index = [samples.index; combinations(1:numNegative,:)];
            samples.ID = [samples.ID; [repmat(personsIDs(i), numNegative, 1) dataset.personID(combinations(1:numNegative,2))']];
            %samples.label = [samples.label; zeros(numNegative,1)];
        end
    end
end
if ~isempty(samples.ID)
    samples.label = samples.ID(:,1) == samples.ID(:,2);
end
end



function [samplesTrain, samplesTest] = getOverlappingSamples(dataset, personsIDsTrain, personsIDsTest, ...
    cameraPair, numSamplesTrain, numSamplesTest, useSameImageIndexForPositiveAndNegativeSamples)
samplesTrain.ID = [];
samplesTrain.index = [];
samplesTrain.label = [];
samplesTest.ID = [];
samplesTest.index = [];
samplesTest.label = [];
allPersonIDs = unique([personsIDsTrain personsIDsTest]);
if (numSamplesTrain(1) < 0 && numSamplesTest(1) < 0)
    numSamplesTrain = abs(numSamplesTrain);
    numSamplesTest = abs(numSamplesTest);
    
    idx1_tr = cell(length(allPersonIDs), 1);
    idx2_tr = cell(length(allPersonIDs), 1);
    idx1_te = cell(length(allPersonIDs), 1);
    idx2_te = cell(length(allPersonIDs), 1);
    for i=1:length(allPersonIDs)
        tmp = dataset.imageIndex(dataset.cam==cameraPair(1) & dataset.personID==allPersonIDs(i));
        tmp = tmp(randperm(length(tmp)));
        idx1_tr{i} = tmp(1:numSamplesTrain(1));
        idx1_te{i} = tmp(numSamplesTrain(1)+1:numSamplesTrain(1)+numSamplesTest(1));
        
        tmp = dataset.imageIndex(dataset.cam==cameraPair(2) & dataset.personID==allPersonIDs(i));
        tmp = tmp(randperm(length(tmp)));
        idx2_tr{i} = tmp(1:numSamplesTrain(1));
        idx2_te{i} = tmp(numSamplesTrain(1)+1:numSamplesTrain(1)+numSamplesTest(1));
    end
    
    for i=1:length(allPersonIDs)
        for j=1:length(allPersonIDs)
            
            % Training set
            if all(ismember([allPersonIDs(i) allPersonIDs(j)], personsIDsTrain))
                combinations = allcomb(idx1_tr{i}, idx2_tr{j});
                samplesTrain.ID = [samplesTrain.ID; repmat([allPersonIDs(i) allPersonIDs(j)], size(combinations, 1), 1)];
                samplesTrain.index = [samplesTrain.index; combinations];
                if i==j
                    samplesTrain.label = [samplesTrain.label; ones(size(combinations,1),1)];
                else
                    samplesTrain.label = [samplesTrain.label; zeros(size(combinations,1),1)];
                end
            end
            
            % Test set
            if all(ismember([allPersonIDs(i) allPersonIDs(j)], personsIDsTest))
                combinations = allcomb(idx1_te{i}, idx2_te{j});
                samplesTest.ID = [samplesTest.ID; repmat([allPersonIDs(i) allPersonIDs(j)], size(combinations, 1), 1)];
                samplesTest.index = [samplesTest.index; combinations];
                if i==j
                    samplesTest.label = [samplesTest.label; ones(size(combinations,1),1)];
                else
                    samplesTest.label = [samplesTest.label; zeros(size(combinations,1),1)];
                end
            end
        end
    end
   
elseif (numSamplesTrain(1) > 0 && numSamplesTest(1) < 0)
    numSamplesTrain = abs(numSamplesTrain);
    numSamplesTest = abs(numSamplesTest);
    %personID_index = [];
    
    idx1_tr = cell(length(allPersonIDs), 1);
    idx2_tr = cell(length(allPersonIDs), 1);
    idx1_te = cell(length(allPersonIDs), 1);
    idx2_te = cell(length(allPersonIDs), 1);
    for i=1:length(allPersonIDs)
        tmp = dataset.imageIndex(dataset.cam==cameraPair(1) & dataset.personID==allPersonIDs(i));
        tmp = tmp(randperm(length(tmp)));
        idx1_tr{i} = tmp(1:numSamplesTrain(1));
        idx1_te{i} = tmp(numSamplesTrain(1)+1:numSamplesTrain(1)+numSamplesTest(1));
        
        tmp = dataset.imageIndex(dataset.cam==cameraPair(2) & dataset.personID==allPersonIDs(i));
        tmp = tmp(randperm(length(tmp)));
        idx2_tr{i} = tmp(1:numSamplesTrain(1));
        idx2_te{i} = tmp(numSamplesTrain(1)+1:numSamplesTrain(1)+numSamplesTest(1));
    end
    
    for i=1:length(allPersonIDs)
        
        % Positivive training set samples
        if ismember(allPersonIDs(i), personsIDsTrain)
            combinations = allcomb(idx1_tr{i}, idx2_tr{i});
            combinations(numSamplesTrain(1)+1:end,:) = [];
            samplesTrain.ID = [samplesTrain.ID; repmat([allPersonIDs(i) allPersonIDs(i)], size(combinations, 1), 1)];
            samplesTrain.index = [samplesTrain.index; combinations];
            samplesTrain.label = [samplesTrain.label; ones(size(combinations,1),1)];
        end
        
        jjs = randperm(length(allPersonIDs));
        randAllPersonIDs = allPersonIDs(jjs);
        numAddedNegativeTrainingSamples = 0;
        for j=1:length(allPersonIDs)
            
            %Negative Training set
            if allPersonIDs(i) ~= randAllPersonIDs(j) && numAddedNegativeTrainingSamples < numSamplesTrain(2) && all(ismember([allPersonIDs(i) randAllPersonIDs(j)], personsIDsTrain))
                jj = jjs(j);
                combinations = allcomb(idx1_tr{i}, idx2_tr{jj});
                combinations(2:end,:) = [];
                samplesTrain.ID = [samplesTrain.ID; repmat([allPersonIDs(i) randAllPersonIDs(j)], size(combinations, 1), 1)];
                samplesTrain.index = [samplesTrain.index; combinations];
                samplesTrain.label = [samplesTrain.label; zeros(size(combinations,1),1)];
                numAddedNegativeTrainingSamples = numAddedNegativeTrainingSamples + 1;
            end
            
            %Positive and Negative Test set
            if all(ismember([allPersonIDs(i) allPersonIDs(j)], personsIDsTest))
                combinations = allcomb(idx1_te{i}, idx2_te{j});
                samplesTest.ID = [samplesTest.ID; repmat([allPersonIDs(i) allPersonIDs(j)], size(combinations, 1), 1)];
                samplesTest.index = [samplesTest.index; combinations];
                if i==j
                    samplesTest.label = [samplesTest.label; ones(size(combinations,1),1)];
                else
                    samplesTest.label = [samplesTest.label; zeros(size(combinations,1),1)];
                end
            end
        end
    end
else
    numSamplesTrain = abs(numSamplesTrain);
    numSamplesTest = abs(numSamplesTest);
    for i=1:length(allPersonIDs)
        idx = dataset.imageIndex(dataset.cam==cameraPair(1) & dataset.personID==allPersonIDs(i));
        idxPos = dataset.imageIndex(dataset.cam==cameraPair(2) & dataset.personID==allPersonIDs(i));
        idxNeg = dataset.imageIndex(dataset.cam==cameraPair(2) & ismember(dataset.personID, allPersonIDs(allPersonIDs~=allPersonIDs(i))));
        
        % Positive combinations
        if ismember(allPersonIDs(i), personsIDsTrain)
            samplesTrain.ID = [samplesTrain.ID; repmat([allPersonIDs(i) allPersonIDs(i)], numSamplesTrain(1), 1)];
        end
        if ismember(allPersonIDs(i), personsIDsTest)
            samplesTest.ID = [samplesTest.ID; repmat([allPersonIDs(i) allPersonIDs(i)], numSamplesTest(1), 1)];
        end
        combinations = allcomb(idx, idxPos);
        combinations = combinations(randperm(size(combinations, 1)),:);
        
        % Chose combinations for training and testing
        if ismember(allPersonIDs(i), personsIDsTrain)
            samplesTrain.index = [samplesTrain.index; combinations(1:numSamplesTrain(1),:)];
            samplesTrain.label = [samplesTrain.label; ones(numSamplesTrain(1),1)];
        end
        if ismember(allPersonIDs(i), personsIDsTest)
            samplesTest.index = [samplesTest.index; combinations(numSamplesTrain(1)+1:numSamplesTrain(1)+numSamplesTest(1),:)];
            samplesTest.label = [samplesTest.label; ones(numSamplesTest(1),1)];
        end
        
        % Negative combinations
        combinations = allcomb(idx, idxNeg);
        combinations = combinations(randperm(size(combinations, 1)),:);
        %if numNegative>length(combinations)
        %    numNegative = length(combinations);
        %end
        
        % Chose combinations for training and testing
        if ismember(allPersonIDs(i), personsIDsTrain)
            samplesTrain.index = [samplesTrain.index; combinations(1:numSamplesTrain(2),:)];
            samplesTrain.ID = [samplesTrain.ID; [repmat(allPersonIDs(i), numSamplesTrain(2), 1) dataset.personID(combinations(1:numSamplesTrain(2),2))']];
            samplesTrain.label = [samplesTrain.label; zeros(numSamplesTrain(2),1)];
        end
        if ismember(allPersonIDs(i), personsIDsTest)
            samplesTest.index = [samplesTest.index; combinations(numSamplesTrain(2)+1:numSamplesTrain(2)+numSamplesTest(2),:)];
            samplesTest.ID = [samplesTest.ID; [repmat(allPersonIDs(i), numSamplesTest(2), 1) dataset.personID(combinations(numSamplesTrain(2)+1:numSamplesTrain(2)+numSamplesTest(2),2))']];
            samplesTest.label = [samplesTest.label; zeros(numSamplesTest(2),1)];
        end
    end
end

% Ensure labels are logical values
samplesTrain.label = logical(samplesTrain.label);
samplesTest.label = logical(samplesTest.label);

end