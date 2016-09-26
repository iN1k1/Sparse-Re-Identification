function [pars] = init_parameters(testID, datasetName, testCams)

fprintf('Initialize parameters...');
t = tic;

%% ========================================================================
%   DATASET
% =========================================================================
pars.dataset.name = datasetName;
pars.dataset.imageWidth = 64;
pars.dataset.imageHeight = 128;
pars.dataset.useMasks = false;
pars.dataset.imageColorSpace = 'RGB';

%% ========================================================================
%   SETTINGS
% =========================================================================
pars.settings.testID = testID;

% Output file on which save test data
pars.settings.outputPrefix = [pars.dataset.name, '_Id', pars.settings.testID];
pars.settings.rootFolder = pwd;
pars.settings.dataFolder = fullfile(pars.settings.rootFolder, 'test', pars.settings.outputPrefix);
if ~exist(fullfile(pars.settings.dataFolder), 'file')
    mkdir(fullfile(pars.settings.dataFolder));
end

pars.settings.testCams = testCams;
pars.settings.numTests = 10;
pars.settings.numPersons = [];
pars.settings.numberImages = 1;
pars.settings.useSameImageIndexForPositiveAndNegativeSamples = true;
pars.settings.numberSamplesPerPersonTraining = [pars.settings.numberImages pars.settings.numberImages];
pars.settings.numberSamplesPerPersonTesting = [-1 -1];
pars.settings.trainAndTestWithSamePersons = false;
pars.settings.testPeopleIDs = [];
pars.settings.learningSets = [0 1];
pars.settings.extendNumberOfImagesPerPersonPerCamera = 0;
pars.settings.availableColorSpaces = {'grayscale', 'HSV', 'HSL', 'HSI', 'RGB', 'LCH', 'Lab', 'Luv', 'YCbCr', 'YPbPr', 'XYZ', 'YUV'};



%% ========================================================================
%   BODY PARTS
% =========================================================================
pars.body.detectTorsoAndLegs = false;
pars.body.splitEachBodyPart = false;
pars.body.torso = [64 64 0 0];
pars.body.legs  = [64 64 0 0];
pars.body.blockSize = [16 16];
pars.body.step = [8 8];
pars.body.mag   = 1;


%% ========================================================================
%   TRACKING REGIONS
% =========================================================================
pars.regions.size = [16 16];
se = strel('disk', floor(pars.regions.size(1)*0.5),0);
pars.regions.mask = se.getnhood;
pars.regions.num = 24;
pars.regions.minNeighDist = 1;
pars.regions.maxNeighDist = 3;
pars.regions.numNeighs = 4;
pars.regions.minClassifierConf = 0.8;

%% ========================================================================
%   FEATURES
% =========================================================================
histEnabled = true;
histBins = [24 24 24];
histNormalize = true;

lbpEnabled = true;
lbpBlockSize = [];
lbpStep = [];

%----------------------------------------------------------------------
% Color Histogram

pars.features.HSV.enabled = histEnabled;
pars.features.HSV.colorSpace = 'HSV';            % Lab, LCH, HSV, HSI, RGB, XYZ
pars.features.HSV.bins = histBins;
pars.features.HSV.weights = [];
pars.features.HSV.excludeRange = [-inf -1];
pars.features.HSV.normalize = histNormalize;

pars.features.Lab.enabled = histEnabled;
pars.features.Lab.colorSpace = 'Lab';            % Lab, LCH, HSV, HSI, RGB, XYZ
pars.features.Lab.bins = histBins;
pars.features.Lab.weights = [];
pars.features.Lab.excludeRange = [-inf -129];
pars.features.Lab.normalize = histNormalize;

pars.features.YUV.enabled = histEnabled;
pars.features.YUV.colorSpace = 'YUV';
pars.features.YUV.bins = histBins;
pars.features.YUV.excludeRange = [-inf -100];
pars.features.YUV.normalize = histNormalize;

%----------------------------------------------------------------------
% Local Binary Patterns
pars.features.lbp.enabled = lbpEnabled;
pars.features.lbp.colorSpace = 'grayscale';
pars.features.lbp.patchSize = [];
pars.features.lbp.step = [];
pars.features.lbp.normalizedHistogram = false;
pars.features.lbp.points = 8;
pars.features.lbp.radius = 1;
pars.features.lbp.mapping = 'u2';
pars.features.lbp.blockSize = lbpBlockSize;
pars.features.lbp.step = lbpStep;

fprintf('done in %.2f(s)\n', toc(t));


end