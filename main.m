function [results] = main()
% Author:    Niki Martinel
% Date:      2016/09/26 
% Revision:  0.1
% Copyright: Niki Martinel, 2016

%% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%   INITALIZE PARAMETERS
pars = init_parameters( '001', 'WARD', [1 2; 1 3; 2 3]);

%% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%   LOAD DATASET
dataset = load_dataset( pars );

%% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%   FEATURE EXTRACTION
features = extract_features( dataset, pars );

%% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%   GET TRAIN AND TEST SETS
[~, ~, test] = split_dataset(dataset, pars);

%% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% REID
[tests] = reidentification(features.patches, features.X, test, pars);

%% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% RESULTS
[results] = compute_results(tests, pars);


end