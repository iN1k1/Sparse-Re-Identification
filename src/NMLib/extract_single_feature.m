function [ bodyParts ] = extract_single_feature( imageRGB, mask, weightImg, featType, featPars, pars )
%EXTRACT_SINGLE_FEATURE Summary of this function goes here
%   Detailed explanation goes here

% Get images
featPars.detectTorsoAndLegs = pars.body.detectTorsoAndLegs;
featPars.splitEachBodyPart = pars.body.splitEachBodyPart;
featPars.torso = pars.body.torso;
featPars.legs = pars.body.legs;
if isempty(featPars.blockSize)
    featPars.blockSize = pars.body.blockSize;
    featPars.step = pars.body.step;
end
[imageUpperTorso, imageLowerTorso, imageUpperLegs, imageLowerLegs] = NM_body_parts(imageRGB, mask, featPars);
weightPatches = [];

% Convert each image in [0,1]
if strcmpi(featType, 'sparseHaar')
    range = NM_color_range(featPars.colorSpace);
    sz = [size(imageUpperTorso{1},1), size(imageUpperTorso{1},2)];
    minVal = cat(3, repmat(range(1,1), sz), repmat(range(2,1), sz), repmat(range(3,1), sz));
    maxVal = cat(3, repmat(range(1,2), sz), repmat(range(2,2), sz), repmat(range(3,2), sz));
    if size(imageUpperTorso{1},3) == 1
        minVal(:,:,2:3) = [];
        maxVal(:,:,2:3) = [];
    end

    imageUpperTorso = cellfun(@(x)(rescale_color_image(x, minVal, maxVal)), imageUpperTorso, 'UniformOutput', false);
    if ~isempty(imageLowerTorso)
        imageLowerTorso = cellfun(@(x)(rescale_color_image(x, minVal, maxVal)), imageLowerTorso, 'UniformOutput', false);
    end
    if ~isempty(imageUpperLegs)
        imageUpperLegs = cellfun(@(x)(rescale_color_image(x, minVal, maxVal)), imageUpperLegs, 'UniformOutput', false);
    end
    if ~isempty(imageLowerLegs)
        imageLowerLegs = cellfun(@(x)(rescale_color_image(x, minVal, maxVal)), imageLowerLegs, 'UniformOutput', false);
    end
elseif strcmpi(featType, 'whist') && ~isempty(weightImg)
    weightPatches = NM_dense_patches(weightImg , featPars.blockSize(1), featPars.blockSize(2), featPars.step(1), featPars.step(2));
end

    
%NM_pose_estimation(im, 'mixture', 'model', NM_pose_estimation_load_model('mixture'));
%bodyParts = [];
%return

% Extract features
if pars.body.detectTorsoAndLegs
    if pars.body.splitEachBodyPart
        for p=1:length(imageUpperTorso)
            upperTorso{p} = [];
            lowerTorso{p} = [];
            upperLegs{p} = [];
            lowerLegs{p} = [];
            for c=1:size(imageUpperTorso{1},3)
                upperTorso{p} = NM_extractFeatures(imageUpperTorso{p}, featType, featPars);
                lowerTorso{p} = NM_extractFeatures(imageLowerTorso{p}, featType, featPars);
                
                upperLegs{p} = NM_extractFeatures(imageUpperLegs{p}, featType, featPars);
                lowerLegs{p} = NM_extractFeatures(imageLowerLegs{p}, featType, featPars);
            end
            
        end
    else
        lowerLegs = [];
        lowerTorso = [];
        for p=1:length(imageUpperTorso)
            upperTorso{p} = [];
            upperLegs{p} = [];
            for c=1:size(imageUpperTorso{1},3)
                upperTorso{p} = NM_extractFeatures(imageUpperTorso{p}, featType, featPars);
                upperLegs{p} = NM_extractFeatures(imageUpperLegs{p}, featType, featPars);
            end
        end
    end
else
    lowerTorso = [];
    upperLegs = [];
    lowerLegs = [];
    for p=1:length(imageUpperTorso)
        if ~isempty(weightPatches)
            featPars.weight = weightPatches{p};
        end
        upperTorso{p} = NM_extractFeatures(imageUpperTorso{p}, featType, featPars);
    end
end

bodyParts{1} = upperTorso;
bodyParts{2} = lowerTorso;
bodyParts{3} = upperLegs;
bodyParts{4} = lowerLegs;

end

function [rescaledImage] = rescale_color_image(image, minVal, maxVal)
rescaledImage = (image-minVal) ./ (maxVal-minVal);
end

