function [ features ] = extract_features( dataset, pars )

fprintf('Extracting features...');
t = tic;

% Try to load data
featuresFile = fullfile(pars.settings.dataFolder, [pars.settings.outputPrefix '_features.mat']);
if exist(featuresFile, 'file')
    load(featuresFile);
else
    
    % Patch size and #
    regionSize = pars.regions.size; % AS [H W]
    numRegions = pars.regions.num;
    
    % Available features
    availableFeatures = fieldnames(pars.features);
    availableFeatures = availableFeatures(structfun(@(x)(x.enabled), pars.features));
    
    % Define random sampling distribution
    sz = [pars.dataset.imageHeight pars.dataset.imageWidth];
    mu = sz([2 1]) ./ [2 2.3];
    Sigma = [sz(2)*2 0; 0 sz(1)*3];
    ymap = NM_mvnpdf(sz(1), sz(2), mu, Sigma );
    
    % Fix random generator
    rng(2);
    
    % Sample patches
    centers = discretesample(ymap,numRegions);
    [i,j] = ind2sub(size(ymap), centers);
    randCoords = [i' j'];
    imzeros = zeros(pars.dataset.imageHeight, pars.dataset.imageWidth);
    randCoords = fix_coords(randCoords, imzeros, regionSize);
    boxPositions = get_patch_rect_from_coord(randCoords, regionSize);
    initialBoxPositions = boxPositions;
    initialBoxPositions = resize_boxes(dataset.images(:,:,:,1), initialBoxPositions, regionSize);
    
    % Init features
    [X, patches] =  extract_features_image(dataset.images(:,:,:,1), initialBoxPositions, availableFeatures, pars);
    X = zeros([dataset.count size(X,1) size(X,2)]);
    patches = zeros([dataset.count size(patches)], 'double');
    
    % Loop over all persons
    parfor ii=1:dataset.count
        
        % Get image & mask
        img = dataset.images(:,:,:,ii);
        mask = dataset.masks(:,:,ii);
        
        % Sample patches
        centers = discretesample(ymap,numRegions);
        [i,j] = ind2sub(size(ymap), centers);
        randCoords = [i' j'];
        imzeros = zeros(pars.dataset.imageHeight, pars.dataset.imageWidth);
        randCoords = fix_coords(randCoords, imzeros, regionSize);
        boxPositions = get_patch_rect_from_coord(randCoords, regionSize);
    
        % Extract features
        [X(ii, :, :), patches(ii,:,:,:,:)] =  extract_features_image(img, initialBoxPositions, availableFeatures, pars);
    end
    
    % Output structure
    features.X = X;
    features.patches = patches;
    
    % Save data
    try
        save(featuresFile, 'features');
    catch ME
        warning('extract_features:save', 'Unable to save features data on file %s.', featuresFile)
    end
end

% Features extraction time
fprintf('done in %.2f(s)\n', toc(t));

end




%% ========================================================================
%       PATCHES OPS
% =========================================================================
function [patches, patchesCoordsAsRowCol] = get_image_patches_from_boxes(boxes, img, regionSize, regionMask)
for p=1:size(boxes,1)
    patchesCoordsAsRowCol(p,:) = fix_coords(round([boxes(p,2)+(boxes(p,4)*0.5), boxes(p,1)+(boxes(p,3)*0.5)]), img, regionSize);
end
patches = get_image_patches(img, regionSize, regionMask, patchesCoordsAsRowCol);
end

function [randPatchesFeats, randPatchesCoords] = get_features_from_patches_around_previous_coords(img, regionSize, regionMask, patchesCoords, numNeighs, minNeighDist, maxNeighDist, availableFeatures, pars)
% Get random patches around previous tracked coords
[randPatches, randPatchesCoords] = get_image_patches_around_previous_coords(img, regionSize, regionMask, patchesCoords, numNeighs, minNeighDist, maxNeighDist);
% Extract features
for p=1:length(randPatches)
    randPatchesFeats{p} = get_features_from_patches(randPatches{p}, availableFeatures, pars);
end
end

function [randPatches, randPatchesCoords] = get_image_patches_around_previous_coords(img, regionSize, regionMask, patchesCoords, numNeighs, minNeighDist, maxNeighDist)
for p=1:size(patchesCoords,1)
    [randPatches{p}, randPatchesCoords{p}] = get_rand_coords_around_coord(img, regionSize, regionMask, patchesCoords(p,:), numNeighs, minNeighDist, maxNeighDist);
end
end

function [patches] = get_image_patches(image, regionSize, regionMask, patchesCoordsAsRowCol)
patches = cell(size(patchesCoordsAsRowCol,1), 1);
%imshow(image);
for p=1:size(patchesCoordsAsRowCol,1)
    rows = round(patchesCoordsAsRowCol(p,1)-(regionSize(1)*0.5));
    rows = rows:rows+regionSize(2)-1;
    cols = round(patchesCoordsAsRowCol(p,2)-(regionSize(2)*0.5));
    cols = cols:cols+regionSize(1)-1;
    try
        patches{p} = NM_immasked(image(rows,cols,:), regionMask, -inf);
    catch e
        a = 0;
    end
    %rectangle('Position', rect);
end
end

function [rect] = get_patch_rect_from_coord(patchesCoordsAsRowCol, regionSize)
rect = zeros(size(patchesCoordsAsRowCol,1),4);
for p=1:size(patchesCoordsAsRowCol,1)
    rect(p,1) = round(patchesCoordsAsRowCol(p,2)-(regionSize(2)*0.5));
    rect(p,2) = round(patchesCoordsAsRowCol(p,1)-(regionSize(1)*0.5));
    rect(p,[3 4]) = [regionSize(2) regionSize(1)];
end
end

function [] = plot_patches_rects(img, rects, color)
if nargin==2
    color = 'b';
end
if ~isempty(img)
    imshow(img);
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
end
for p=1:size(rects,1)
    rectangle('Position', rects(p,:), 'EdgeColor', color);
end
drawnow;
end


%% ========================================================================
%       FEATURE OPS
% =========================================================================
function type = get_feature_type(feature)
type = feature;
if any(strcmpi(feature, {'HSV', 'HSL', 'HSI', 'RGB', 'normRGB', 'LCH', 'Lab', 'Luv', 'YCbCr', 'YPbPr', 'XYZ', 'YUV'}))
    type = 'hist';
    %type = 'hist';
elseif strfind(feature, 'colorMean')
    type = 'colorMean';
elseif strfind(feature, 'sift')
    type = 'sift';
elseif strcmpi(feature, 'lmf')
    type = 'lm_filter';
elseif strcmpi(feature, 'mrf')
    type = 'mr_filter';
end
end

function [feat] = get_features_from_patches(patches,availableFeatures,pars)
feat = cell(1,length(patches));
for p=1:length(patches)
    processedPatch = NM_reid_image_preprocessing( patches{p}, ones(size(patches{p},1), size(patches{p},2)), pars.dataset.imageColorSpace, ...
        'HSV', 'eq', -inf);
    feat{p} = reshape(processedPatch(:,:,1), [], 1);
    
    %     for f=1:length(availableFeatures)
    %         type = get_feature_type(availableFeatures{f});
    %         processedPatch = NM_reid_image_preprocessing( patches{p}, ones(size(patches{p},1), size(patches{p},2)), pars.dataset.imageColorSpace, ...
    %             pars.features.(availableFeatures{f}).colorSpace, pars.features.(availableFeatures{f}).histOp, -inf);
    %         ff = NM_extractFeatures(processedPatch, type, pars.features.(availableFeatures{f}));
    %         if iscell(ff)
    %             feat{p} = [feat{p}; vertcat(ff{:})];
    %         else
    %             feat{p} = [feat{p}; ff];
    %         end
    %
    %     end
end
end


%% ========================================================================
%       COORDS OPS
% =========================================================================
function [coords] = get_rand_coords(image, regionSize, numPoints)
pad = regionSize + 2;
imCoords = zeros(size(image,1) - pad(1), size(image,2) - pad(2));
[r,c] = ind2sub(size(imCoords), randi(prod(size(imCoords)), numPoints,1));
coords = bsxfun(@plus, [r c], round(pad*0.5));
% Fix coordinates so teh patches may not lie out of the image
% boundaries
coords = fix_coords(coords, image, regionSize);
end

function [randPatches, randPatchesCoordsAsRowCol] = get_rand_coords_around_coord(img, regionSize, regionMask, coordAsRowCol, numNeighs, minNeighDist, maxNeighDist)
% Take random coords
randPatchesCoordsAsRowCol = repmat(coordAsRowCol, numNeighs,1) + randi([minNeighDist maxNeighDist], numNeighs, 2);
% Fix random coords so the patch is not out of image bounds
randPatchesCoordsAsRowCol = fix_coords(randPatchesCoordsAsRowCol, img, regionSize);
% Extract imag patches
randPatches = get_image_patches(img, regionSize, regionMask, randPatchesCoordsAsRowCol);
end

function [coordAsRowCol] = fix_coords(coordAsRowCol, img, regionSize)
for p=1:size(coordAsRowCol,1)
    % Fix rows
    if coordAsRowCol(p,1) - ceil(regionSize(1)*0.5) < 1
        coordAsRowCol(p,1) = ceil(regionSize(1)*0.5) + 1;
    end
    if coordAsRowCol(p,1) + ceil(regionSize(1)*0.5) > size(img,1)
        coordAsRowCol(p,1) = size(img,1)-ceil(regionSize(1)*0.5)-1;
    end
    
    % Fix cols
    if coordAsRowCol(p,2) - ceil(regionSize(2)*0.5) < 1
        coordAsRowCol(p,2) = ceil(regionSize(2)*0.5) + 1;
    end
    if coordAsRowCol(p,2)+ ceil(regionSize(2)*0.5) > size(img,2)
        coordAsRowCol(p,2) = size(img,2)-ceil(regionSize(2)*0.5)-1;
    end
end
end

function [fixedBoxes] = resize_boxes(img, boxes, regionSize)

coordsAsRowCol = round(boxes(:,[2 1])+(boxes(:,[4 3])*0.5));
fixedCoords = fix_coords(coordsAsRowCol, img, regionSize);
fixedBoxes = get_patch_rect_from_coord(fixedCoords, regionSize);

end



















