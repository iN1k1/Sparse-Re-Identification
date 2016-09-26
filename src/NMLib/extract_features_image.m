function [features, patches] = extract_features_image(image, regions, availableFeatures, pars)

% Convert to double!
if ~NM_isdouble(image)
    image = im2double(image);
end

%-------------------------
% Init tracked patch accumulated features array
[ff, pp] = get_features_from_regions(image, regions(1,:), availableFeatures, pars);
features = NaN*ones(size(regions,1), length(ff));
patches = zeros([size(pp) size(regions,1)], 'double');
clear ff pp;

% For all patches..
for p=1:size(regions,1)
    [features(p,:), patches(:,:,:,p)] = get_features_from_regions(image, regions(p,:), availableFeatures, pars);
end

end


function [feats, patches] = get_features_from_regions(img, regions, availableFeatures, pars)
patches = get_image_patches_from_regions(regions, img);
% Extract features
for p=1:size(regions,1)
    feats{p} = get_features_from_patch(patches(:,:,:,p), availableFeatures, pars);
end
feats = [feats{:}];
end   

function [patches, patchesCoordsAsRowCol] = get_image_patches_from_regions(regions, img)
for p=1:size(regions,1)
    patchesCoordsAsRowCol(p,:) = fix_coords(round([regions(p,2)+(regions(p,4)*0.5), regions(p,1)+(regions(p,3)*0.5)]), img, regions(p,[4 3]));
end
patches = get_image_patches(img, regions(:,[4 3]), patchesCoordsAsRowCol);
end

function [patches] = get_image_patches(image, regionSize, patchesCoordsAsRowCol)
patches = zeros([regionSize size(image,3), size(patchesCoordsAsRowCol,1)], 'like', image);
%imshow(image);
for p=1:size(patchesCoordsAsRowCol,1)
    rows = round(patchesCoordsAsRowCol(p,1)-(regionSize(p,1)*0.5));
    rows = rows:rows+regionSize(p,1)-1;
    cols = round(patchesCoordsAsRowCol(p,2)-(regionSize(p,2)*0.5));
    cols = cols:cols+regionSize(p,2)-1;  
    try
        patches(:,:,:,p) = image(rows,cols,:);
    catch e
        warning(e);
    end
    %rectangle('Position', rect);
end
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

function [feat] = get_features_from_patch(patch,availableFeatures,pars)
feat = [];
for f=1:length(availableFeatures)
    type = get_feature_type(availableFeatures{f});
    
    % Color conversion
    processedPatch = NM_reid_image_preprocessing( patch, pars.dataset.imageColorSpace, pars.features.(availableFeatures{f}).colorSpace);
    
    % Feature extraction
    ff = NM_extractFeatures(processedPatch, type, pars.features.(availableFeatures{f}));
    if iscell(ff)
        feat = [feat; vertcat(ff{:})];
    else
        feat = [feat; ff];
    end
end
end

function type = get_feature_type(feature)
type = feature;
if any(strcmpi(feature, {'HSV', 'HSL', 'HSI', 'RGB', 'normRGB', 'LCH', 'Lab', 'Luv', 'YCbCr', 'YPbPr', 'XYZ', 'YUV', 'grayscale'}))
    type = 'hist';
end
end
