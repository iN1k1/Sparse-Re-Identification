function [image] = NM_reid_image_preprocessing( image,  inputColorSpace, ...
                                                outputColorSpace )
% Author:    Niki Martinel
% Date:      2012/10/10 08:30:47
% Revision:  0.1
% Copyright: Niki Martinel, 2012


if NM_isdouble(image) || NM_issingle(image) || isempty(find(image(:)>1,1))
    ;
else
    image = im2double(image);
end

% Convert image
image = NM_colorconverter(image, inputColorSpace, outputColorSpace);

end
