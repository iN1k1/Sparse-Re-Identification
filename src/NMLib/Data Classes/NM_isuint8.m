function [isuint8] = NM_isuint8(data)
% NM_ISUINT8 Check if input data is of UINT8 class type
% 
%   Copyright: Niki Martinel
%   Date: 09/08/2011
%   Return Data: True if input data is of UINT8 class, false otherwise
%   Parameters: 1. any data
%               
%   [ISDOUBLE] = NM_ISUINT8(DATA) takes input data, DATA, and
%   check if it is of UINT8 class
%
%   DATA should be any kind of data
%   
%   OUT is a TRUE/FALSE if input data is of UINT8 class type or not
%

isuint8 = strcmpi(class(data), 'uint8');

end