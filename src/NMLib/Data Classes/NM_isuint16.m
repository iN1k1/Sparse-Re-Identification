function [isuint16] = NM_isuint16(data)
% NM_ISUINT16 Check if input data is of UINT16 class type
% 
%   Copyright: Niki Martinel
%   Date: 09/08/2011
%   Return Data: True if input data is of UINT16 class, false otherwise
%   Parameters: 1. any data
%               
%   [ISDOUBLE] = NM_ISUINT16(DATA) takes input data, DATA, and
%   check if it is of UINT16 class
%
%   DATA should be any kind of data
%   
%   OUT is a TRUE/FALSE if input data is of UINT16 class type or not
%

isuint16 = strcmpi(class(data), 'uint16');

end