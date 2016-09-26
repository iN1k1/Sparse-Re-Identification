function [isdouble] = NM_isdouble(data)
% NM_ISDOUBLE Check if input data is of double class type
% 
%   Copyright: Niki Martinel
%   Date: 09/08/2011
%   Return Data: True if input data is of double class, false otherwise
%   Parameters: 1. any data
%               
%   [ISDOUBLE] = NM_ISDOUBLE(DATA) takes input data, DATA, and
%   check if it is of double class
%
%   DATA should be any kind of data
%   
%   OUT is a TRUE/FALSE if input data is of double class type or not
%

isdouble = strcmpi(class(data), 'double');

end