function [Ymap] = NM_mvnpdf(rows, cols, mu, sigma)
% NM_MVNPDF Compute a multivariate gaussian PDF
% 
%   Copyright: Niki Martinel
%   Date: 01/13/2012
%   Output Data: multivariate gaussian PDF of size rows x cols
%   Parameters: 1. Number of rows
%               2. Number of columns
%               3. Mean (centre of PDF)
%               3. Covariance matrix (shape of PDF)
%
%   YMAP = NM_MVNPDF(ROWS, COLS,MU,SIGMA) returns the probability density 
%   of the multivariate normal distribution with mean MU and covariance 
%   SIGMA, evaluated at each row of X.
%   SIGMA is a D-by-D matrix, or an D-by-D-by-N array, in which case
%   the density is evaluated for each row of X with the corresponding page
%   of SIGMA, i.e., as MVNPDF matlab built-in function 
%   NM_MVGPD computes Y(I) using X(I,:) and SIGMA(:,:,I).
%   If the covariance matrix is diagonal, containing variances along the 
%   diagonal and zero covariances off the diagonal, SIGMA may also be
%   specified as a 1-by-D matrix or a 1-by-D-by-N array, containing 
%   just the diagonal. Pass in the empty matrix for MU to use its default 
%   value when you want to only specify SIGMA.
% 
%   If X is a 1-by-D vector, MVNPDF replicates it to match the leading
%   dimension of MU or the trailing dimension of SIGMA
%
%  Example:
%  
%       mu = [25 25]; Sigma = [.9 .4; .4 .3];
%       Ymap = NM_mvnpdf(50, 50, mu, Sigma );
%       [X1,X2] = meshgrid(1:50, 1:50);
%       surf(X1,X2,Ymap);
%      

% Meshgrid
[x1,x2] = meshgrid(1:cols,1:rows);
X = [x1(:) x2(:)];

% Compute multivariate gaussian PDF
Ymap = mvnpdf(X, mu, sigma);

% Reshape mvnpdf
Ymap = reshape(Ymap, rows, cols);

end