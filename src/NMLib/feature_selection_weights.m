function [W] = feature_selection_weights(A,B,th)
%#codegen
coder.inline('never')

%W = zeros(size(A,1), size(B,1));
%for ii=1:size(A,1)
%    for jj=1:size(B,1)
%        rho = mycorr(A(ii,:)',B(jj,:)');
%        W(ii,jj) = mean(rho(:));
%    end
%end

%tic
W = abs(corr(A',B'));
%toc

%tic
%mycorr(A,B);
%toc
W(isnan(W)) = 0;
W(W < th) = 0;

% Make it a markov matrix (so each column sums to 1)
%W = bsxfun(@rdivide, W, sum(W,1));

% Find the steady state (princip eigenvector)
%[AL,iteri] = principalEigenvectorRaw( Wbar , 0.001 );

end


function [rho] = mycorr(A,B)
% %#codegen
%coder.inline('never')
%Am = (A-mean(A))./std(A);
%Bm = (B-mean(B))./std(B);
%rho = (Am' * Bm) / sum(Am.^2);
for ii=1:size(A,1)
	for jj=1:size(B,1)
        rho(ii,jj) = fastcorrelation(A(ii,:)',B(jj,:)');
    end
end
end
