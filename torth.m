function [Q] = torth(Q0, klist, A, maxiter)

% Q0: initialization matrix of size p x m
% klist: the list of integers for cardinality constriants. length(klist)=m
% A: positive semidefinite matrix. For sparse PCA, A is the covariance
% matrix
% maxiter: default is 200

if nargin < 4

  maxiter = 200;

end

Q = Q0;

for j=1:maxiter
    Z=A*Q;    
    for i=1:size(Q0,2)
        Z(:,i)=truncate_operator(Z(:,i), klist(i));        
    end

    [Qnew, R]=qr(Z,0);

    
    if (norm(Qnew-Q)<1e-12)
        break
    end
    Q = Qnew;
end

end

