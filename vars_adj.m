function [vars_adj_rate] = vars_adj(u, A)
%% adjusted variance explained
% B = u'*A*u; 
% [U D V] = svd(B);
% Z = U*D.^0.5;  %% ZZ^T=B
% 
% [Q R] = qr(Z');
% r = diag(R);
% var_adj = norm(r)^2;
% vars_adj_rate = var_adj / trace(A);



% CEPV
% [Q, R]=qr(u,0);
% vars_adj_rate = trace(Q'*A*Q) / trace(A);

Am = u*(inv(u'*u))'*u'*A*u*(inv(u'*u))*u';
vars_adj_rate = trace(Am)/trace(A);
