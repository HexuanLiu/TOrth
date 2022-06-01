clear all; close all; clc

sst = ncread('sst.wkmean.1990-present.nc','sst');
%%
p = 360*180;
n = 1455;

X=zeros(p,n);
for j=1:n
   X(:,j)=reshape(sst(:,:,j),[p,1]);   
end

%%
A = mean(X');
X = X'-repmat(A, n, 1);

%%
tic
[U,S,V]=svd(X,'econ');
toc
%%
V4 = reshape(V(:,4), [360, 180]);
masked = V4.*mask;
imagesc(masked')
%%
% A = mean(X');
% X = X'-repmat(A, n, 1);

% X = X';
%% TPower
r = [p,p,p,800];
maxiter = 100;
Q = zeros(p,4); 
tic
Xn= X;
for i=1:4
    Xop = @(y)Xn*y;
    Xtop = @(y)Xn'*y;
    
    x0 = randn(p,1);
    x0 = x0/norm(x0);        
    x = x0;
    for j = 1:maxiter
        z = Xop(x);
        xnew = Xtop(z);
        temp = truncate_operator(xnew, r(i));
        xnew = temp/norm(temp); 
        if (norm(xnew-x)<1e-4)
            
            break
        end
        x = xnew;
        Q(:,i)=x;
    end
        
    Xn = Xn - (Xn*x)*x';     
end
toc

%% truncated orthogonal iteration
r = [p, p, p, 800]; % number of nonzeros in each component
Q=randn(p,4);
maxiter = 100;
Xop = @(y)X*y;
Xtop = @(y)X'*y;
tic
for j=1:maxiter
    y = Xop(Q);
    Z = Xtop(y);   
    for k=1:length(r)
        temp=truncate_operator(Z(:,k), r(k));        
        Z(:,k)=temp;
    end
    [Qnew, R]=qr(Z,0);
    if (norm(Qnew-Q)<1e-4)
        break;
    end
    
    Q = Qnew;
    %Q = sparse(Q);
end
toc
% for k=1:length(r)
%     temp=truncate_operator(Q(:,k), r(k));        
%     Q(:,k)=temp/norm(temp);
% end

%% relax and split
tic
[B, A, obj]=RS_PCA(X,4, 1e-4,1e-5, 1e-6,200);

toc

%% GPower block
tic
m=4;
gamma=0.1*ones(1,m);  
[Q] = GPower(X,[0,0,0,0.11],m,'l1',0,[1,1,1,1]);
toc


%% zou hastie spca
tic
[Q SD] = spca_zouhastie(X, [], 4, 1e-4, -[p, p, p, 800], 2, 1e-3, false);
toc
%%
mask = ncread('lsmask.nc','mask');
q1 = Q(:,1);
q1 = reshape(q1, [360,180]);
masked = q1.*mask;
imagesc(masked')

%%
q2 = Q(:,2);
q2 = reshape(q2, [360,180]);
masked = q2.*mask;
imagesc(masked')

%%
q3 = Q(:,3);
q3 = reshape(q3, [360,180]);
masked = q3.*mask;
imagesc(masked')

%%
Q = SL;
q4 = Q(:,4);
coeff = X*q4;
figure(1)
plot(coeff)
hold on
plot(zeros(n,1))
XTick = [1:45:1455];
set(gca,'xtick',XTick)
%%

q4 = Q(:,4);
q4 = reshape(q4, [360,180]);
masked = q4.*mask;
figure(2)
imagesc(masked')

%%

q4 = truncate_operator(V(:,4), 1000);

q4 = reshape(q4, [360,180]);
masked = q4.*mask;
figure(1)
imagesc(masked')
