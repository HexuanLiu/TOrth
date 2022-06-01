clear all; close all; clc

p = 1000;

D = zeros(p,1);
D(1) = 1;
D(2) = 0.9;
D(3) = 0.8;


for j=4:p
    D(j)=0.1;
end

D = diag(D);
V = randn(p,p);

v1 = zeros(p,1);
v2 = zeros(p,1);
v3 = zeros(p,1);


%% case 1: complete overlap

% v1(1:10) = randn(10,1);
% v2(1:10) = randn(10,1);
% v3(1:10) = randn(10,1);

%case 2: partially overlap
v1(1:9)=ones(9,1);
v1(10)=-1;

v2(9) =-1;
v2(10:18)=-1;

v3(17)=-1;
v3(18:26)=1;

%     v1(1:9)=ones(9,1);
%     v1(10)=-1;
% 
%     v2(9) =-1;
%     v2(10)=1;
%     v2(10:18)=1;
% 
%     v3(10)=-1;
%     v3(11:19)=1;


    


%  % case 3: non-overlap
% v1(1:10)=ones(10,1);
% v2(11:20)=ones(10,1);
% v3(21:30)=ones(10,1);
% 
V(:,1)=v1/norm(v1);
V(:,2)=v2/norm(v2);
V(:,3)=v3/norm(v3);
%  
% orthonormalize
for j=2:p
    vj=V(:,j);
    for i=1:j-1
        rij = V(:,i)'*V(:,j);
        vj=vj-rij*V(:,i);
    end
    rij = norm(vj);
    qj=vj/rij;
    V(:,j)=qj;
end

Sig = V*D*V';
Sig = (Sig+Sig')/2;

% X = mvnrnd(zeros(p,1),Sig,10);
% Xm = mean(X);
% X = X-repmat(Xm, 10, 1);

E = 0.015*randn(p,p);
E = E'*E;
c = Sig+E;

%% standard+simple truncation
maxiter = 10;
success = 0;
recovery = 0;
bl = zeros(3,maxiter);
m = 3;
k = 10;

tic
for rep=1:maxiter
    A = c;
    Q0 = randn(p,3);
    [Q0, R] = qr(Q0, 0);
%     [evec, eval] = eig(A);
%     [~, ind] = sort(diag(eval));
%     Q1 = evec(:, ind(end-2:end));
    Q1 = torth(Q0, [p,p,p], A);
%     Q2 = torth(Q1, [4*k, 4*k, 4*k], A);
%     Q3 = torth(Q2, [2*k, 2*k, 2*k], A);
    Q4 = Q1;
    for i=[1:3]
        Q4(:,i)=truncate_operator(Q1(:,i),10);
    end
    

%      Q4 = torth(Q0, [p, p, p], A);
    
    bl(1, rep)= abs(Q4(:,1)'*V(:,1));    
    bl(2, rep)= abs(Q4(:,2)'*V(:,2));
    bl(3, rep)= abs(Q4(:,3)'*V(:,3));
    
    if (isequal(find(Q4(:,1)), find(V(:,1)))&& isequal(find(Q4(:,2)), find(V(:,2)))&& isequal(find(Q4(:,3)), find(V(:,3))))
        recovery = recovery+1;
    end
    if min(bl(:,rep)>0.99)
        success = success+1;
    end
end
mean(bl(1,:))
mean(bl(2,:))
mean(bl(3,:))

success
recovery

toc

%% Tpower
maxiter = 1000;
recovery = 0;
success = 0;
al = zeros(3,maxiter);
m = 3;
k = 10;

tic

for rep=1:maxiter
    
    A = c;
    flag = 1;
    for i=1:m
        x0 = randn(p,1);
        x0 = x0/norm(x0);
        
        x1 = tp(x0, p, A);
%         x2 = tp(x1, 4*k, A);
%         x3 = tp(x2, 2*k, A);
        x4 = tp(x1, k, A);

        acc = abs(x4'*V(:,i));
        al(i,rep)=acc;
        P = eye(p)-x4*x4';

        A = P*A*P;
       
        if (~isequal(find(x4), find(V(:,i))))
            flag=0;
        end
    end
    if min(al(:,rep)>0.99)
        success = success+1;
    end
    if(flag==1)
       recovery = recovery+1;
    end
end

% for rep = 1:maxiter
%     for ck = [8*k, 4*k 2*k, k]
%         A = c;
%         for i = 1:m
%             x0 = randn(p,1);
%             x0 = x0/norm(x0);
%             x1 = tp(x0, ck, A);
%             P = eye(p)-x4*x4';
%             A = P*A*P;
%         end
%     end
% end
mean(al(1,:))
mean(al(2,:))
mean(al(3,:))

success
recovery
toc
%% Torth 1
maxiter = 10;
success = 0;
recovery = 0;
bl = zeros(3,maxiter);
m = 3;
k = 10;

tic
for rep=1:maxiter
    A = c;
    Q0 = randn(p,3);
    [Q0, R] = qr(Q0, 0);
    Q1 = TOrthT(Q0, [p,p,p], A);
    Q2 = TOrthT(Q1, [4*k, 4*k, 4*k], A);
    Q3 = TOrthT(Q2, [2*k, 2*k, 2*k], A);
    Q4 = TOrthT(Q3, [k, k, k], A);

%      Q4 = torth(Q0, [p, p, p], A);
    
    bl(1, rep)= abs(Q4(:,1)'*V(:,1));    
    bl(2, rep)= abs(Q4(:,2)'*V(:,2));
    bl(3, rep)= abs(Q4(:,3)'*V(:,3));
    
    if (isequal(find(Q4(:,1)), find(V(:,1)))&& isequal(find(Q4(:,2)), find(V(:,2)))&& isequal(find(Q4(:,3)), find(V(:,3))))
        recovery = recovery+1;
    end
    if min(bl(:,rep)>0.99)
        success = success+1;
    end
end
mean(bl(1,:))
mean(bl(2,:))
mean(bl(3,:))

success
recovery

toc
%% torth 2
success = 0;
cl = zeros(3,100);
m = 3;
k = 10;

tic
for rep=1:10
    A = c;
    Q0 = randn(p,3);
%    Q1 = ortht(Q0, [p, p, p], A);
%     Q1 = ortht(Q0, [8*k, 8*k, 8*k], A);
%     Q2 = ortht(Q1, [4*k, 4*k, 4*k], A);
%     Q3 = ortht(Q2, [2*k, 2*k, 2*k], A);
    Q4 = ortht(Q0, [8*k, 8*k, 8*k],[k, k, k], A);
    
    cl(1, rep)= abs(Q4(:,1)'*V(:,1));    
    cl(2, rep)= abs(Q4(:,2)'*V(:,2));
    cl(3, rep)= abs(Q4(:,3)'*V(:,3));
    
%     if(isequal(find(Q4(:,1)), find(V(:,1)))&& isequal(find(Q4(:,2)), find(V(:,2)))&& isequal(find(Q4(:,3)), find(V(:,3))))
%         success = success+1;
%    end
    if min(cl(:,rep)>0.9)
        success = success+1;
    end
end
toc


%% Pitprops
load('pitprops.mat')
klist = [7,2,4,3,5,4];
%klist = [6,2,1,2,1,1];
%klist = [7,2,1,1,1,1];
A = covX;
V = zeros(13,6);

%% zou and hastie spca
tic
[Q SD] = spca_zouhastie([], A, 6, 1e-4, -klist, 100, 1e-6, true);
toc

%%
%klist = [7,2,4,3,5,4];
klist = [6,2,1,2,1,1];
p = 13;
exp_variance = 0;
restp=[];
results=[];
for rep=1:1
    A = covX;
    Q0 = zeros(13, 6);
    for i=1:6
%             [val,idx]=max(diag(A));
%             x0 = zeros(size(A,1),1);
%             x0(idx) = 1;
%             Q0(:,i)=x0;
    x0 = randn(p,1);
    x0 = x0/norm(x0);
    Q0(:,i)=x0;
    %xp = tp(x0, p, A);
%     x1 = tp(x0, min(p,8*klist(i)), A);
%     x2 = tp(x1, min(p,4*klist(i)), A);
%     x3 = tp(x2, min(p,2*klist(i)), A);
%     x4 = tp(x3, klist(i), A);
%     
    
    x1 = tp(x0, p, A);
    x2 = tp(x1, max(round(p/2), klist(i)), A);
    x3 = tp(x2, max(round(p/4), klist(i)), A);
    x4 = tp(x3, max(round(p/8), klist(i)), A);
    x5 = tp(x4, klist(i), A);
    
    
    
    V(:,i)=x5;    
    P = eye(p)-x5*x5';
    A = P*A*P;

end
exp_variance=vars_adj(V, covX)
restp=[restp, exp_variance];


% A = covX;
% % Q0 = rand(p,6);
% % [Q0, R]=qr(Q0, 0);
% Q4 = ortht(Q0, klist,klist, A);
% expv=vars_adj(Q4,covX)
% results = [results, expv];
end
% % figure
% % plot(restp, 'bo')
% % hold on
% % plot(results, 'r*')
% figure
% histogram(restp, 'BinWidth',0.01, 'FaceColor','blue', 'FaceAlpha',0.3)
% hold on
% histogram(results, 'BinWidth',0.01, 'FaceColor','red', 'FaceAlpha',0.3)
% xlabel('Proportion of explained variance')
% ylabel('Number of trials')
% legend('TPower', 'TOrth')
% set(gca,'fontsize',20);
%%
p=13;
%klist = [7,2,4,3,5,4];
klist = [6,2,1,2,1,1];
%klist = [1,7,2,1,1,1];

A = covX;
Q0 = randn(p,6);
[Q0, R0]=qr(Q0,0);

% Q1 = torth(Q0, [13, 13, 13, 13, 13, 13], A);
% Q2 = torth(Q1, [7, 7, 7, 7, 7, 7], A);
% Q3 = torth(Q2, [6, 4, 4, 4, 4, 4], A);
% Q4 = torth(Q3, [6, 2, 2, 2, 2, 2], A);
% Q5 = torth(Q4, klist, A);
% 
% vars_adj(Q5,covX)


% Q0 = zeros(p,6);
% for i=[1:1]
%     Q0(i,i)=1/100;
% end

% % init in GPower
% S = chol(covX);
% n = p;
% m = 6;
% norm_a_i=zeros(p,1);
%     for i=1:n,
%         norm_a_i(i)=norm(S(:,i));
%     end
%     [ignore,i_max]=max(norm_a_i);
%     [Q0,rho_max]=qr([S(:,i_max)/norm_a_i(i_max), randn(p,m-1)],0); %initialization point
    
    
%    results=[];
% for rep=1
% 
% %Q0 = rand(p,6);
% 
%Q0 = torth(Q0,[13,13,13,13,13,13], A);
% Q1 = torth(Q0,[13,13,8,13,8,8], A);
% Q2 = torth(Q1,[13,8,4,8,4,4], A);
% Q3 = torth(Q2, [13,4,2,4,2,2], A);
% Q4 = TOrthT(Q0, klist, A);
Q4 = ortht(Q0, [26,26,26,26,26,26],klist, A);

expv=vars_adj(Q4,covX)
% results = [results, expv];
% end



%% result reported in rSVD by Shen and Huang
PC = zeros(13, 6);
PC(:,1)=[-0.449,-0.460,0,0,0,-0.199, -0.399, -0.279, -0.380, -0.407, 0,0,0];
PC(3:4,2)=[-0.707, -0.707];
PC(5:7,3)=[0.550, 0.546, 0.366];
PC(13,3)=-0.515;
PC(:,4)=[-0.114, -0.102, 0,0,0, -0.176, 0, 0.422,0,0.283, 0, -0.785, -0.265];
PC(10:11,5)=[0.231, -0.973];
PC(5, 6)=-0.744;
PC(12:13,6)=[0.161, -0.648];
vars_adj(PC, covX)

%% ITSPCA by MA
[PC, k, nstep, d, sigma2_hat, Q0] = itspca(A, 5.5, 1.5,'hard_thresholding',6, 1e-4);
vars_adj(PC, covX)

%% SPCA
PC = zeros(13, 6);
PC(:,1)=[-0.477, -0.476, 0, 0, 0.177,0, -0.250, -0.344, -0.416, -0.400, 0,0,0];
PC(:,2)=[0, 0, 0.785, 0.620, 0, 0,0, -0.021, 0, 0, 0, 0.013, 0];
PC(5:7,3)=[0.640, 0.589, 0.492];
PC(13,3)=-0.015;
PC(11,4)=-1;
PC(12,5)=-1;
PC(13,6)=1;
vars_adj(PC(:,1:6), covX)

%% TPower-13

PC = zeros(13,6);
PC(:,1)=[0.4444, 0.4534, 0, 0, 0, 0, 0.3779, 0.3415, 0.4032, 0.4183, 0, 0, 0];
PC(3:4, 2)=[0.7071, 0.7071];
PC(5,3)=1;
PC(6:7,4)=[0.8569, 0.5154];
PC(11,5)=1;
PC(12,6)=1;
vars_adj(PC(:,1:6), covX)