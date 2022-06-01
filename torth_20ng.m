clear all; clc
load('20news_w100.mat')
X = double(documents);
y = newsgroups;
n = size(X,2);
Xm = mean(X');
X = X'-repmat(Xm, n, 1);
A = X'*X;

%%
%A = X*X';
Q0 = randn(100,4);
[Q0, R0] = qr(Q0, 0);

klist = [10, 10, 10, 10];
% for j=1:200
%     Z=A*(A'*Q);    
%     for i=1:size(Q0,2)
%         Z(:,i)=truncate_operator(Z(:,i), klist(i));        
%     end
% 
%     [Qnew, R]=qr(Z,0);
%     for i=1:size(Q0,2)
%         Qnew(:,i)=truncate_operator(Qnew(:,i), klist(i));        
%     end
%     Q = Qnew;
% end
% %     for i=1:size(Q0,2)
% %         Q(:,i)=truncate_operator(Q(:,i), klist(i));        
% %     end


%Q3 = ortht(Q0,  [100,100,100], A);

Q4 = ortht(Q0, [100,100, 100, 100], klist, A);
%     for i=1:size(Q4,2)
%         Qnew(:,i)=truncate_operator(Q4(:,i), klist(i));        
%     end
sparse(Q4)
projected = X*Q4;
xx = projected(:,1);
yy = projected(:,2);
zz = projected(:,3);
dd = projected(:,4);
ind1 = find(y==1);
ind2 = find(y==2);
ind3 = find(y==3);
ind4 = find(y==4);

figure
% scatter3(xx(ind1), yy(ind1), zz(ind1), 'co')
% hold on
% scatter3(xx(ind2), yy(ind2), zz(ind2), 'r<')
% scatter3(xx(ind3), yy(ind3), zz(ind3), 'm*')
% scatter3(xx(ind4), yy(ind4), zz(ind4), 'gs')

scatter(xx(ind1), yy(ind1), 'bo')
% hold on
% scatter(xx(ind2), yy(ind2), 'r<')
% scatter(xx(ind3), yy(ind3), 'm*')
% scatter(xx(ind4), yy(ind4), 'm<')
%legend('comp', 'rec', 'sci', 'talk')
legend('comp')
xlabel('PC1')
ylabel('PC2')
set(gca,'fontsize',20);
xlim([-0.5 3])
ylim([-3 0.5])
% xlim([-5 1.5])
% ylim([-3 4])
grid on
% set(gca,'XTick',[1 1 100])
% set(gca,'YTick',[1:1:100])
%%
clc
ind = find(Q4(:,2));
wordlist{ind}

%variance_rate = vars_adj(Q4,A)