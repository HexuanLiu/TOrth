clear all; close all; clc

p = 100;

D = zeros(p,1);
D(1) = 1;
for j=2:8
    D(j)=D(j-1)-0.05;
end




for j=9:p
    D(j) = max(D(j-1)-0.01,0);
    %D(j)=0.5;
end

gamma = D(9)/D(8);

D = diag(D);
V = randn(p,p);

% orthonormalize
for j=1:p
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
E = 0.0000*randn(p,p);
E = (E+E')/2;
c = Sig+E;

%% block power method


A = Sig;
m = 8;
Q0 = randn(p,m);
[Q0, R0] = qr(Q0,0);
Vtrue = V(:, 1:m);

bound = [m-norm(Q0'*V(:,1:m),'fro')^2];
actual = [m-norm(Q0'*V(:,1:m),'fro')^2];
bound2= [m-norm(Q0'*V(:,1:m),'fro')^2];
V0 = V(:,1:m);
maxiter = 20;
d0 = norm(Q0*Q0'-V0*V0');
bound3 = [(m-norm(Q0'*V(:,1:m),'fro')^2)/(1-d0^2)];
theta = norm(Q0'*V(:,1:m),'fro')^2;

r = size(Q0,2);
Q=Q0;
for j=1:50
    Z=A*Q;    
%     for j=1:3
%         temp=truncate_operator(Z(:,j), d(j));        
%         Z(:,j)=temp;
%     end
    [Qnew, R]=qr(Z,0);
%    pr = m-norm(Q'*V(:,1:m),'fro')^2;
    
%     norm(Qnew'*V(:,1:3),'fro')^2
%     norm(Z'*V(:,1:3),'fro')^2/norm(R,2)^2
%    
     pr = bound(end);
     actual = [actual, m-norm(Qnew'*V(:,1:m),'fro')^2];
    
    bound = [bound, m*gamma^2*pr/((1-gamma^2)*(m-pr)+m*gamma^2)];
    
    pr2 = bound2(end);
    c = actual(end-1)*gamma^2/(1-(1-gamma^2)*norm(Q*Q'-Vtrue*Vtrue', 2)^2);
    
    
    cospr = norm(Q'*V(:,1:m),'fro')^2;
    bound2 = [bound2, c];
    %bound2 = [bound2, min(m, actual(end-1)*gamma^2/(1-norm(Qnew*Qnew'-Vtrue*Vtrue', 2)^2))];
    %bound2 = [bound2, m-((cospr)/((1-gamma^2)*(cospr)+m*gamma^2))];
    
    pr3 = bound3(end);
    bound3 = [bound3, pr3*gamma^2];
    Q = Qnew;
    
end

%% 
semilogy(actual,'b-+','linewidth',0.5)
hold on
semilogy(bound,'r-o','linewidth',0.5)
semilogy(bound2,'k-*','linewidth',0.5)
xlabel('Iteration t')
ylabel('||sin\Theta(P, Q_t)||_F^2')
legend('Observed', 'Approximated','Theoretical')
set(gca,'fontsize',20);



%%
p = 100;

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


% % case 1: complete overlap
% v1(1:10) = randn(10,1);
% v2(1:10) = randn(10,1);
% v3(1:10) = randn(10,1);

% case 2: partially overlap
v1(1:9)=ones(9,1);
v1(10)=-1;

v2(9) =-1;
v2(10:18)=-1;

v3(17)=-1;
v3(18:26)=1;

% %  % % case 3: non-overlap
% v1(1:10)=ones(10,1);
% v2(11:20)=ones(10,1);
% v3(21:30)=ones(10,1);

V(:,1)=v1/norm(v1);
V(:,2)=v2/norm(v2);
V(:,3)=v3/norm(v3);
 
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

E = 0.075*randn(p,p); %0.014, 0.035, 0.075
E = (E+E')/2;
c = Sig+E;

%% truncation after orthogonalization

Q0 = randn(p,3);
[Q0, R]=qr(Q0,0);
Q = Q0;

Qlist1 = [];

for j=1:20
    Z=c*Q;    
    for j=1:3
        temp=truncate_operator(Z(:,j), p);        
        Z(:,j)=temp;
    end
    
    [Qnew, R]=qr(Z,0);
    
    for j=1:3
        temp=truncate_operator(Qnew(:,j), p);        
        Q(:,j)=temp;
    end
    
    angle = norm(Q-Qnew,'fro')^2;

    Qlist1 = [Qlist1, angle];

end

for j=1:20
    Z=c*Q;    
    for j=1:3
        temp=truncate_operator(Z(:,j), 50);        
        Z(:,j)=temp;
    end
    
    [Qnew, R]=qr(Z,0);
    
    for j=1:3
        temp=truncate_operator(Qnew(:,j), 50);        
        Q(:,j)=temp;
    end
    
    angle = norm(Q-Qnew,'fro')^2;

    Qlist1 = [Qlist1, angle];

end

for j=1:20
    Z=c*Q;    
    for j=1:3
        temp=truncate_operator(Z(:,j), 25);        
        Z(:,j)=temp;
    end
    
    [Qnew, R]=qr(Z,0);
    
    for j=1:3
        temp=truncate_operator(Qnew(:,j), 25);        
        Q(:,j)=temp;
    end
    
    angle = norm(Q-Qnew,'fro')^2;

    Qlist1 = [Qlist1, angle];

end

for j=1:20
    Z=c*Q;    
    for j=1:3
        temp=truncate_operator(Z(:,j), 10);        
        Z(:,j)=temp;
    end
    
    [Qnew, R]=qr(Z,0);
    nnz(Qnew)
    
    for j=1:3
        temp=truncate_operator(Qnew(:,j), 10);        
        Q(:,j)=temp;
    end
    
    angle = norm(Q-Qnew,'fro')^2;

    Qlist1 = [Qlist1, angle];

end
plot([1:80],Qlist1,'r-o','linewidth',0.5)

xlabel('Iteration t')
xlim([1 80])
set(gca,'XTick',[1 20 40 60 80])
ylabel('||Q_{truncate}-Q||_F^2')

set(gca,'fontsize',20);
%% whether truncation at each step help with convergence and accuracy

Q0 = randn(p,3);
[Q0, R]=qr(Q0,0);
Q = Q0;

Qlist1 = [];

for j=1:20
    Z=c*Q;    
    for j=1:3
        temp=truncate_operator(Z(:,j), p);        
        Z(:,j)=temp;
    end
    
    [Qnew, R]=qr(Z,0);
    angle = 3-norm(Qnew'*Q, 'fro')^2;
    angle_true = 3-norm(Qnew'*V(:,1:3), 'fro')^2;
    Qlist1 = [Qlist1, angle_true];

    Q = Qnew;
end

Q = Q0;
Qlist2 = [];
for j=1:20
    Z=c*Q;    
    for j=1:3
        temp=truncate_operator(Z(:,j), 10);        
        Z(:,j)=temp;
    end
    
    [Qnew, R]=qr(Z,0);
    angle = 3-norm(Qnew'*Q, 'fro')^2;
    angle_true = 3-norm(Qnew'*V(:,1:3), 'fro')^2;
    Qlist2 = [Qlist2, angle_true];

    Q = Qnew;
end

plot(Qlist1(1:end), 'r')
hold on
plot(Qlist2(1:end), 'g')

%% different truncation level
Q = randn(p,3);
[Q,R]=qr(Q,0);
Qlist = [];
Qtrue = [];
for i=1:20
    Z=c*Q;    
    for j=1:3
        temp=truncate_operator(Z(:,j), p);        
        Z(:,j)=temp;
    end
    
    [Qnew, R]=qr(Z,0);
    angle = 3-norm(Qnew'*Q, 'fro')^2;
    angle_true = 3-norm(Qnew'*V(:,1:3), 'fro')^2;
    Qlist = [Qlist, angle];
    Qtrue = [Qtrue, angle_true];
    
%     if(norm(Qnew-Q)<1e-6)
%         i
%         break
%     end
    Q = Qnew;
end

for i=1:20
    Z=c*Q;    
    for j=1:3
        temp=truncate_operator(Z(:,j), p/2);        
        Z(:,j)=temp;
    end
    
    [Qnew, R]=qr(Z,0);
    angle = 3-norm(Qnew'*Q, 'fro')^2;
    angle_true = 3-norm(Qnew'*V(:,1:3), 'fro')^2;
    Qlist = [Qlist, angle];
    Qtrue = [Qtrue, angle_true];
    
    Q = Qnew;
end

for i=1:20
    Z=c*Q;    
    for j=1:3
        temp=truncate_operator(Z(:,j), 25);        
        Z(:,j)=temp;
    end
    
    [Qnew, R]=qr(Z,0);
    angle = 3-norm(Qnew'*Q, 'fro')^2;
    angle_true = 3-norm(Qnew'*V(:,1:3), 'fro')^2;
    Qlist = [Qlist, angle];
    Qtrue = [Qtrue, angle_true];
    
    Q = Qnew;
end

for i=1:20
    Z=c*Q;    
    for j=1:3
        temp=truncate_operator(Z(:,j), 10);        
        Z(:,j)=temp;
    end
    
    [Qnew, R]=qr(Z,0);
    angle = 3-norm(Qnew'*Q, 'fro')^2;
    angle_true = 3-norm(Qnew'*V(:,1:3), 'fro')^2;
    Qlist = [Qlist, angle];
    Qtrue = [Qtrue, angle_true];
    
    Q = Qnew;
end

for i=1:20
    Z=c*Q;    
    for j=1:3
        temp=truncate_operator(Z(:,j), 5);        
        Z(:,j)=temp;
    end
    
    [Qnew, R]=qr(Z,0);
    angle = 3-norm(Qnew'*Q, 'fro')^2;
    angle_true = 3-norm(Qnew'*V(:,1:3), 'fro')^2;
    Qlist = [Qlist, angle];
    Qtrue = [Qtrue, angle_true];
    
    Q = Qnew;
end

rhoE = norm(E)*ones(length(Qlist)-4,1);
xlist = [5:length(Qlist)];

plot(xlist,Qlist(5:end),'b-+','linewidth',0.5)
hold on
plot(xlist,Qtrue(5:end),'r-o','linewidth',0.5)
plot(xlist, rhoE, 'k-', 'linewidth',0.5)
xlabel('Iteration t')
xlim([5 100])
set(gca,'XTick',[5 20 40 60 80 100])
ylabel('||sin\Theta(P, Q_t)||_F^2')
legend('P=Q_{t-1}', 'P=V','\rho(E)')
set(gca,'fontsize',20);