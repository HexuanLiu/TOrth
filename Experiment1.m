clear all; close all; clc

p = 1000;
m = 3;
k = 10;

D = zeros(p,1);
D(1) = 1;
D(2) = 0.9;
D(3) = 0.8;


for j=4:p
    D(j)=0.1;
end

D = diag(D);

maxiter = 100;
success0 = 0;
recovery0 = 0;
bl0 = zeros(3,maxiter);
   
success1 = 0;
recovery1 = 0;
bl1 = zeros(3,maxiter);

success2 = 0;
recovery2 = 0;
bl2 = zeros(3,maxiter);

success3 = 0;
recovery3 = 0;
bl3 = zeros(3,maxiter);

success4 = 0;
recovery4 = 0;
bl4 = zeros(3,maxiter);

truepos = zeros(5,1);
falsepos = zeros(5,1);
falseneg = zeros(5,1);

for rep=1:maxiter
    V = randn(p,p);

    v1 = zeros(p,1);
    v2 = zeros(p,1);
    v3 = zeros(p,1);
    
    % % case 1: complete overlap
    v1(1:10) = randn(10,1);
    v2(1:10) = randn(10,1);
    v3(1:10) = randn(10,1);
%     
%     % %case 2: partially overlap
%     v1(1:9)=ones(9,1);
%     v1(10)=-1;
% 
%     v2(9) =-1;
%     v2(10:18)=-1;
% 
%     v3(17)=-1;
%     v3(18:26)=1;

%     % case 3: non-overlap
%     v1(1:10)=ones(10,1);
%     v2(11:20)=ones(10,1);
%     v3(21:30)=ones(10,1);
%     
    
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
    E = 0.015*randn(p,p);
    E = E'*E;
    c = Sig+E;
    
    % standard

    A = c;
    Q0 = randn(p,3);
    [Q0, R] = qr(Q0, 0);
    Q4 = torth(Q0, [p,p,p], A);
    
    
    bl0(1, rep)= abs(Q4(:,1)'*V(:,1));    
    bl0(2, rep)= abs(Q4(:,2)'*V(:,2));
    bl0(3, rep)= abs(Q4(:,3)'*V(:,3));
    
    list1 = find(V(:,1:3));
    list2 = find(Q4(:,1:3));
    list3 = find(~V(:,1:3));
    
    
    for kindex =1:length(list1)       
        if ismember(list1(kindex), list2)==1 % hit
            truepos(1)=truepos(1)+1;
        else
            falseneg(1) = falseneg(1)+1;
        end
    end
    
    % if a true zero is identified as nonzero, fp
    for kindex = 1:length(list3)
        
        if ismember(list3(kindex), list2)==1
            falsepos(1)= falsepos(1)+1;
        end
    end
    
        
%     if (isequal(find(Q4(:,1)), find(V(:,1)))&& isequal(find(Q4(:,2)), find(V(:,2)))&& isequal(find(Q4(:,3)), find(V(:,3))))
%         recovery0 = recovery0+1;
%    end
    if min(bl0(:,rep)>0.99)
        success0 = success0+1;
    end
    
    % standard+simple truncation
    A = c;
    Q1 = torth(Q0, [p,p,p], A);
    Q4 = Q1;
    for i=[1:3]
        Q4(:,i)=truncate_operator(Q1(:,i),k);
    end
    
    
    bl1(1, rep)= abs(Q4(:,1)'*V(:,1));    
    bl1(2, rep)= abs(Q4(:,2)'*V(:,2));
    bl1(3, rep)= abs(Q4(:,3)'*V(:,3));
    
    list1 = find(V(:,1:3));
    list2 = find(Q4(:,1:3));
    list3 = find(~V(:,1:3));
    
    
    for kindex =1:length(list1)       
        if ismember(list1(kindex), list2)==1 % hit
            truepos(2)=truepos(2)+1;
        else
            falseneg(2) = falseneg(2)+1;
        end
    end
    
    % if a true zero is identified as nonzero, fp
    for kindex = 1:length(list3)
        if ismember(list3(kindex), list2)==1
            falsepos(2)= falsepos(2)+1;
        end
    end
    
    if (isequal(find(Q4(:,1)), find(V(:,1)))&& isequal(find(Q4(:,2)), find(V(:,2)))&& isequal(find(Q4(:,3)), find(V(:,3))))
        recovery1 = recovery1+1;
    end
    if min(bl1(:,rep)>0.99)
        success1 = success1+1;
    end
    
    % TPower
    A = c;
    flag = 1;
    Q4=zeros(p,m);
    for i=1:m
        x0 = randn(p,1);
        x0 = x0/norm(x0);
        
        x1 = tp(x0, p, A);
%         x2 = tp(x1, 4*k, A);
%         x3 = tp(x2, 2*k, A);
        x4 = tp(x1, k, A);

        acc = abs(x4'*V(:,i));
        Q4(:,i)=x4;
        bl2(i,rep)=acc;
        P = eye(p)-x4*x4';

        A = P*A*P;
       
        if (~isequal(find(x4), find(V(:,i))))
            flag=0;
        end
    end
    if min(bl2(:,rep)>0.99)
        success2 = success2+1;
    end
    if(flag==1)
       recovery2 = recovery2+1;
    end
    
    list1 = find(V(:,1:3));
    list2 = find(Q4(:,1:3));
    list3 = find(~V(:,1:3));
    
    
    for kindex =1:length(list1)       
        if ismember(list1(kindex), list2)==1 % hit
            truepos(3)=truepos(3)+1;
        else
            falseneg(3) = falseneg(3)+1;
        end
    end
    
    % if a true zero is identified as nonzero, fp
    for kindex = 1:length(list3)
        if ismember(list3(kindex), list2)==1
            falsepos(3)= falsepos(3)+1;
        end
    end
    
   % TOrth
    A = c;
    Q0 = randn(p,3);
    [Q0, R] = qr(Q0, 0);
    Q1 = torth(Q0, [p,p,p], A);
%     Q2 = torth(Q1, [4*k, 4*k, 4*k], A);
%     Q3 = torth(Q2, [2*k, 2*k, 2*k], A);
    Q4 = torth(Q1, [k, k, k], A);

    
    bl3(1, rep)= abs(Q4(:,1)'*V(:,1));    
    bl3(2, rep)= abs(Q4(:,2)'*V(:,2));
    bl3(3, rep)= abs(Q4(:,3)'*V(:,3));
    
    list1 = find(V(:,1:3));
    list2 = find(Q4(:,1:3));
    list3 = find(~V(:,1:3));
    
    
    for kindex =1:length(list1)       
        if ismember(list1(kindex), list2)==1 % hit
            truepos(4)=truepos(4)+1;
        else
            falseneg(4) = falseneg(4)+1;
        end
    end
    
    % if a true zero is identified as nonzero, fp
    for kindex = 1:length(list3)
        if ismember(list3(kindex), list2)==1
            falsepos(4)= falsepos(4)+1;
        end
    end
    
    if (isequal(find(Q4(:,1)), find(V(:,1)))&& isequal(find(Q4(:,2)), find(V(:,2)))&& isequal(find(Q4(:,3)), find(V(:,3))))
        recovery3 = recovery3+1;
    end
    if min(bl3(:,rep)>0.99)
        success3 = success3+1;
    end
    
    % TOrthT
    A = c;
    Q0 = randn(p,3);
    [Q0, R] = qr(Q0, 0);
    Q1 = TOrthT(Q0, [p,p,p], A);
%     Q2 = TOrthT(Q1, [4*k, 4*k, 4*k], A);
%     Q3 = TOrthT(Q2, [2*k, 2*k, 2*k], A);
    Q4 = TOrthT(Q1, [k, k, k], A);
    
    bl4(1, rep)= abs(Q4(:,1)'*V(:,1));    
    bl4(2, rep)= abs(Q4(:,2)'*V(:,2));
    bl4(3, rep)= abs(Q4(:,3)'*V(:,3));
    
    
    list1 = find(V(:,1:3));
    list2 = find(Q4(:,1:3));
    list3 = find(~V(:,1:3));
    
    
    for kindex =1:length(list1)       
        if ismember(list1(kindex), list2)==1 % hit
            truepos(5)=truepos(5)+1;
        else
            falseneg(5) = falseneg(5)+1;
        end
    end
    
    % if a true zero is identified as nonzero, fp
    for kindex = 1:length(list3)
        if ismember(list3(kindex), list2)==1
            falsepos(5)= falsepos(5)+1;
        end
    end
    
    if (isequal(find(Q4(:,1)), find(V(:,1)))&& isequal(find(Q4(:,2)), find(V(:,2)))&& isequal(find(Q4(:,3)), find(V(:,3))))
        recovery4 = recovery4+1;
    end
    if min(bl4(:,rep)>0.99)
        success4 = success4+1;
    end
    
end
mean(bl0(1,:))
mean(bl0(2,:))
mean(bl0(3,:))

success0
recovery0

mean(bl1(1,:))
mean(bl1(2,:))
mean(bl1(3,:))

success1
recovery1

mean(bl2(1,:))
mean(bl2(2,:))
mean(bl2(3,:))

success2
recovery2

mean(bl3(1,:))
mean(bl3(2,:))
mean(bl3(3,:))

success3
recovery3

mean(bl4(1,:))
mean(bl4(2,:))
mean(bl4(3,:))

success4
recovery4


for i=1:5
    Fscore = truepos(i)/(truepos(i)+0.5*(falseneg(i)+falsepos(i)))
end