function [Q] = torth(Q0, klist, A)
maxiter = 200;
Q = Q0;

for j=1:maxiter
    Z=A*Q;    
    for i=1:size(Q0,2)
        Z(:,i)=truncate_operator(Z(:,i), klist(i));        
    end

    [Qnew, R]=qr(Z,0);
    for i=1:size(Q,2)
         temp=truncate_operator(Qnew(:,i), klist(i));   
         Qnew(:,i)=temp/norm(temp);
    end
    
    if (norm(Qnew-Q)<1e-12)
        break
    end
    Q = Qnew;
end

end