function [Q] = ortht(Q0, klist, kfinal, A)
maxiter = 800;
Q = Q0;

for j=1:maxiter+1
    Z=A*Q; 
    
    for i=1:size(Q0,2)
        Z(:,i)=truncate_operator(Z(:,i), min(klist(i), size(A,2)));        
    end
    
    [Z, R]=qr(Z,0);


    for i=1:size(Q0,2)
        temp=truncate_operator(Z(:,i), min(klist(i), size(A,2)));        
        Qnew(:,i)=temp/norm(temp);
    end
    Qnew = full(Qnew);
    
    % reduce klist norm(Z(:,i)-Qnew(:,i))<1e-4 ||
                
   if (rem(j,150)==0 || norm(Qnew-Q)<1e-12)
        for i = 1:length(klist)
            if (klist(i)>kfinal(i))
                klist(i)= max(round(klist(i)/2), kfinal(i));
            end
        end

   end
        
   
%     if (norm(Qnew-Q)<1e-12)
%         break
%     end
    Q = Qnew;
end
%     for i=1:size(Q0,2)
%         temp=truncate_operator(Q(:,i), 10);        
%         Q(:,i)=temp/norm(temp);
%     end

end

