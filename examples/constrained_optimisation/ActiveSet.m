function [v,step,f] = ActiveSet(P,q,AI,b,f0)
% script to find minimum of quadratic function with linear inequality
% constraints, using active sets
% func to be minimised is
%     f'P f + q'f
% constraints are
%     A'*f >= b
% Assume x0 is feasible
k = 1; f(:,1) = f0;v(1) = 0.5*f0'*P*f0 + f0'*q;step(1) = 0;
n = size(P,1);  % the number of unknowns
m = size(AI,2);  % the number of constraints 
Aset = find(AI'*f(:,k)-b >= 0); % active set
KeepGoing = true; tau = 1;
htol = 1e-6; % tolerance for step length
while(KeepGoing && k < 10)
    disp(['Iteration ',num2str(k),' active set ']); Aset
    na = size(Aset,1); % number of active constraints
    Iset = setdiff([1:m],Aset); % inactive set
    % find if active set is empty
    asz = size(Aset); nda = length(asz); 
    nsz = 1;
    for nn= 1:nda
        nsz = nsz*asz(nn);
    end
    if(nsz > 0)
        AA = [P,AI(:,Aset); AI(:,Aset)', zeros(na)];
        g = -q -P*f(:,k); % gradient of Phi(f)
        x = AA\[g;zeros(na,1)];
        h = x(1:n);lambda = x(n+1:end);
    else
        AA = P;
        g = -q -P*f(:,k); % gradient of Phi(f)
        x = AA\g;
        h = x(1:n);
        
    end
    if(norm(h) < htol)
        rset = find(lambda<0);      % this set is switched off
        if(size(rset,1) == 0)
            KeepGoing = false;
        else
            lmin = min(lambda); rset = find(lambda == lmin);
            Aset = setdiff(Aset,Aset(rset));  % this set is new active set
%            Iset = setdiff([1:m],Aset); % inactive set
        end
    end
    if KeepGoing
        na = length(Aset); % number of active constraints
        AA = [P,AI(:,Aset); AI(:,Aset)', zeros(na)];
        g = -q -P*f(:,k); % gradient of Phi(f)
        x = AA\[g;zeros(na,1)];
        h = x(1:n);lambda = x(n+1:end);
        if size(Iset,1) > 0
            proj = (b(Iset) - AI(:,Iset)'*f(:,k))./(AI(:,Iset)'*h); % orthogonal distances to inactive sets
            vset = find(proj >= 0); % only positive steps allowed.
            tau = min([1;proj(vset)]);
        else
            tau = 1;
        end
        if(tau < 1)
            p = find(tau==proj);
            Aset = [Aset;Iset(p)]; % add to active set
        end
    end
    f(:,k+1) = f(:,k) + tau*h;
    k = k+1;
    v(k) = 0.5*f(:,k)'*P*f(:,k) + f(:,k)'*q;
    step(k) = norm(h);
end
