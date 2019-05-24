function [svec, rho, q] = solveSubpartADMM(D, phi, J, v0,knownInds,unknownInds, startConductances, startq)

if nargin == 0
    nverts = randi(50)+10;
    nedges = randi(50)+10;
    nmeas = randi(50)+10;
    D = rand(nedges,nverts)-.5;
    phi = rand(nverts,nmeas);
    knownInds = randsample(nedges,10);
    unknownInds = find(~ismember([1:nedges],knownInds));
    condGT = rand(nedges,1); condGT(knownInds)=1;
    v0 = sum(condGT);
    J = D'*diag(condGT)*D*phi + 10*rand(size(phi));
    
    [svecMosek, E] = solveByCVX(D, phi, J, v0,knownInds,unknownInds);
    verifySub = 0;
    verify = 1;
else
    verifySub = 0;
    verify = 0;
end
Dphi = D*phi;
Dt = D';
nV = size(Dt,2);
A = [ones(1,nV); speye(nV)];
B = [zeros(1,nV); -speye(nV)];
c = [v0; zeros(nV,1)];
if exist('startq','var')
    q = startq;
else
    q = rand(1+nV, 1);
end
if exist('startConductances','var')
    svec = startConductances;
else
    svec = sparse(rand(nV,1));
end
s2vec = sparse(rand(nV,1));
s2vecpre = s2vec;
rho = 132.531;
rhomin = .001;
rhomax = 1e4;
e_abs = 1e-1;
e_rel = 1e-1;
tau_incr = 2;
tau_decr = 2;
mu = 10;
K = (Dphi*Dphi').*(Dt'*Dt);
h = sum((Dt'*J).*Dphi,2);
svec(knownInds,:) = 1;
rho_updated = true;

converged = 0; counter = 1;
while ~converged
    %% dual update
    q = q + rho*(A*svec+B*s2vec-c);
    
    %% solve svec
    if rho_updated
        Aprime = 2*K + rho * A'*A;
        Aprime2 =  Aprime(unknownInds,unknownInds);
        %Aprime2 = 2*K(unknownInds,unknownInds) + rho * A(:,unknownInds)'*A(:,unknownInds); 
        [Lm, Dm, pm] = ldl(Aprime2, 'vector');
    end
    bprime = (2*h' - q'*A - rho*(B*s2vec-c)'*A)';
    bprime2 = (bprime(unknownInds,:) - sum(Aprime(knownInds,unknownInds))');
    clear x; x(pm,:) = Lm'\(Dm\(Lm\(bprime2(pm,:))));
    svec(unknownInds,:) = x; % Aprime2\bprime2
   
    if verifySub 
        cvx_begin
            cvx_solver mosek
            variable conductances(nV,1);
            physFeas = pow_pos(norm(D' * diag(sparse(conductances)) * D ...
                    *phi - J,'fro'),2);
            minimize physFeas + rho*pow_pos(norm(A*conductances+B*s2vec-c),2) ...
             + transpose(q)*(A*conductances+B*s2vec-c) 
            subject to
                 conductances(knownInds) == 1
        cvx_end
        norm(D' * diag(sparse(svec)) * D * phi - J,'fro')^2 + ...
                q'*(A*svec + B*s2vec - c) + rho*pow_pos(norm(A * svec + B * s2vec - c),2)
    end
    
    %% solve s2vec
    q0 = q(2:end);
    s2vec = svec + q0/rho;
    s2vec(s2vec>1) = 1;
    s2vec(s2vec<0) = 0;
    
    if verifySub 
        cvx_begin
            cvx_precision best
            cvx_solver mosek
            variable conductances2(nV,1);
            %minimize transpose(q)*(A*svec+B*conductances2-c) + rho/2*pow_pos(norm(A*svec+B*conductances2-c),2)
            minimize transpose(q0)*(svec-conductances2) + rho/2*pow_pos(norm(svec - conductances2),2)
            subject to
                conductances2 >= 0;
                conductances2 <= 1;
        cvx_end
        transpose(q0)*(svec-s2vec) + rho/2*pow_pos(norm(svec - s2vec),2)
        % (q'*(A*svec+B*s2vec-c) + rho*pow_pos(norm(A*svec+B*s2vec-c),2))
    end

    %% check convergence stuff
    obj0(counter) = norm(Dt*diag(svec)*Dphi-J,'fro')^2;
        
    e_pri = e_abs*sqrt(numel(c)) + e_rel*max([norm(A*svec),norm(B*s2vec),norm(c)]);
    e_dual = e_abs*sqrt(size(B,2)) + e_rel*norm(A'*q);
    resid_primal = norm(A*svec+B*s2vec-c);
    resid_dual = rho*norm(A'*B*(s2vec-s2vecpre));
    rprimal(:,counter)=resid_primal; rdual(:,counter)=resid_dual;
    if (resid_primal <= e_pri && resid_dual < e_dual) || (max(resid_primal,resid_dual)<1e-10)
        converged = true;
    end

    rho_updated = false;
    if resid_primal > mu*resid_dual && rho <= rhomax
        rho = rho * tau_incr;
        rho_updated = true;
    elseif resid_dual > mu*resid_primal && rho >= rhomin
        rho = rho / tau_decr;
        rho_updated = true;
    end

    s2vecpre = s2vec;
    counter = counter + 1;
    if rho_updated
        fprintf('   ADMM: %d obj:%f rp:%g, rd:%g rho:%g (rhoupdated!)\n',counter,obj0(end),resid_primal,resid_dual,rho);
    else
        fprintf('   ADMM: %d obj:%f rp:%g, rd:%g rho:%g\n',counter,obj0(end),resid_primal,resid_dual,rho);
    end
end

if verify
    EADMM = norm(D' * diag(sparse(svec)) * D * phi - J,'fro')^2;
    EMOSEK = norm(D' * diag(sparse(svecMosek)) * D * phi - J,'fro')^2;

    assert(all(svec(knownInds)==1))
    assert(all(svec<=1+1e-3))
    assert(all(svec>=-1e-3))
    assert(norm(sum(svec)-v0)<=1e-3)
    [EADMM EMOSEK]
end
end

function [out, E] = solveByCVX(D, phi, J, v0,knownInds,unknownInds)
    cvx_begin
        cvx_solver mosek
        variable conductances0(size(D,1),1);
        physFeas = pow_pos(norm(D' * diag(sparse(conductances0)) * D * phi - J,'fro'),2);

        minimize physFeas;
        subject to
            conductances0 <= 1;
            conductances0 >= 0;
            sum(conductances0) == v0;
            conductances0(knownInds)==1
    cvx_end
    out = conductances0;
    E = pow_pos(norm(D' * diag(sparse(conductances0)) * D * phi - J,'fro'),2);
end
