function svec = solveSubpartADMM(D, phi, J, v0,knownInds,unknownInds)

verify = 1;

Dphi = D*phi;
Dt = D';
nV = size(Dt,2);
A = [ones(1,nV); speye(nV)];
B = [zeros(1,nV); -speye(nV)];
c = [v0; zeros(nV,1)];
q = rand(1+nV, 1);
svec = sparse(rand(nV,1));
s2vec = sparse(rand(nV,1));
s2vecpre = s2vec;
rho = rand*100+100;
rhomin = .001;
rhomax = 1e4;
e_abs = 1e-10;
e_rel = 10e-10;
tau_incr = 2;
tau_decr = 2;
mu = 10;

converged = 0; counter = 1;
while ~converged
    %% dual update
    q = q + rho*(A*svec+B*s2vec-c);
    
    %% solve svec
    K = (Dphi*Dphi').*(Dt'*Dt);
    h = sum((Dt'*J).*Dphi,2);
    Aprime = 2*K + rho * A'*A;
    bprime = (2*h' - q'*A - rho*(B*s2vec-c)'*A)';
    svec = Aprime\bprime;
    svec(knownInds,:) = 1;
    svec(unknownInds,:) = Aprime(unknownInds,unknownInds)\(bprime(unknownInds,:) - sum(Aprime(knownInds,unknownInds))');
   
    if verify
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
    s2vec = svec - q0/rho;
    s2vec(s2vec>1) = 1;
    s2vec(s2vec<0) = 0;
    
    if verify
        cvx_begin
            cvx_precision best
            cvx_solver mosek
            variable conductances2(nV,1);
            minimize transpose(q)*(A*svec+B*conductances2-c) + rho*pow_pos(norm(A*svec+B*conductances2-c),2)
            subject to
                conductances2 >= 0;
                conductances2 <= 1;
        cvx_end
        q'*(A*svec+B*s2vec-c) + rho*pow_pos(norm(A*svec+B*s2vec-c),2)
    end

    %% check convergence stuff
    obj0(counter) = norm(Dt*diag(svec)*Dphi-J,'fro');
        
    e_pri = e_abs*sqrt(numel(c)) + e_rel*max([norm(A*svec),norm(B*s2vec),norm(c)]);
    e_dual = e_abs*sqrt(size(B,2)) + e_rel*norm(A'*q);
    resid_primal = norm(A*svec+B*s2vec-c);
    resid_dual = rho*norm(A'*B*(s2vec-s2vecpre));
    rprimal(:,counter)=resid_primal; rdual(:,counter)=resid_dual;
    if (resid_primal <= e_pri && resid_dual < e_dual) || (max(resid_primal,resid_dual)<1e-10)
        converged = true;
    end

    if resid_primal > mu*resid_dual && rho <= rhomax
        rho = rho * tau_incr;
    elseif resid_dual > mu*resid_primal && rho >= rhomin
        rho = rho / tau_decr;
    end

    s2vecpre = s2vec;
    counter = counter + 1;
    fprintf('   ADMM: %d obj:%f rp:%g, rd:%g\n',counter,obj0(end),resid_primal,resid_dual);
end



end
