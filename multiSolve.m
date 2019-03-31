% solves many small linear systems in one go. Surely faster than
% for-looping, but kind of silly conceptually.
% As:nxn2xm 
% bs: nxm
% solutions: n2xm
% TIMING
% blockdiag version is faster if each subproblem is actually square.
% otherwise forloop is faster.

%%
function solutions = multiSolve(As,bs)

% if nargin == 0
%     n = 4; n2=3; m=1000000;
%     As = rand(n,n2,m)-.5;
%     bs = rand(n,m)-.5;
% end

n = size(As,1);
n2 = size(As,2);
m = size(As,3);
assert(size(bs,1)==n);
assert(size(bs,2)==m);

if n==n2
    A_n2m_n = reshape(permute(As,[1 3 2]),n2*m,n);
    II = repmat([1:n*m]',n2,1);
    Jt = reshape(repmat((0:(m-1))*n2,n,1),1,[]);
    JJ = reshape((repmat(Jt,n2,1)+[1:n2]')',[],1)';
    KK = A_n2m_n(:);
    bulkMat = sparse(II,JJ,KK,n*m,n2*m);
    bulkB = bs(:);
    solutions = reshape(bulkMat\bulkB,n2,m);
else
    % if it's not perfectly diagonal blocks, the sparse block diag version
    % is slower than just for looping.
    solutions = zeros(n2,m);
    for i = 1:m
        solutions(:,i) = As(:,:,i)\bs(:,i);
    end
end

% if nargin == 0
%     for i = 1:m
%         assert(norm(As(:,:,i)\bs(:,i)-solutions(:,i))<1e-4);
%     end
% end


end