% Evaluates polynomial(s) on 2D meshrid

% INPUT: cp      :  a vector of coefficients of the polynomial OR a cell  array of vectors of coefficients (indexed in Yalmip basis)
%        X1, X2  :  meshrid

% OUTPUT: p_val  :  values of the polynomial on the meshgrid OR a cell array of values of the polynomials on the meshgrid
function p_val = polyval_grid2D(cp,X1,X2)

ISCELL = 1;
if(~iscell(cp))
    ISCELL = 0;
    cp = {cp};
end


% Determine the largest total degree of the polynomials
nk = 0;
for k = 1:numel(cp)
    nk = max(nk,numel(cp{k}));
end
d = 0;
while(nchoosek(2+d,2) < nk)
    d = d+1;
end



X1pows = cell(1,d+1);
X2pows = cell(1,d+1);

X1pows{1} = ones(size(X1)); X2pows{1} = ones(size(X1));
for i = 1:d
    X1pows{i+1} = X1pows{i}.*X1;
    X2pows{i+1} = X2pows{i}.*X2;
end

pows = monpowers(2,d);

p_val = num2cell(zeros(1,numel(cp)));
for i = 1:size(pows,1)
    temp = X1pows{pows(i,1)+1}.*X2pows{pows(i,2)+1};
    
    for k = 1:numel(cp)
        if(i > numel(cp{k})) % Takes care of non-equal degrees of p's
            mult = 0;
        else
            mult = cp{k}(i);
        end
        p_val{k} = p_val{k} + mult*temp;
    end
end


if(ISCELL == 0)
    p_val = p_val{1};
end


