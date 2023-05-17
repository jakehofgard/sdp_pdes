% Evaluates polynomials V1 and V1 with coefficients CV1 and CV2 of degree d on the 3D
% meshgrid given by X1, X2, X3
function [V1_val,V2_val] = eval_v1v2_grid3D(CV1,CV2,d,X1,X2,X3)

V1_val = zeros(size(X1));
V2_val = zeros(size(X1));

for k = 1:size(X1,3)
    if(numel(X1) > 300^3 && mod(k,50) == 0)
        fprintf('Evaluating v1 and v2: %f %% complete \n', 100 * k / size(X1,3))
    end
    
    X1pows = cell(1,d+1);
    X2pows = cell(1,d+1);
    X3pows = cell(1,d+1);
    
    X1pows{1} = ones(size(X1,1),size(X1,2)); X2pows{1} = ones(size(X1,1),size(X1,2)); X3pows{1} = ones(size(X1,1),size(X1,2));
    for i = 1:d
        X1pows{i+1} = X1pows{i}.*X1(:,:,k);
        X2pows{i+1} = X2pows{i}.*X2(:,:,k);
        X3pows{i+1} = X3pows{i}.*X3(:,:,k);
    end
    
    
    pows = monpowers(3,d);
    
    V1_val_k = 0;
    V2_val_k = 0;
    
    for i = 1:size(pows,1)
        temp = X1pows{pows(i,1)+1}.*X2pows{pows(i,2)+1}.*X3pows{pows(i,3)+1};
        V1_val_k = V1_val_k + CV1(i)*temp;
        V2_val_k = V2_val_k + CV2(i)*temp;
    end
    
    V1_val(:,:,k) = V1_val_k;
    V2_val(:,:,k) = V2_val_k;
end