function svd_power(data_file,M,iterations)
% data_file = 'D:\CSE 6363 Machine Learning\Assignment\Assignment 9\input1.txt';
% M = 4;
% iterations = 100;

delimiterIn = ' ';
train = importdata(data_file,delimiterIn);
[tr, tc] = size(train);

AAT = (train * train');
[ar, ac] = size(AAT);

D=ac;
X = zeros(ar,M,D);
U = zeros(D,M);
S = zeros(M,M);
for i = 1:ar
    X(i,1,:) = AAT(i,1:D)';
end
for d = 1:M
    %1)Covariance Matrix
    %(xn-x_mean)
    tran = squeeze(X(:,d,:));
    A = tran;
    %cov(tran);
    
    bk=ones(D,1);
    for k = 1:iterations
        normAbk = norm(A * bk);
        bk = (A * bk)/normAbk;
    end
    
    for i = 1:ar
        val = squeeze(X(i,d,:));
        X(i,d+1,:) = val - (bk' * val * bk);
    end

    U(:,d) = bk;
    S(d,d) = sqrt(normAbk);

end


ATA = (train' * train);
[ar, ac] = size(ATA);

D=ac;
X = zeros(ar,M,D);
V = zeros(D,M);
for i = 1:ar
    X(i,1,:) = ATA(i,1:D)';
end
for d = 1:M
    tran = squeeze(X(:,d,:));
    A = tran;
    %cov(tran);
    
    bk=ones(D,1);
    for k = 1:iterations
        bk = (A * bk)/norm(A * bk);
    end
    
    for i = 1:ar
        val = squeeze(X(i,d,:));
        X(i,d+1,:) = val - (bk' * val * bk);
    end

    V(:,d) = bk;

end
fprintf('Matrix U: \n');
printMatrix(U);
fprintf('Matrix S: \n');
printMatrix(S);
fprintf('Matrix V: \n');
printMatrix(V);
Output =  U * S * V';
fprintf('Reconstruction (U*S*V''): \n');
printMatrix(Output);
end

function printMatrix(data)
[tr, tc] = size(data);
    for i = 1:tr
        fprintf('  Row   %3d: ',i);
        for j = 1:tc
            fprintf('%8.4f ',data(i,j));
        end
        fprintf('\n');
    end
end
function A = covarianceMatrix(data)
    [tr, tc] = size(data);
    D=tc;
    Meean = zeros(tr,D);
    for di = 1:D
        x_mean = mean(data(:,di));
        for i = 1:tr
            Meean(i,di) = data(i,di) - x_mean;
        end
    end
    A = zeros(D,D);
    %summation covariance matrix
    A = Meean(1,:)' * Meean(1,:); 
    for i = 2:tr
        A = A + (Meean(i,:)' * Meean(i,:)); 
    end
    A = A/tr;
end