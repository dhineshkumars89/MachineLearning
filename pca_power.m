function pca_power(training_file,test_file,M,iterations)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% training_file='D:\CSE 6363 Machine Learning\Assignment\Assignment 2\yeast_training.txt';
% test_file='D:\CSE 6363 Machine Learning\Assignment\Assignment 2\yeast_test.txt';
% M=3;
% iterations=30;

delimiterIn = ' ';
train = importdata(training_file,delimiterIn);
test = importdata(test_file,delimiterIn);

[tr, tc] = size(train);

D=tc-1;
X = zeros(tr,M,D);
BK = zeros(M,D);
for i = 1:tr
    X(i,1,:) = train(i,1:D)';
end
for d = 1:M
    %1)Covariance Matrix
    %(xn-x_mean)
    tran = squeeze(X(:,d,:));
    A = covarianceMatrix(tran);
    
    %bk=ones(D,1);
    bk=rand(D,1);
    for k = 1:iterations
        bk = (A * bk)/norm(A * bk);
    end
    
    for i = 1:tr
        val = squeeze(X(i,d,:));
        X(i,d+1,:) = val - (bk' * val * bk);
    end

    BK(d,:) = bk;

end


%print
for i = 1:M
    fprintf('Eigenvector %3d \n',i);
    for j = 1:D
        fprintf('%3d: %.4f \n',j,BK(i,j));
    end
end

%test data projection
fprintf('\n');
[sr, sc] = size(test);
sd = sc-1;
for i = 1:sr
    dat = test(i,1:sd)';
    F = zeros(M,1);
    U = zeros(M,D);
    for p = 1:M
        U(p,:) = BK(p,:)';
        %F(p,1) = BK(p,:)' * dat;
    end
    F = U * dat;
    fprintf('Test object %3d \n',i-1);
    for j = 1:M
        fprintf('%3d: %.4f \n',j,F(j,1));
    end
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

