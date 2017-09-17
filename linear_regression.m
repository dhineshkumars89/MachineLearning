function output_args = linear_regression(filepath,degree,lambda)
%degree=1
%lambda=0
M=degree+1; %degree since it starts from 1 , appending 1
filepath='D:\CSE 6363 Machine Learning\Assignment\Assignment 4\sample_data1.txt';
delimiterIn = ' ';
A = importdata(filepath,delimiterIn);
X=A(1:20,1);
T=A(1:20,2);

sz=size(X);
W = zeros(M,1);
sigma = zeros(sz(1),M);
I = eye(M);
weights = zeros(3,1);

for m = 1:M
    for i = 1:sz(1)
        sigma(i,m) = (X(i)).^(m-1);
    end
end

weights = double(inv((lambda*I) + (vpa(sigma' * sigma))) * sigma' * T);
%class(weights)
%size(weights)
%set zero for w2 if it doesn't exist
try
  weights(3);
catch
  weights(3)=0;
end


fprintf('w0=%.4f\n', weights(1));
fprintf('w1=%.4f\n', weights(2));
fprintf('w2=%.4f\n', weights(3));

end

