function logistic_regression(training_file,degree,test_file)
%training_file='D:\CSE 6363 Machine Learning\Assignment\Assignment 2\pendigits_training.txt';
%degree=1;
%test_file='D:\CSE 6363 Machine Learning\Assignment\Assignment 2\pendigits_test.txt';
delimiterIn = ' ';
train = importdata(training_file,delimiterIn);
test = importdata(test_file,delimiterIn);

%yeast 2nd degree
% train(:,6)=[];
% test(:,6)=[];
%yeast 2nd Degree

[tr,tc] = size(train);
[sr,sc] = size(test);

range = find(train(:,end)~=1);
for m = range(:,1)
    train(m,end)=0;
end

range = find(test(:,end)~=1);
for m = range(:,1)
    test(m,end)=0;
end

T=train(:,end);
lim = (degree*(tc-1))+1;
W=zeros(lim,1);
cross_entro_err=0;
cross_entro_err_old=100;

flag=1;
while(flag==1)
    
    bas=ones(lim,1);
    Y=zeros(tr,1);
    sigma=ones(tr,lim);
    for rw = 1:tr
        for m = 2:lim
            if(degree==1)
                bas(m,1)= train(rw,m-1);
            else
                if(mod(m,2)==0)
                    bas(m,1)= train(rw,round(m/2));
                else
                    bas(m,1)=bas(m-1,1)*bas(m-1,1); 
                end
                %disp(m+'---'+bas(m,1));
            end
            sigma(rw,m)=bas(m,1);
        end
        %sigma(rw)= bas';
        %sigmodial function
        y = W' * bas;
        Y(rw) = 1/(1+(exp(-y)));
        %disp(Y);
        cross_entro_err = cross_entro_err + T(rw)*log(Y(rw)) + (1-T(rw))*log(1-Y(rw));
    end
    
    %computer cross-entrophy error
   
    if(abs(cross_entro_err-cross_entro_err_old) < 0.001)
        flag=0;
    end
    
    %Calculate R 
    R=zeros(tr,tr);
    for i = 1:tr
        for j = 1:tr
            if(i==j)
                R(i,j)=Y(i,1) * (1-Y(i,1));
            end
        end
    end

    Wold=W;
    W=Wold-(pinv(sigma'*R*sigma) * sigma' * (Y-T));
    cross_entro_err_old = cross_entro_err;
    cross_entro_err=0;
    if(sumabs(W-Wold)<0.001)
        flag=0;
    end
    %flag=0;
end
for i = 1:length(W(:,1))
    %disp(W(i,1));
    fprintf('w%d=%.4f\n',i-1, W(i,1));
end


%Classification
class_accuracy=0;
Y=zeros(sr,1);
T=test(:,end);
lim = (degree*(tc-1))+1;
bas=ones(lim,1);
for rw = 1:sr
    for m = 2:lim
            if(degree==1)
                bas(m,1)= test(rw,m-1);
            else
                if(mod(m,2)==0)
                    bas(m,1)= test(rw,round(m/2));
                else
                    bas(m,1)=bas(m-1,1)*bas(m-1,1); 
                end
            end
    end
    y = W' * bas;
    Y(rw) = 1/(1+(exp(-y)));
    
    accuracy=0;
    if(Y(rw)>0.5)
        predicted = 1;
        if(predicted==T(rw,1))
            accuracy=1;
        end
    elseif(Y(rw)<0.5)
        predicted = 0;
        if(predicted==T(rw,1))
            accuracy=1;
        end
    else
        %clashing condition
        predicted=2;
        accuracy=0.5;
    end
    
    prob = (Y(rw).^T(rw,1)) * ((1-Y(rw)).^(1-T(rw,1)));
    if prob < 0.5
        prob= 1-prob;
    end
    class_accuracy=class_accuracy+accuracy;
    fprintf('ID=%5d, predicted=%3d, probability = %.4f, true=%3d, accuracy=%4.2f\n',rw-1, predicted, prob, T(rw,1), accuracy);
end
%disp(class_accuracy);
classification_accuracy = (class_accuracy./sr);
fprintf('classification accuracy=%6.4f\n', classification_accuracy);
end