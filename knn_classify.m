function knn_classify (training_file,test_file,k)

%input
% training_file='D:\CSE 6363 Machine Learning\Assignment\Assignment 2\pendigits_training.txt';
% test_file='D:\CSE 6363 Machine Learning\Assignment\Assignment 2\pendigits_test.txt';
% k=1;

delimiterIn = ' ';
train = importdata(training_file,delimiterIn);
test = importdata(test_file,delimiterIn);

[traindata,testdata] = normalise(train,test);
%testdata = normalise(test);

%finding minimum distance train data

[tr tc] = size(traindata);
[tsr tsc] = size(testdata);
classes = unique(traindata(:,end));
totalacc=0;

%testdatum = testdata(2,:);
for r = 1:tsr
    testdatum = testdata(r,:);
    %minimum = 1000000;
    result = zeros(tr,2);
    for i=1:tr
        comp = euclideandist(testdatum,traindata(i,:));
        %fprintf('%8.6f - %d \n', comp,traindata(i,end)); 
        result(i,1) = comp; 
        result(i,2) = traindata(i,end);
    end
    output = sortrows(result,1);
    
    %mapObj = containers.Map;
    mapObj = containers.Map('KeyType','int64','ValueType','double');
    for s = 1:k
        match = output(s,:);
        if(isKey(mapObj,match(2)))
            mapObj(match(2)) = mapObj(match(2)) + 1;
        else
            mapObj(match(2)) = 1;
        end
%         if(match(2) == crtclass)
%             count = count + 1;
%             predicted_class = match(2);
%         end
    end
    acc = 0;
    crtclass = testdata(r,end);
    count = 0;
    clscount = 0;
    valuess = cell2mat(values(mapObj));
    keyss = cell2mat(keys(mapObj));
    %loop thro mapObj
    for stri = 1:length(keyss)
        strng = keyss(stri);
        if((mapObj(strng) == max(valuess)) && (length(unique(valuess)) == length(valuess)))
            if(strng == crtclass)
                acc = 1;
            else
                acc = 0;
            end
            predicted_class = strng;
            break;
        elseif(mapObj(strng) == max(valuess))
            count = count + 1;
            predicted_class = strng;
            if(strng == crtclass)
                clscount = clscount+1;
            end
        else
            
        end
    end
    if(count>0)
        if(clscount>0)
            acc = 1/count;
        else
            acc = 0;
        end
    end
    totalacc = totalacc + acc;
    fprintf('ID=%5d, predicted=%3d, true=%3d, accuracy=%4.2f  \n', r-1, predicted_class, crtclass, acc);
end

fprintf('classification accuracy=%6.4f\n', totalacc/tsr);

end

function dist = euclideandist(testdatum,traindatum)
summ = 0;
[tr tc] = size(traindatum);
col = tc-1;
for j = 1:col
    resu = power((testdatum(j) - traindatum(j)),2);
    summ = summ + resu;
end
dist=sqrt(summ);
end

function [data, tdata] = normalise(data,tdata)
[tr tc] = size(data);

%obj = zeros(tr,tc);
MeanD = zeros(tc-1);
stdD = zeros(tc-1);
for j = 1:tc-1
    MeanD(j) = mean(data(:,j));
    Summation = 0;
    for i = 1:tr
        res = power(abs(data(i,j) - MeanD(j)),2);
        Summation = Summation + res;
    end
    stdD(j) = sqrt(Summation/tr);
%     stdD(j)= std(data(:,j));
end

for i = 1:tr
    for j = 1:tc-1
        data(i,j) = (data(i,j)-MeanD(j))/stdD(j);
    end
end

[tsr, tsc] = size(tdata);
for i = 1:tsr
    for j = 1:tsc-1
        tdata(i,j) = (tdata(i,j)-MeanD(j))/stdD(j);
    end
end

end