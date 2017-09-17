function dtw_classify( training_file , test_file)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% delimiterIn = ' ';
% train = importdata(training_file,delimiterIn);
% test = importdata(test_file,delimiterIn);

% fileID = fopen(training_file,'r');
% %formatSpec = '%d %f';
% 'object ID:'
% content = fscanf(fileID,formatSpec);

%input
%training_file='D:\CSE 6363 Machine Learning\Assignment\Assignment 11\asl_training.txt';
%test_file='D:\CSE 6363 Machine Learning\Assignment\Assignment 11\asl_test.txt';

traindata = extractObj(training_file);
testdata = extractObj(test_file);

totalacc = 0;
%tstval = cell2mat(values(testdata));
tstkey = cell2mat(keys(testdata));
for i = 1:length(tstkey)
    strkey = tstkey(i);
    testdatum = testdata(strkey);
    tobjectID = testdatum.objectID;
    tclassID = str2double(testdatum.classID);
    costmap = containers.Map('KeyType','int64','ValueType','double');
    mincost=intmax;
    mincost_class=0;
    acc=0;
    %do against every object in training data
    trnkey = cell2mat(keys(traindata));
    for j = 1:length(trnkey)
       trnstrkey = trnkey(j); 
       traindatum = traindata(trnstrkey);
       cost = findCost(testdatum,traindatum);
       costmap(traindatum.objectID) = cost;
       if(cost<mincost)
           mincost = cost;
           mincost_class = str2double(char(traindatum.classID));
           %mincost_obj = traindatum.objectID;
       end
    end
%     disp(mincost)
%     disp(mincost_class);
    if(mincost_class==tclassID)
        acc = 1;
    else
        acc = 0;
    end
    fprintf('ID=%5d, predicted=%3d, true=%3d, accuracy=%4.2f, distance = %.2f\n', tobjectID, mincost_class, tclassID, acc, mincost);
    totalacc = totalacc + acc;
end


fprintf('classification accuracy=%6.4f\n', totalacc/length(tstkey));

end

function costi = findCost(tst,trn)
[a m] = size(tst.X);
[b n] = size(trn.X);
C = zeros(m,n);

C(1,1) = euclidean(tst.X(1,1),tst.Y(1,1),trn.X(1,1),trn.Y(1,1));
for i = 2:m
    C(i,1) = C(i-1,1) + euclidean(tst.X(1,i),tst.Y(1,i),trn.X(1,1),trn.Y(1,1));
end
for j = 2:n
    C(1,j) = C(1,j-1) + euclidean(tst.X(1,1),tst.Y(1,1),trn.X(1,j),trn.Y(1,j));
end

for i = 2:m
    for j = 2:n
        C(i,j) = min([C(i-1,j) ,C(i,j-1) , C(i-1,j-1)]) + euclidean(tst.X(1,i),tst.Y(1,i),trn.X(1,j),trn.Y(1,j));
    end
end
costi = C(m,n);
end

function dist = euclidean(x,y,a,b)
dist = sqrt(power((x - a),2) + power((y - b),2)); 
end

function mapObj = extractObj(file)
fid = fopen(file);
tline = fgetl(fid);
count = 0;
flag = 1;
X = [];
Y = [];
mapObj = containers.Map('KeyType','int64','ValueType','any');
while ischar(tline)
    if(contains(tline,'-------------------------------------------------'))
        %do something
        count=0;
        if(flag==1)
            flag=0;
        else
            %save X and Y from file
            data.X = X;
            data.Y = Y;
            %store value in map 
            mapObj(data.objectID) = data;
            clear data;
            X = [];
            Y = [];
        end
    else
        count=count+1;
    end
    if(count==1)
        objectID = regexp(tline,'object ID: (\d*)','tokens');
        %disp(objectID{:});
        data.objectID = str2double(objectID{1,1});
    elseif(count==2)
        classID = regexp(tline,'class label: (\d*)','tokens');
        %disp(classID{:});
        data.classID = char(classID{1,1});
    elseif(count==3)
        sign = regexp(tline,'sign meaning: (\S*)','tokens');
        %disp(sign{:});
        data.sign = char(sign{1,1});
    elseif(count > 5)
        %splitstring = regexp(tline,'\d*','split')
        strtxt = strtrim(tline);
        x = split(strtxt);
        %disp(x(1,:));
        %disp(x(2,:));
        X = [X, str2double(x(1,1))];
        Y = [Y, str2double(x(2,1))];
        %disp(y);
    else
        %empty line
    end
    
    tline = fgetl(fid);
end
%final object
data.X = X;
data.Y = Y;
%store value in map 
mapObj(data.objectID) = data;
fclose(fid);
end


        %find low cost and their corresponding class
%         crtclass = testdatum.classID;
%         coskey = cell2mat(keys(costmap));
%         cosval = cell2mat(values(costmap));
%         if((length(unique(cosval)) == length(cosval)))
%             [x idx]= min(cosval);
%             if(traindata(idx).classID == crtclass)
%                 acc = 1;
%             else
%                 acc = 0;
%             end
%             predicted_class = traindata(idx).classID;
%             break;
%         elseif(costmap(coskey) == max(valuess))
%             count = count + 1;
%             predicted_class = strng;
%             if(strng == crtclass)
%                 clscount = clscount+1;
%             end
%         else
%             
%         end