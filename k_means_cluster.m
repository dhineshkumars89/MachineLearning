function k_means_cluster(data_file,k,iterations)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% k_means_cluster yeast_test.txt 2 5
% k_means_cluster yeast_test.txt 3 5

%data_file = 'D:\CSE 6363 Machine Learning\Assignment\Assignment 2\yeast_test.txt';
%k=3;
%iterations=5;


%intialization
delimiterIn = ' ';
data_set = importdata(data_file,delimiterIn);

[tr tc] = size(data_set);

trow = randperm(tr);

%data selection
clusters = java.util.HashMap();
for g = 1:k
    list = java.util.ArrayList();
    if g==1
        offset=0;
    else
        offset=endval;
    end
    start = offset+1;
    if(g==k)
        endval = tr;
    else
        endval = offset+floor(tr/k);
    end
    for h = start:endval
        %list.add(data_set(h,1:tc-1));
        list.add(data_set(trow(h),1:tc-1));
    end
    clusters.put(g,list);
end


for i = 1:iterations+1
    
    meanhash = java.util.HashMap();

    %Mean function
    for g = 1:k
        val = cusmean(clusters.get(g),tc-1);
        %val
        meanhash.put(g,val);
    end

    %Error Function
    newhash = java.util.HashMap();
    for g = 1:k
        list = java.util.ArrayList();
        newhash.put(g,list);
    end

    total = 0;
    for g = 1:k
        lst = clusters.get(g);
        sz = lst.size();
        for j = 1:sz
            curval = lst.get(j-1);
            curdis = euclideanFucn(curval,meanhash.get(g));
            min = intmax;
            for b = 1:k
                eucDis = euclideanFucn(curval,meanhash.get(b));
                if(min > eucDis)
                    min = eucDis;
                    label = b;
                end
            end

            lsts = newhash.get(label);
            lsts.add(curval);
            newhash.put(label,lsts);

            total = total + curdis;
        end
    end
    
    if(i==1)
        fprintf('After initialization: error = %.4f\n', total);
    else
        fprintf('After iteration %d: error = %.4f\n', i-1, total);
    end

    clusters = newhash;

    %do iteration
    %re-assign elements to clusters
    %re-calculated cluster mean
    %calculate error
end

% for v = 1:tr
%     
%     list.add(data_set(v,:));
% end
%vertcat()

end

function dis = euclideanFucn(X,Y)
    sz = size(X,1);
    summ = 0;
    for i = 1:sz
        val = power((X(i,1) - Y(i,1)),2);
        summ = summ + val;
    end
    dis = sqrt(summ);
end

function val = cusmean(list,col)
val = zeros(1,col);
sz = list.size();
mat = zeros(sz,col);
for y = 1:sz
    mat(y,:) = list.get(y-1);
end

for cl = 1:col
    val(1,cl) = mean(mat(:,cl));
end

end
