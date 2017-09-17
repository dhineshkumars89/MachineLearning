function neural_network(training_file,test_file,layers,units_per_layer,rounds)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%input
% training_file='D:\CSE 6363 Machine Learning\Assignment\Assignment 2\pendigits_training.txt';
% test_file='D:\CSE 6363 Machine Learning\Assignment\Assignment 2\pendigits_test.txt';
% layers=2;
% units_per_layer=20;
% rounds=50;

delimiterIn = ' ';
train = importdata(training_file,delimiterIn);
test = importdata(test_file,delimiterIn);

%learning rate
learn_rate=1;
%-0.05+rand(4,4)*(0.05+0.05)  Wji

%classes
classes = unique(train(:,end));

%number of classes
[num_classes,coli] = size(classes);

%size of data set
[tr,tc] = size(train);
[sr,sc] = size(test);
d=(tc-1);

%maximum value in matrix train
cmax=max(max(train));

%number of units
if(layers>2)
    U = (tc-1)+1+num_classes+((layers-2)*units_per_layer);
else
    U = (tc-1)+1+num_classes;
end

%normalising the value using maximum value
X = zeros(tr,tc-1);
for i = 1:tr
    for j = 1:(tc-1)
        X(i,j) = (train(i,j)/cmax);
    end
end

%cln = d+1 + ((layers-2) * units_per_layer);
%rw = ((layers-2) * units_per_layer) + num_classes;
rxw = d+1 + ((layers-2) * units_per_layer) + num_classes;
xmin = -0.05;
xmax = 0.05;
%rand('seed', 100);
W = xmin+rand(rxw,rxw)*(xmax-xmin);
%W = ones(rxw,rxw)*0.05;
W_org = W;


for r = 1:rounds
    %disp(r);
    %round r learning_rate=learning_rate * power(0.98,r-1)
    learning_rate = learn_rate * power(0.98,r-1);
    T = zeros(tr,U);
    error = 0;
    
    for n = 1:tr
        A = zeros(U,1);
        G = zeros(U,1);
        Z = zeros(U,1);
        Z(1)=1;
        for j = 1:(tc-1)
            Z(j+1) = X(n,j) ; %Xnj
        end
        
        %set tn of particular class to 1
        curr_classindex = find(classes==train(n,end),1);
        T(n,U-num_classes+curr_classindex) = 1;
        
        
        for l = 2:layers
            
            if((l==2) && (l~=layers))
                input_units = d+1;
                output_units = units_per_layer;
                input_offset = 0;
                output_offset = d+1; 
            elseif(l==layers)
                if(l==2)
                    input_units = d+1;
                    input_offset = 0;
                else
                    input_units = units_per_layer;
                    input_offset = d+1+(units_per_layer * (l-3));
                end
                output_units = num_classes;
                output_offset = d+1+(units_per_layer * (l-2));
            else
                input_units = units_per_layer;
                output_units = units_per_layer;
                input_offset = d+1+(units_per_layer * (l-2));
                output_offset = d+1+(units_per_layer * (l-2));
            end
            outstart = output_offset+1;
            outstop = output_offset+output_units;
            for j = outstart:outstop    %handle output units case
                temp = W(j,1)*Z(1);
                
                instart = input_offset+1;
                instop = input_offset+input_units;
                for i = instart:instop
                    if(i~=1)
                        temp = temp + (W(j,i) * Z(i));
                    end
                end
                A(j)=temp;
                Z(j) = 1/(1+(exp(-A(j))));
                
            end            
        end
        %updation of weight
        for l = layers:-1:2
            if(l==layers)
                out_offset = d + 1 + (units_per_layer * (l-2));
                if(l==2)
                    in_offset = out_offset - (d+1);
                else
                    in_offset = out_offset - units_per_layer;
                end
                for j = out_offset+1:out_offset+num_classes    %handle output units case
                    G(j) = (Z(j)-T(n,j)) * Z(j) * (1-Z(j));
                    W(j,1) = W(j,1) - (learning_rate * G(j) * Z(1));
                    for i = in_offset+1:out_offset
                        if(i~=1)
                            W(j,i) = W(j,i) - (learning_rate * G(j) * Z(i));
                        end
                    end
                end
            else
                ind_end = d + 1 + (units_per_layer * (l-2));

                if(l==2)
                    in_offset = 0;
                else
                    in_offset = ind_end - units_per_layer;
                    
                end
                units = units_per_layer;
                lm=l+1;
                if(lm==layers)
                    curr_output = d + 1 + (units_per_layer * (lm-2));
                    curr_units = num_classes;
                else
                    curr_output = d + 1 + (units_per_layer * (lm-2));
                    curr_units = units_per_layer;
                end
                for j = ind_end+1:ind_end+units    
                    Gof = 0;
                    for u = curr_output+1:curr_output+curr_units
                        Gof = Gof + (G(u) * W(u,j));
                    end
                    G(j) = (Gof) * Z(j) * (1-Z(j));
                    W(j,1) = W(j,1) - (learning_rate * G(j) * Z(1));
                    for i = in_offset+1:ind_end
                        if(i~=1)
                            W(j,i) = W(j,i) - (learning_rate * G(j) * Z(i));
                        end
                    end
                end
            end

        end
        error = error + (0.5*(power(T(n)-Z,2)));
    end
    
    if(r~=1)
            %fprintf('error=%6.4f\n',abs(diff(error-old_error)));
            if(abs(diff(error-old_error)) < -1)
                %disp('BREAK');
            end
    end
    old_error = error;
end


wdiff = W-W_org;
%disp(W-W_org);

%classification
%normalising the value using maximum value
tst = zeros(sr,sc-1);
total_accuracy=0;
cmax=max(max(test));
for i = 1:sr
    for j = 1:(sc-1)
        tst(i,j) = (test(i,j)/cmax);
    end
end

for n = 1:sr
    A = zeros(U,1);
    Z = zeros(U,1);
    Z(1)=1;
    for j = 1:(sc-1)
        Z(j+1) = tst(n,j) ; 
    end

    for l = 2:layers
        if((l==2) && (l~=layers))
            input_units = d+1;
            output_units = units_per_layer;
            input_offset = 0;
            output_offset = d+1; 
        elseif(l==layers)
            if(l==2)
                input_units = d+1;
                input_offset = 0;
            else
                input_units = units_per_layer;
                input_offset = d+1+(units_per_layer * (l-3));
            end
            output_units = num_classes;
            output_offset = d+1+(units_per_layer * (l-2));
        else
            input_units = units_per_layer;
            output_units = units_per_layer;
            input_offset = d+1+(units_per_layer * (l-2));
            output_offset = d+1+(units_per_layer * (l-2));
        end
        for j = output_offset+1:output_offset+output_units    
            A(j) = W(j,1)*Z(1);
            for i = input_offset+1:input_offset+input_units
                if(i~=1)
                    A(j) = A(j) + (W(j,i) * Z(i));
                end
            end
            Z(j) = 1/(1+(exp(-A(j))));
        end
    end
    ary = Z(end-num_classes+1:end);
    cls_idx = [];
    max_val=0;
    for e = 1:num_classes
        if(max_val < ary(e))
            max_val = ary(e);
            idx = e;
            cls_idx=[];
        elseif(max_val == ary(e))
            cls_idx = [cls_idx, e];
        end    
    end
    %[probability index] = max(ary);
    probability = max_val;
    if(isempty(cls_idx))
        predicted_class = classes(idx);
    end
    %predicted_class = classes(index);
    true_class=test(n,end);
    if(predicted_class == true_class)
        accuracy=1;
    elseif(~isempty(cls_idx))
        for f = 1:length(cls_idx)
            pred_class = classes(cls_idx(f));
            if(pred_class==true_class)
               accuracy=1/length(cls_idx);
            else
               accuracy=0; 
            end
        end
    else
        accuracy=0;
    end
    
    total_accuracy = total_accuracy + accuracy;
    fprintf('ID=%5d, predicted=%3d, true=%3d, accuracy=%4.2f\n', n-1, predicted_class, true_class, accuracy);
end
classification_accuracy = total_accuracy/sr;
fprintf('classification accuracy=%6.4f\n', classification_accuracy);

end

    