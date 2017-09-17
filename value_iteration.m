function value_iteration( environment_file, non_terminal_reward, gam, K )

%environment_file = 'D:\CSE 6363 Machine Learning\Assignment\Assignment 13\environment2.txt';
%gam = 1;
defval = non_terminal_reward;
%K = 20;

realmx = realmax;

%intialization
R = getData(environment_file,realmx);

%R = [[0,0,0,1.0];[0,-realmx,0,-1.0]; [0,0,0,0]];
[N, M] = size(R);
maxl = max(max(R));

UT = zeros(N,M);
%marking block
for i = 1:N
        for j = 1:M
            if(R(i,j) == 0)
                R(i,j) = defval;
                UT(i,j) = defval;       
            else
                %do nothing 
                UT(i,j) = R(i,j);
            end
        end
end
%condition = true;
count = 0;
while(count <= K)
    U = UT;
    %g = 0;
    for i = 1:N
        for j = 1:M
           if((R(i,j) == maxl) || (R(i,j) == -maxl) || (R(i,j) == -realmx))
               %do nothing to terminal states
           else
               %evaluate all possible directions        
               UT(i,j) = R(i,j) + gam * valfunc(U,i,j,realmx);
           end
        end
    end
    count = count + 1;
end
%fprintf('Output:\n');
for i = 1:N
    strg = '';
    for j = 1:M
        vaal = UT(i,j);
        if(vaal == -realmx)
            vaal = 0.0;
        end
        if (j == M)
            %fprintf('%6.3f',UT(i,j))
            txt = sprintf('%6.3f',vaal);
            strg = strcat(strg,txt);
        else
            txt = sprintf('%6.3f,',vaal);
            strg = strcat(strg,txt);
        end
    end
    fprintf('%s\n',strg)
end

end

function maxval = valfunc(UT,i,j,realmx)
   %update the maxval, if new val is greater than max val
   %up 1 %down 2 %right 3 %left 4
   val = zeros(1,4);
   val(1,1) = (0.8 * getUT(UT,i+1,j,i,j,realmx)) + (0.1 * getUT(UT,i,j-1,i,j,realmx)) + (0.1 * getUT(UT,i,j+1,i,j,realmx));
   val(1,2) = (0.8 * getUT(UT,i-1,j,i,j,realmx)) + (0.1 * getUT(UT,i,j-1,i,j,realmx)) + (0.1 * getUT(UT,i,j+1,i,j,realmx));
   val(1,3) = (0.8 * getUT(UT,i,j+1,i,j,realmx)) + (0.1 * getUT(UT,i+1,j,i,j,realmx)) + (0.1 * getUT(UT,i-1,j,i,j,realmx));
   val(1,4) = (0.8 * getUT(UT,i,j-1,i,j,realmx)) + (0.1 * getUT(UT,i+1,j,i,j,realmx)) + (0.1 * getUT(UT,i-1,j,i,j,realmx));
   maxval = max(val);
end

function res = getUT(UT,a,b,oi,oj,realmx)
    
    [N, M] = size(UT);
    %check banging wall
    if(((a<1) || (b<1)) || ((a>N) || (b>M)))
        res = UT(oi,oj);
    else
       %check for block
       if(UT(a,b) == -realmx)
           %do nothing
           res = UT(oi,oj);
       else
           %check for proper cell
           res = UT(a,b);
       end
    end
    
    
end

function R = getData(environment_file,realmx)
    delimiterIn = ',';
    fid = fopen(environment_file);
    tline = fgetl(fid);
    mainArr = [];
    while ischar(tline)
       arr = strsplit(tline,delimiterIn);
       tline = fgetl(fid);
       mainArr = [mainArr; arr];
    end

    [p,q] = size(mainArr);
    R = zeros(p,q);
    for x = 1:p
        for y = 1:q
            if(mainArr{x,y} == '.')
                R(x,y) = 0;
            elseif(mainArr{x,y} == 'X')
                R(x,y) = -realmx;
            else
                R(x,y) = str2double(mainArr(x,y));
            end
        end
    end
end

