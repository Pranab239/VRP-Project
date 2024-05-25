% Optimizing the Vehicle Routing Problem (One depot)using GA.

% Taking the number of instances to be created.
prompt = "Enter the size of the population: "; 
population = input(prompt);
prompt = "Enter the Number of vehicles: ";
vehicle = input(prompt);
prompt = "Enter the Crossover Probability: ";
crossP = input(prompt);
prompt = "Enter the Mutation Probability: ";
mutP = input(prompt);
prompt = "Enter the vehicle capacity: ";
vCapacity = input(prompt);
fitarr = zeros(2*population,1,'double');
count = 0;
convarr=[];
c=[];
d=[];
flag = 1;
%Demand Matrix of the citiies.
demandMat = [
1 0 
2 19 
3 21 
4 6 
5 19 
6 7 
7 12 
8 16 
9 6 
10 16 
11 8 
12 14 
13 21 
14 16 
15 3 
16 22 
17 18 
18 19 
19 1 
20 24 
21 8 
22 12 
23 4 
24 8 
25 24 
26 24 
27 2 
28 20 
29 15 
30 2 
31 14 
32 9 
];
%co-ordinates of the cities. first co-ordinate is of the depot.
distmat = [
 1 82 76
 2 96 44
 3 50 5
 4 49 8
 5 13 7
 6 29 89
 7 58 30
 8 84 39
 9 14 24
 10 2 39
 11 3 82
 12 5 10
 13 98 52
 14 84 25
 15 61 59
 16 1 65
 17 88 51
 18 91 2
 19 19 32
 20 93 3
 21 50 93
 22 98 14
 23 5 42
 24 42 9
 25 61 62
 26 9 97
 27 80 55
 28 57 69
 29 23 15
 30 20 70
 31 85 60
 32 98 5
];
tic
[row, ~] = size(distmat);
adjMatrix = zeros(row, row, 'double');
adjMatrix = findAdj(distmat,adjMatrix);
adjMatrix(1,1) = Inf;
fprintf("The Adjacency Matrix: \n");
disp(adjMatrix);

% Initializing the current population.

[row,col] = size(adjMatrix);
len = col + (vehicle - 1);
currentpop = zeros(population, len);
temp = ones(1,vehicle - 1);

%Implementing sweep algorithm.
chromosome = sweepAlgo(adjMatrix, len, vCapacity, demandMat, col, vehicle);
checkZero = 0;
for i = 1:len
    if(chromosome(i) == 0)
        checkZero = 1;
        break;
    end
end
% disp(chromosome);
% fitness = fitnessValue(chromosome, len, adjMatrix);
% fitness


% %generating the chromosomes and check if valid.
% %m is #rows and n is #cols.
% 
[m,n] = size(currentpop);
if(checkZero == 0)
    currentpop(1,:) = chromosome;
    i = 2;
else
    i = 1;
end
while(i <= m)
    %cities are randomized and no of vechicles are added.
    cities = randperm(col);
    currentpop(i,:) = [cities,temp];
    
    %the below code is to make the currentpop shuffle.
    currentpop(i,:) = shuffle(currentpop(i,:));
    
    %if the last value is 1 in chromosome.
    if(currentpop(i,n) == 1 && currentpop(i,1) ~= 1)
        currentpop(i,:) = removeLastOne((currentpop(i,:)));
    end
    
    %if last and first both are 1 in chromosome.
    if((currentpop(i,n) == 1) && (currentpop(i,1) == 1))
        currentpop(i,:) = removeLastAndFirstOne(currentpop(i,:), n);
    end
    
    %if the first value if not 1 in chromosome.
    if(currentpop(i,1) ~= 1)
        currentpop(i,:) = makeFirstPlaceOne(currentpop(i,:), n);
    end
    
    %if two 1 occurs simultenously.
    if(checkConsecutiveOne(currentpop(i,:), n) == 1)
        currentpop(i,:) = removeConsecutiveOne(currentpop(i,:), n);
    end
    
    %if still one occurs at the end.
    if(currentpop(i,n) == 1)
        currentpop(i,:) = moveOneToMid(currentpop(i,:), n);
    end
    
    %checking if the chromosome satisfy the weight constraints.
    if(weightValidation(currentpop(i,:), demandMat, n, vCapacity) == 1)
        i = i + 1;
        continue;
    end
end
fprintf("The current population: \n");
disp(currentpop);


% Creating mating Pool.

matingPool = zeros(population,len);
while(flag ~= 0)
    prevParent = [];
    pplen = 0;
    for j=1:population
        parent1=1; 
        parent2=1;
        while(parent1 == parent2)
            parent1 = randi(population);
            parent2 = randi(population);
        end
        parent1
        parent2
        p1chromosome = currentpop(parent1,:);
        p2chromosome = currentpop(parent2,:);
        p1fitval = fitnessValue(p1chromosome,len,adjMatrix);
        p2fitval = fitnessValue(p2chromosome,len,adjMatrix);
        if p1fitval < p2fitval
            matingPool(j,:) = currentpop(parent1,:);
        else
            matingPool(j,:) = currentpop(parent2,:);
        end
     end
%     fprintf("The mating pool: \n");
%     disp(matingPool);
    
% Performing crossover operation.
    i = 1;
    while(i <= population)
        prob = unifrnd(0,1);
        if prob <= crossP
            parent1=0; 
            parent2=0;
            while(parent1==parent2)
                parent1 = randi(population);
                parent2 = randi(population);
            end 
            p1chromosome = matingPool(parent1,:);
            p2chromosome = matingPool(parent2,:);
            offspring1 = exchangeCross(p2chromosome, n, vehicle);
            if(checkConsecutiveOne(offspring1, n) == 1)
                offspring1 = removeConsecutiveOne(offspring1, n);
            end
            if(offspring1(n) == 1)
                offspring1 = removeLastAndFirstOne(offspring1, n);
            end
            offspring2 = exchangeCross(p1chromosome, n, vehicle);
            if(checkConsecutiveOne(offspring2, n) == 1)
                offspring2 = removeConsecutiveOne(offspring2, n);
            end
            if(offspring2(n) == 1)
                offspring2 = removeLastAndFirstOne(offspring2, n);
            end
            if(weightValidation(offspring1, demandMat, n, vCapacity) == 1 && weightValidation(offspring2, demandMat, n, vCapacity) == 1)
                matingPool(parent1,:) = offspring1;
                matingPool(parent2,:) = offspring2;
                i = i + 1;
                continue;
            else
                continue;
            end
        end
        i = i + 1;
    end
%     fprintf("Mating Pool after crossover: \n");
%     disp(matingPool);
    
% Performing Mutation operation.
    
    k = 1;
    while(k <= population)
        prob = unifrnd(0,1);
        if prob <= mutP
            parent = randi(population);
            p1chromosome = matingPool(parent,:);
            point1 = 1;
            point2 = 1;
            while(p1chromosome(point1) == 1 || p1chromosome(point2) == 1 || point1 == point2)
                point1 = round(unifrnd(2,len));
                point2 = round(unifrnd(2,len));
            end
            p1chromosome = swapArrayEl(p1chromosome, point1, point2);
            if(weightValidation(p1chromosome, demandMat, n, vCapacity) == 1)
                matingPool(parent,:) = p1chromosome;
                k = k + 1;
                continue;
            else 
                continue;
            end
        end
        i = i + 1;
    end
%     fprintf("MatingPool after the mutation: \n");
%     disp(matingPool);
    
    % Creating  a temporary population with current population and matingPool
    tempPool = currentpop;
    j=1;
    for i = (population+1):(2*population)
        tempPool(i,:) = matingPool(j,:);
        j=j+1;
    end
    
    % Computing the fitness value.
    for i=1:(2*population)
        fitarr(i,:)=fitnessValue(tempPool(i,:),len,adjMatrix);
    end
    count=count+1;
%     avg=findmean(fitarr,(2*population));
%     d=[d,avg];
%     c=[c,count];
    fprintf("The minimum distance is: %d(%d)\n",min(fitarr),count); 
    convarr = [convarr,min(fitarr)];
%     plot(c,convarr,c,d,'.');drawnow
%     title('Blue: Minimum            Green: Average');
%     xlabel("Generation");
%     ylabel("Min Dist")
    
    %converging the solution.
    if(count == 10)
        break;
    end
    if (count > 100)
        for i=count:-1:(count-100)
            if min(fitarr)==convarr(i)
                flag=0;
            else
                flag=1;
            end
        end
        if flag == 0
            break;
        end
    end
    temp = 0;
    for i=1:(population)
        [M,I] = min(fitarr);
        if(temp == 0)
            bestChromosome = I;
            bestfitness = M;
            temp = 1;
        end
        currentpop(i,:)=tempPool(I,:);
        fitarr(I,1) = Inf;
    end
    currentpop
end

fprintf("The optimal solution value is:\n");
disp(bestfitness);
fprintf("The optimal solution is:\n");
disp(tempPool(bestChromosome,:));
fprintf("No of generations taken: \n");
disp(count);

toc
currentpop


%all required user defined functions.


%function to remove 1 from the last position.
function[chromosome] = removeLastOne(chromosome)
    chromosome = flip(chromosome);
end


%function to remove if one occurs at the first place and last place also.
function[chromosome] = removeLastAndFirstOne(chromosome, n)
    replacedPos = randi(floor(n/2));
    count = 1;
    for j = (n-1):-1:1
        %swapping with the randomly choosen non one element from the last.
        if(chromosome(j) ~= 1 && count == replacedPos)
            chromosome = swapArrayEl(chromosome, j, n);
            break;
        end 
        if(chromosome(j) ~= 1)
            count = count + 1;
        end
    end
end


%function to make the first place 1.
function[chromosome] = makeFirstPlaceOne(chromosome, n)
    for j = 2:n
        if(chromosome(j) == 1)
            %swapping first one from left with the first element.
            chromosome = swapArrayEl(chromosome, j, 1);
            break;
         end
    end
end


%function to check consecutive ones.
function[flag] = checkConsecutiveOne(chromosome, n)
    flag = 0;
    for ind = 1:n-1
        if(chromosome(ind) == chromosome(ind + 1))
            flag = 1;
            break;
        end
    end
end


%function to remove two consecutive ones.
function[chromosome] = removeConsecutiveOne(chromosome, n)
    for k = 1:(n-2)
        if(chromosome(k) == 1 && chromosome(k+1) == 1)
            ind = k + 1;
            %if more than one 1 is present go to the 1st non-one value.
            while(ind < n && chromosome(ind) == 1)
                ind = ind + 1;
            end
            chromosome = swapArrayEl(chromosome, ind, k+1);
        end      
    end
end


%function to move last one to the mid.
function[chromosome] = moveOneToMid(chromosome, n)
    mid = floor(n/2);
    chromosome = swapArrayEl(chromosome, mid, n);
end


%function to exchange chromosome in crossover.
function[tempChromosome] = exchangeCross(p2chromosome, len, vehicle)
    count = vehicle;
    tempChromosome = zeros(1,len);
    tempChromosome(1) = 1;
    while(count > 1)
        pos = randi([3,len - 1], 1, 1);
        if(tempChromosome(pos) == 1)
            continue;
        end
        tempChromosome(pos) = 1;
        count = count - 1;
    end
    %all ones are positioned in the random places now.
    p1left = 2;
    p2left = 2;
    while(p1left <= len)
        if(p2chromosome(p2left) ~= 1 && tempChromosome(p1left) ~= 1)
            tempChromosome(p1left) = p2chromosome(p2left);
            p1left = p1left + 1;
            p2left = p2left + 1;
        elseif (tempChromosome(p1left) == 1)
            p1left = p1left + 1;
        elseif(p2chromosome(p2left) == 1)
            p2left = p2left + 1;
        end
    end
end


%function to swap elements in 1D array.
function[chromosome] = swapArrayEl(chromosome, point1, point2)
    temp = chromosome(point1);
    chromosome(point1) = chromosome(point2);
    chromosome(point2) = temp;
end


%function to calculate adjacency matrix.
function [adjMatrix] = findAdj(distmat, adjMatrix)
    [m,~]=size(distmat);
    r=1; 
    for i=1:m
        for j=1:m
            adjMatrix(i,j) = power(power((distmat(r,2)-distmat(j,2)),2)+power((distmat(r,3)-distmat(j,3)),2),.5);
        end
        r=r+1;
    end
end


%calculating fitness value of each chromosome.
function [fitval]=fitnessValue(chromosome,length,adjMatrix)
    fitval=0;
    for i=1:length-1
        fitval=fitval+adjMatrix(chromosome(i),chromosome(i+1));
    end
    fitval=fitval+adjMatrix(chromosome(length),chromosome(1));
end


%finding the mean of the array.
function [mean]=findmean(array,length)
    sum=0;
    for i=1:length
        sum = sum + array(i,1);
    end
    mean = sum/length;
end


%function to weight validation of chromosome.
function [flag] = weightValidation(chromosome, demandMat, n, vCapacity)
    demand = 0;
    for i = 2:n
        if(chromosome(i) ~= 1)
            demand = demand + demandMat(chromosome(i),2);
            if(demand > vCapacity)
                %vehicle capacity exceeds.
                flag = 0;
                return;
            end
        else
            %new car so initialize demand.
            demand = 0;
        end
    end
    flag = 1;
end


%function to make 1D array shuffle.
function [b]=shuffle(array)
    [~,n]=size(array);
    idx=randperm(n);
    b=array;
    b(1,idx)=array;
end
  

%Sweep algorithm to initialize a greedy chromosome.
function[chromosome] = sweepAlgo(adjMatrix, len, vCapacity, demandMat, col, vehicle)
    chromosome = zeros(1,len);
    %first place should be 1 always.
    chromosome(1) = 1;
    count = 1; %counting the no of ones in chromosome.
    i = 2;
    weight = 0;
    %nn = nearest neighbour.
    while(i <= len)
        nn = findMin(chromosome, adjMatrix, col, i);
        weight = weight + demandMat(nn,2);
        if(weight <= vCapacity)
            chromosome(i) = nn;
        else
            count = count + 1;
            if(count > vehicle) 
                break;
            end
            chromosome(i) = 1;
            weight = 0;
        end
        i = i + 1;
    end
end


%Function to calculate the minimum in an array except 0 given row and
%matrix & eliminate if already present in chromosome.

function[nn] = findMin(chromosome, adjMatrix, col, curlen)
    row = chromosome(curlen - 1);
    min = Inf;
    i = 1;
    while(i <= col)
        if(adjMatrix(row,i) ~= 0 && adjMatrix(row,i) < min && notExists(i, chromosome, curlen) == 1)
            min = adjMatrix(row,i);
            nn = i;
        end
        i = i + 1;
    end
end


%Function to check if city already exists in chromosome.

function [flag] = notExists(wanted, chromosome, curlen)
   flag = 1;
   for j = 1 : (curlen-1)
       if( wanted == chromosome(j))
           flag = 0;
       end
   end
end