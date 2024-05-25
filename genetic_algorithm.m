prompt = "Enter the size of the population: ";
population = input(prompt);
prompt = "Enter the length of the chromosome: ";
length = input(prompt);
prompt = "Enter the crossover probability[0-1]: ";
crossP = input(prompt); %Taking the crossover probability.
prompt = "Enter the mutation probability[0-1]: ";
mutP = input(prompt); %Taking the mutation probability.
fitarr = zeros(population,1); %Initialize fitness array.
max = 2^length -1;
min = 0;

%Initializing the current population & mating pool.

currentpop = rand(population,length);
currentpop = round(currentpop);

while (min ~= -5)
    %Creating the mating pool.
    matingPool = zeros(population,length);
    for j=1:population
        parent1 = randi(population);
        parent2 = randi(population);
        p1chromosome = num2str(currentpop(parent1,:));
        p2chromosome = num2str(currentpop(parent2,:));
        p1value = bin2dec(p1chromosome);
        p2value = bin2dec(p2chromosome);
        p1mapValue = mapping(p1value);
        p2mapValue = mapping(p2value);
        p1fitval = fitnessValue(p1mapValue);
        p2fitval = fitnessValue(p2mapValue);
        if p1fitval < p2fitval
            matingPool(j,:) = currentpop(parent2,:);
        else
            matingPool(j,:) = currentpop(parent1,:);
        end
    
    end


    %Performing crossover(single point) operation.
    for k=1:population
        prob = unifrnd(0,1);
        if prob <= crossP
            parent1 = randi(population);
            parent2 = randi(population);
            p1chromosome = currentpop(parent1,:);
            p2chromosome = currentpop(parent2,:);
            crossPoint = round(unifrnd(1,(length-1)));
            child1 = [p1chromosome(1:crossPoint),p2chromosome(crossPoint+1:end)];
            child2 = [p2chromosome(1:crossPoint),p1chromosome(crossPoint+1:end)];
            matingPool(parent1,:) = child1;
            matingPool(parent2,:) = child2;
        end
    end


    %Performing mutation operation in mating pool.
    for l=1:population
        for m=1:length
            r = rand;
            if r < mutP
                matingPool(l,m) = 1-matingPool(l,m);
            end
        end
    end


    %Calculating the fitness of the current population.
    for i=1:population
        string = num2str(currentpop(i,:));
        value = bin2dec(string);
        temp = mapping(value);
        fitval = fitnessValue(temp);
        fitarr(i,1) = fitval;
    end

    min=1000;
    for i=1:population
        value = fitarr(i,1);
        if value <= min
        min = value;
        end
    end

    %Calculating the fitness of the matingPool.
    for i=1:population
        string = num2str(matingPool(i,:));
        value = bin2dec(string);
        temp = mapping(value);
        fitval = fitnessValue(temp);
        fitarr(i,1) = fitval;
    end

    min=1000;
    index=0;
    for i=1:population
        value = fitarr(i,1);
        if value <= min
           min = value;
           index=i;
        end
    end
    fprintf("The minimum value of the function is: %d\n",min);
    currentpop = matingPool;
end
fprintf("\nThe minimum value of the function is: %d",min);
string = num2str(matingPool(index,:));
value = bin2dec(string);
fprintf("\nThe value of x for which the minimum is obtained: %d\n",value);


%Bringing the range of the value in range 0 - 10.
function [mapVal] = mapping(value)
    mapVal = (10/255)*value;
end

%Computing fitness value.
function [fitVal]= fitnessValue(val)
    fitVal = (val^2-5);
end

