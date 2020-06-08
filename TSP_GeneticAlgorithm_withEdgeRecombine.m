% 1. generate random population
% 2. parent selection / tournamentSize
% 3. Tournament
% 4. Crossover
% 5. Mutation
% 6. Offspring
% 7. Survivor selection
% 8. goto 1 until iteration size = n
% 9. termination will hapen after defined iteration size

% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
clear;
clc;
sizePop = 100; % define here
iterSize = 200; % define here to determine termination
sizeTournament = 6; % tournament size should not be greater than size of population
% load route file
    [FileName,PathName] = uigetfile('*.txt','Select the text file');
    direct=strcat(PathName,FileName); % get path file
    [matrix,delimiterOut] = importdata(direct); % get the matrix data
    [sizeRow, colSize] = size(matrix);

% 1. generate random population (random possibility of path)
    population = inf(sizePop, sizeRow); % initialize place for pop
    for i = 1:sizePop
        population(i,:) = randperm(sizeRow);
    end
    % population contains random path of tour
    
    % initialization for GA cycle

    optRoute=inf(1,sizeRow);
    offspring=inf(sizeRow+2,sizeRow); % 2 is 1 pair or size of parents
    fit = inf(1,sizePop);
    distancePop=inf(1,sizePop);
    parents = inf(2,sizeRow); % 2 is 1 pair or size of parents
    avgFit = inf(1,iterSize);
    fitLocalBest = inf(1,iterSize);
    target = inf;
    figure(1);
    title('Convergence of the best-fit solution each generation');
    xlabel('Generation');
    ylabel('fitness value');
    
    h = animatedline('Marker', '.', 'Color', 'r');
    axis([1 iterSize 1 2000]);

    
    % #####for dynamic parent size allocation
    % selectSizeParents=0; % should be even (pair)
    % determine how many population as parents
    %     tmp = ceil(sizePop/2);
    %     if mod(tmp,2) == 0
    %         selectSizeParents = tmp;
    %     else
    %         selectSizeParents = tmp+1;
    %     end
    % #####
    %disp('initial population')
    
% ************* GA cycle started here *********
    for iter = 1:iterSize
% 2. parent selection / tournamentSize
        % calculate the distance of each population (per row), distance of each city.
        % to do parent selection with tournament way

        % remove duplicate path if exist
        population = unique(population, 'rows', 'stable'); 
        distAndPopPair=containers.Map('KeyType','int32','ValueType','any');
        [sizePop, ~] = size(population);
        %disp(population);
        for p = 1:sizePop
           dist=0;
            for i = 1:sizeRow-1
                dist = dist + matrix(population(p,i),population(p,i+1));
                  % sum per population or per row. get the data from matrix
                  % of distance. to get total distance per population
            end
            % count dist between the end of pop with the beginning of pop
            dist = dist + matrix(population(p,sizeRow),population(p,1)); 
            
            % create pair of total distance and the chromosom.
            % total distance as the key
            distAndPopPair(dist) = population(p,:); 
            
            distancePop(p) = dist;
            % when storing new key-value pair, it will be sorted ascending
            % automatically
        end
  
% 3. Calculate the fitness value of each population
    minDist = min(distancePop);
    if iter==1
        target = ceil(minDist/4);
    end
    
    for i=1:sizePop
          fit(1,i)=distancePop(1,i)-target; % delta value as the fitness
          % the less fitness value means the closest to target
    end

% 4. Tournament
        match = inf(1,sizeTournament);
        winner = inf(1,sizeTournament);
        % generate random participants
        for i=1:sizePop
            for j=1:sizeTournament
                randPart = ceil(rand()*sizePop);
                match(j) = fit(randPart);
            end
            winner(i) = min(match); % min of fitness value
            % the less fitness value, the closest to optimal solution of
            % current population. selected as parents
        end

        % count win
        countWin=containers.Map('KeyType','double','ValueType','int32');
        for i=1:sizePop
            if isKey(countWin,winner(i)) == 0
                countWin(winner(i)) = 1;
            else
                countWin(winner(i)) = countWin(winner(i))+1;
            end
        end
       
        % sort based on count
        k = keys(countWin);
        v = values(countWin);
        countNfit = inf(length(countWin), 2); % 2 here is column count and column fit
        sortCountNfit = inf(length(countWin), 2); % 2 here is column count and column fit
        for i=1:length(countWin)
            countNfit(i,1) = v{i};
            countNfit(i,2) = k{i};
        end
        sortCountNfit = countNfit;
        sortCountNfit = sortrows(sortCountNfit, -1); % descending to get highest count
         
        % get parents. assumed only 1 pair of parents
        idx=1;
        [sizeWinner,~] = size(sortCountNfit);

        iterParents=inf;
        if sizeWinner < 2
            iterParents=1;
        else
            iterParents = sizeWinner-1; % -1 depends on size of parents (2)
        end
        
        numParents = sizeWinner-iterParents+1;
        parents = inf(numParents,sizeRow);
        p=inf(1,numParents); % 2 is number of parents
        for i=1:numParents
            p(idx)=sortCountNfit(i,2); % take the fit, fit is in column 2
            idx=idx+1;
        end

        % get the chromosome of parents. assumed only 1 pair of parents
        for i=1:numParents
            [~,idxChromosome] = find(fit==p(i));
            % idxChromosome to get the index of chromosome in distAndPopPair

            distKey = distancePop(idxChromosome);
            distKey = unique(distKey, 'stable');
            parents(i,:) = distAndPopPair(distKey);
        end
     
 % 4. Crossover
    [sizeParents,~] = size(parents);
        
    % Edge recombination crossover
    crossover=0;
    crossProb=.5;
    if sizeParents == 2
        crossover=1;
        neighborList=containers.Map('KeyType','int32','ValueType','any');
        parent1 = parents(1,:);
        parent2 = parents(2,:);
        child = inf(1,sizeRow);
        neighborList = generateNeighborList(parent1, parent2);
        temp=[];
        % produce child
        randSelectFirst = ceil(rand()*2);
        if randSelectFirst == 1
            child(1) = parent1(1);
            temp=parent1;
            neighborList = updateNeighborList(neighborList, child(1));
        else
            child(1) = parent2(1);
            temp=parent2;
            neighborList = updateNeighborList(neighborList, child(1));
        end
        
        % steps
        for i=2:length(parent1)
            neighbor = [];
            if length(neighborList(child(i-1))) > 0
                neighbor = neighborList(child(i-1));
            else
                restCollection=[];
                for x=1:length(neighborList)
                    % here collecting all rest nodes which aren't included
                    % in child
                    if length(neighborList(x)) > 0
                        op = neighborList(x);
                        for p=1:length(op)
                            restCollection(p) = op(p);
                        end
                    end
                end
                restCollection = unique(restCollection, 'stable');
                if length(restCollection)
                    child(i) = restCollection(1);
                    neighborList = updateNeighborList(neighborList, child(i));
                end
            end
            
            %if the neighbor more than 1
            if length(neighbor) > 1
                n_neighbor=[];
                for y=1:length(neighbor)
                    n_neighbor(y) = length(neighborList(neighbor(y)));
                end
                n_neighbor = unique(n_neighbor, 'stable');
                idex=inf;
                if n_neighbor ==1
                    if isequal(neighbor(1),temp(i))
                        idex = 2;
                    else
                        idex = 1;
                    end
                else
                    [~,least] = min(n_neighbor);
                    idex = least;
                end
                
                child(i)=neighbor(idex);
                neighborList = updateNeighborList(neighborList, child(i));
            elseif length(neighbor) == 1
                child(i) = neighbor(1);
                neighborList = updateNeighborList(neighborList, child(i));
            end
        end
    end
    randCross = rand;
    if randCross > crossProb
        crossover = 1;
    else
        crossover = 0;
    end
    
%5. Mutation
        %mutation against selected parents to generate child
        mutRatio=.05;
        if crossover == 0
            [sizeChild,~] = size(parents);
        else
            [sizeChild,~] = size(child);
        end
        childs = inf(sizeChild,sizeRow);
        if crossover == 0
            childs = parents;
        else
            childs = child;
        end
        for i=1:sizeChild % 2 is 1 pair or size of parents
            idxSwappedAllele = ceil(rand(1,2)*sizeRow);
            % generate random value between 0 and 1. matrix 1x2
            % rand value*sizeRow to generate val between 1-sizeRow
            x = idxSwappedAllele(1);
            y = idxSwappedAllele(2);
            
            % not all chromosome will be mutated. so, generating
            % probability to determine whether it will be mutated or not
            mutRand = rand;
            if(mutRand > mutRatio)
                % mutate, a certain point of data row -> colum and vice versa
                childs(i,[y x]) = childs(i,[x y]);
            end
        end
       
%  6. Offspring
        offspring = [population;childs];
        % parents here are new generation population after crossover and
        % mutation
        
     
%  7. Survivor selection % select only best fitness as much as population size 
        % 'Fitness based' survivor selection
        % to keep the population size remain the same and preventing 'out of memory'
        % remove duplicate path if exist. with unique, it will sort the value
        offspring = unique(offspring, 'rows', 'stable'); 
        [sizeOffspring,~]=size(offspring);
        distNoffspring = inf(sizeOffspring,sizeRow+1);
 
        for j = 1:sizeOffspring
            dist=0;
            for i = 1:sizeRow-1
                dist = dist + matrix(offspring(j,i),offspring(j,i+1));
                % sum per offSpring or per row. get the data from matrix
                % of distance. to get total distance per offSpring
            end
            % count dist between the end of offspring with the beginning of offSpring
            dist = dist + matrix(offspring(j,sizeRow),offspring(j,1)); 
            
            % create pair of total distance and the chromosom.
            % total distance as the key
            distNoffspring(j,:) = [dist offspring(j,:)];
        end

    %  Calculate the fitness value of each offspring
%         maxDis = max(distNoffspring(:,1));
%         minDis = min(distNoffspring(:,1));
        fitOfofspring = inf(1,sizeOffspring);
        for i=1:sizeOffspring
%             fitOfofspring(1,i)=(maxDis - distNoffspring(i,1)+0.000001) / (maxDis-minDis+0.00000001);
             fitOfofspring(1,i) = distNoffspring(i,1)-target;
        end
%         fitOfofspring = log2(fitOfofspring);
        %disp('fitness of offspring')
%         fitOfofspring = fitOfofspring+1; % for 0.0000 problem
        %disp(fitOfofspring);
     
        % average fit
        su=0;
%         avgFit = inf(1,length(fitOfofspring));
        for i=1:length(fitOfofspring)
            su = su + (fitOfofspring(i));
        end
        av = su/length(fitOfofspring);
        avgFit(1,iter) = av;

        % replace the dist with fitness value to be sorted
        fitNoffspring = inf(sizeOffspring,sizeRow+1);
        fitNoffspring = distNoffspring;
        fitNoffspring(:,1) = fitOfofspring;
        sortedOffspring = sortrows(fitNoffspring); % ascending to get the lowest fit

    % filter offspring based on the best fitnes as much as initial population size
        sortedOffspring = sortedOffspring(:,2:end); % 2 here is for removing fitness value in column 1
        [sizeSorted, ~] = size(sortedOffspring);
        population = inf(sizeSorted-sizePop, sizeRow);
        idx=1;
        limitIdx = sizeSorted-sizePop;
        if limitIdx == 0
            limitIdx = 1;
        else
            limitIdx = limitIdx+1;
        end
        
        for i=sizeSorted:-1:limitIdx
            population(idx,:) = sortedOffspring(idx,:);
            idx=idx+1;
        end
        
        tmp = population(1,:);
        [row,~]=find(ismember(fitNoffspring(:,2:end),tmp,'rows'));
        localBest = fitNoffspring(row,1);
        fitLocalBest(1,iter) = localBest;
%         disp(localBest);
        
%         hold all;
%         plot(iter,localBest,'--*r');
%         hold off;
%         drawnow;
       
        addpoints(h,iter,localBest);
        drawnow
        
%         pause(0.05);
%         disp('survival selection')
%         disp(population);
%         disp('============================================================');
% 8. goto 2 until iteration size = n
    end
    fprintf('terminated with %i iteration\n', iterSize);
    
    %get the optimal route and show the route and total distance
    optRoute = population(1,:); % the top row of sorted population based on fitness (minimization)
    disp('optimal route');
    disp(optRoute);
    disp('with total distance');
    [row,~]=find(ismember(distNoffspring(:,2:end),optRoute,'rows'));
    disp(distNoffspring(row,1));
    
    figure(2);
    plot(1:iterSize,fitLocalBest,'-r',1:iterSize,avgFit,'-b');
    legend('local best','average');%,'Location','southwest');
    xlabel('Generation');
    ylabel('Fitness value');
    disp('with average fitness');
    disp(avgFit(end));
    disp('best fitness');
    disp(fitLocalBest(end));
    ylim([min(fitLocalBest)-100 max(avgFit)+100]);
    title('Average fitness and fitness with edge recombine crossover');
%     
%     final figure of shortest path TSP
%     figure(3);
%     distance_to_position(FileName, optRoute);
