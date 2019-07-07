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
iterSize = 1000; % define here to determine termination
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

        %fprintf('============== iteration %i ===================\n', iter);
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
        %disp('distance of initial population');
        %disp(distancePop);
%         
% 3. Calculate the fitness value of each population
%     distancePop = sort(distancePop);
%     maxDist = max(distancePop);
    minDist = min(distancePop);
    if iter==1
        target = ceil(minDist/4);
%         axis([1 iterSize target minDist]);
    end
    
    for i=1:sizePop
%         fit(1,i)=(maxDist - distancePop(1,i)+0.000001) / (maxDist-minDist+0.00000001);
          fit(1,i)=distancePop(1,i)-target; % delta value as the fitness
          % the less fitness value means the closest to target
    end
%     fit=log2(fit);
    
    %disp('fitness value of initial population')
    %disp(fit)
    % to address error reading of value 0.0000
%     fit = fit+1;
    %fprintf('the less fitness value, the closest to optimal solution of current population\n\n');
%     %fprintf('the fitness value for further is added by 1 to address error interpretation of value 0.0000\n');
% 
% 4. Tournament
        match = inf(1,sizeTournament);
        winner = inf(1,sizeTournament);
        % generate random participants
        for i=1:sizePop
            for j=1:sizeTournament
                randPart = ceil(rand()*sizePop);
                match(j) = fit(randPart);
            end
            %fprintf('match %i',i);
            %disp(match);
            winner(i) = min(match); % min of fitness value
            % the less fitness value, the closest to optimal solution of
            % current population. selected as parents
        end
        %disp('winner per match');
        %disp(winner);

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
        %disp('count winner')
        countNfit = inf(length(countWin), 2); % 2 here is column count and column fit
        sortCountNfit = inf(length(countWin), 2); % 2 here is column count and column fit
        for i=1:length(countWin)
            countNfit(i,1) = v{i};
            countNfit(i,2) = k{i};
        end
        sortCountNfit = countNfit;
        sortCountNfit = sortrows(sortCountNfit, -1); % descending to get highest count
        %disp(countNfit);
         
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
%         for i=sizeWinner:-1:iterParents % length(countNfit)-1 depends on size of parents
        for i=1:numParents
            p(idx)=sortCountNfit(i,2); % take the fit, fit is in column 2
            idx=idx+1;
        end

        %disp('dist and pop');
        %disp(keys(distAndPopPair));

        % get the chromosome of parents. assumed only 1 pair of parents
        for i=1:numParents
            [~,idxChromosome] = find(fit==p(i));
            % idxChromosome to get the index of chromosome in distAndPopPair

            distKey = distancePop(idxChromosome);
            distKey = unique(distKey, 'stable');
            parents(i,:) = distAndPopPair(distKey);
        end
%         disp('selected parents');
%         disp(parents)
     
 % 4. Crossover
    % Sub-tour exchange crossover
    [sizeParents,~] = size(parents);
    crossover=0;
    crossProb=.5;
    if sizeParents == 2 % 2 here is size of parents as defined (1 pair)
%         crossover=1;
        child = inf(6,sizeRow);
        parent1 = parents(1,:);
        parent2 = parents(2,:);
        for i=1:3
            child(i,:) = parent1;
            child(i+3,:) = parent2;
        end
        
        [found, idxA2, idxB2, idxA3, idxB3] = getCommonSubtour(parent1, parent2);
        if found == 1 % not found common genes
%             disp('found');
            if length(union(idxA2,idxA3)) < 5 % means only one of them (2 or 3 genes common)
                if length(unique(idxA3)) > 1 % will overwrite found 2 with 3 common genes
%                     disp('found 3');
                    child(1,idxA3) = parent2(1,idxB3);
                    child(2,idxB3) = parent1(1,idxA3);
                elseif length(unique(idxA2)) > 1 % if found 2 common genes
%                     disp('found 2');
                    child(1,idxA2) = parent2(1,idxB2);
                    child(2,idxB2) = parent1(1,idxA2);
                end
                
                for i=1:4
                    child(3,:)=[]; % remove 3rd - 6th child
                end
            else
%                 disp('found both of them');
                % 'if found both 2 and 3 genes'
                child(1,idxA2) = parent2(1,idxB2);
                child(2,idxA3) = parent2(1,idxB3);
                child(3,idxA2) = parent2(1,idxB2);
                child(3,idxA3) = parent2(1,idxB3);
                
                child(4,idxB2) = parent1(1,idxA2);
                child(5,idxB3) = parent1(1,idxA3);
                child(6,idxB2) = parent1(1,idxA2);
                child(6,idxB3) = parent1(1,idxA3);
            end
        end
        
        child = unique(child, 'rows', 'stable');
        
    end
    randCross = rand;
    if randCross > crossProb
        crossover = 1;
    else
        crossover = 0;
    end
    
        
    % Edge recombination crossover

    
%5. Mutation
        %mutation against selected parents to generate child
        mutRatio=.02;
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
        %disp('childs');
        %disp(childs);
       
%  6. Offspring
        offspring = [population;childs];
        % parents here are new generation population after crossover and
        % mutation
        
     
%  7. Survivor selection % select only best fitness as much as population size 
        % 'Fitness based' survivor selection
        % to keep the population size remain the same and preventing 'out of memory'
        % remove duplicate path if exist. with unique, it will sort the value
        offspring = unique(offspring, 'rows', 'stable'); 
        %disp('offspring after deduplication')
        %disp(offspring)
        [sizeOffspring,~]=size(offspring);
        distNoffspring = inf(sizeOffspring,sizeRow+1);
%         
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
        %disp('distance of offspring')
        %disp(distNoffspring);
% 
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
%         
        % average fit
        su=0;
%         avgFit = inf(1,length(fitOfofspring));
        for i=1:length(fitOfofspring)
            su = su + (fitOfofspring(i));
        end
        av = su/length(fitOfofspring);
        avgFit(1,iter) = av;
% 
        % replace the dist with fitness value to be sorted
        fitNoffspring = inf(sizeOffspring,sizeRow+1);
        fitNoffspring = distNoffspring;
        fitNoffspring(:,1) = fitOfofspring;
        sortedOffspring = sortrows(fitNoffspring); % ascending to get the lowest fit
        %disp('sorted offspring');
        %disp(sortedOffspring);
%         
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
%         for i=1:limitIdx
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
        
%         pause(0.000000000001);
%         disp('survival selection')
%         disp(population);
%         disp('============================================================');
%8. goto 2 until iteration size = n
    end
    fprintf('terminated with %i iteration\n', iterSize);
    
    %get the optimal route and show the route and total distance
    optRoute = population(1,:); % the top row of sorted population based on fitness (minimization)
    disp('optimal route');
    disp(optRoute);
    disp('with average fitness');
    [row,~]=find(ismember(distNoffspring(:,2:end),optRoute,'rows'));
    %disp(distNoffspring(row,1));
    figure(2);
    disp(avgFit(end));
    disp('best fitness');
    disp(fitLocalBest(end));
    plot(1:iterSize,fitLocalBest,'-r',1:iterSize,avgFit,'-b');
    legend('local best','average');%,'Location','southwest');
    xlabel('Generation');
    ylabel('Fitness value');
    ylim([min(fitLocalBest)-100 max(avgFit)+100]);
    title('Tournament size=6, Pc=0.5, Pm=0.02');
    
    % final figure of shortest path TSP
%     figure(3);
%     distance_to_position(FileName, optRoute);
