function [ neighborList ] = generateNeighborList( parent1, parent2 )
% parent1=[4 1 3 5 7 6 2];
% parent2=[7 4 6 1 3 2 5];
    len = length(parent1);
    neighborList=containers.Map('KeyType','int32','ValueType','any');

    sorted1=sort(parent1);
    neighborParent1=inf(len,2);
    sorted2=sort(parent2);
    neighborParent2=inf(len,2);

    % generating parents neighbor of each node
    for i=1:len
        [~,idx]=find(parent1==sorted1(i));
        if idx == 1
            neighborParent1(i,:) = [parent1(len) parent1(2)];
        elseif idx == len
            neighborParent1(i,:) = [parent1(len-1) parent1(1)];
        else
            neighborParent1(i,:) = [parent1(idx-1) parent1(idx+1)];
        end
    end
    for i=1:len
        [~,idx]=find(parent2==sorted2(i));
        if idx == 1
            neighborParent2(i,:) = [parent2(len) parent2(2)];
        elseif idx == len
            neighborParent2(i,:) = [parent2(len-1) parent2(1)];
        else
            neighborParent2(i,:) = [parent2(idx-1) parent2(idx+1)];
        end
    end

    % combining neighbor to get all parent's neighbors of each node
    for i=1:len
        tmp = union(neighborParent1(i,:), neighborParent2(i,:));
%         neighborList(i,1:length(tmp)) = tmp;
        neighborList(i) = tmp;
    end
end

