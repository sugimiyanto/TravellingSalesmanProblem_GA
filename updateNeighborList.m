function [ updated_neighborList ] = updateNeighborList(in_neighborList, removeElement)
updated_neighborList=containers.Map('KeyType','int32','ValueType','any');
    for i=1:length(in_neighborList)
        tmp = inf(1,length(in_neighborList(i)));
        tmp = in_neighborList(i);
        for j=1:length(tmp)
            tmp = tmp(tmp~=removeElement);
        end
        updated_neighborList(i) = tmp;
    end
end