function [ found, idxA2, idxB2, idxA3, idxB3 ] = getCommonSubtour( chromoA, chromoB )
% a=[12 14 8 6 4 7 1 3 4 0];
% b=[12 4 6 14 9 0 4 3 1 5];
tmp2=inf(1,2);
tmp2b=inf(1,2);
idx2=inf(1,2);
idx2b=inf(1,2);

tmp3=inf(1,3);
tmp3b=inf(1,3);
idx3=inf(1,3);
idx3b=inf(1,3);

% disp(a);
% disp(b);
found=0;
found2=0;
for i=1:length(chromoB)-1
    for j=1:length(chromoB)-1
        if isequal(chromoB(i:i+1), chromoA(j:j+1)) == 0
            if isequal(sort(chromoB(i:i+1)), sort(chromoA(j:j+1))) == 1
                found2=1;
                idx2=[j j+1];
                idx2b=[i i+1];
                tmp2=chromoA(j:j+1);
                tmp2b=chromoB(i:i+1);
                break;
            end
        end
    end
    if found2 == 1
        break;
    end
end
if found2 == 1
%     disp('found at chromoA');
%     disp(tmp2);
%     disp('found at chromoB');
%     disp(tmp2b);
%     disp('chromoA at index');
%     disp(idx2);
%     disp('chromoB at index');
%     disp(idx2b);
    idxA2 = idx2;
    idxB2 = idx2b;
else
    idxA2 = 0;
    idxB2 = 0;
end

found3=0;
for i=1:length(chromoB)-2
    for j=1:length(chromoB)-2
%         disp(b(i:i+2));
%         disp(a(j:j+2));
%         fprintf('\n');
        if isequal(chromoB(i:i+2), chromoA(j:j+2)) == 0
            if isequal(sort(chromoB(i:i+2)), sort(chromoA(j:j+2))) == 1
                found3=1;
                idx3=[j:j+2];
                idx3b=[i:i+2];
                tmp3b=chromoB(i:i+2);
                tmp3=chromoA(j:j+2);
                break;
            end
        end
    end
    if found3 == 1
        break;
    end
%     disp('=================');
end

if found3 == 1
%     disp('found at chromoA');
%     disp(tmp3);
%     disp('found at chromoB');
%     disp(tmp3b);
%     disp('chromoA at index');
%     disp(idx3);
%     disp('chromoB at index');
%     disp(idx3b);
    idxA3 = idx3;
    idxB3 = idx3b;
else
    idxA3 = 0;
    idxB3 = 0;
end

if found2 || found3
    found=1;
end
end

