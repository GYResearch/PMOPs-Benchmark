function Score = KIGD(PopObj,PF)
% <metric> <min>
% Knee-based Inverted generational distance
    PF = unique(PF,'rows');
    [n,~] = size(PopObj);
    if n>0
        Distance = min(pdist2(PF,PopObj),[],2);
        Score    = mean(Distance);
    end
end