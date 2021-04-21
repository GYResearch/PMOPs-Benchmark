function [Population] = TransferFun(pop,Union)
nPop=numel(pop);
PopObjs = Union.objs;
[NS,~]  = size(PopObjs);   
TMP = [];
Next = false(1,NS);
for i=1:nPop   
    TMP = [TMP;pop(i).Cost];
end
[~,d] =intersect(Union.objs,TMP,'rows');
Next(d) = 1;
Population = Union(Next);
end