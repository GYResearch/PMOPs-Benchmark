function [pop] = TransferFun(Union,pindex)

 PopObjs = Union.objs;
 PopDecs = Union.decs;
 [nPop,M]  = size(PopObjs);
 
empty_individual.Position=[];
empty_individual.Cost=[];
empty_individual.Rank=[];
empty_individual.DominationSet=[];
empty_individual.DominatedCount=[];
empty_individual.CrowdingDistance=[];
empty_individual.Associate = [];

pop=repmat(empty_individual,nPop,1);

for i=1:nPop   
    pop(i).Position= PopDecs(i,:);
    pop(i).Cost= PopObjs(i,:);
    pop(i).Associate = pindex(i);
end



end