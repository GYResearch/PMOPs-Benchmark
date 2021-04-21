function [Population,V] = LAEnvironmentalSelection(Union,V,N0,alpha)
% The environmental selection of LA-MOEA
%--------------------------------------------------------------------------

    PopObj = Union.objs;
    [N,M]  = size(PopObj);
    NV     = size(V,1);
    
    %% Translate the population
    PopObj = PopObj - repmat(min(PopObj,[],1),N,1);
    
    %% Associate each solution with one reference point
    %% Calculate the distance of each solution to each reference vector
    Cosine   = 1 - pdist2(PopObj,V,'cosine');
    Distance = repmat(sqrt(sum(PopObj.^2,2)),1,NV).*sqrt(1-Cosine.^2);
    %% Associate each solution with its nearest reference point
    [dis,pindex] = min(Distance',[],1);
    
    %% Calculate the number of associated solutions of each reference point
    rho = hist(pindex(1:N),1:NV);
    
    %% Environmental selection choose N0 from N solutions by using the localized
      % alpha-dominance to select solutions in each subregion.  
    Population = LASort(Union,pindex,N0,V,alpha);
    
    %% adjust the empty referece vector
    mark = [];emp = 0;
    for i=1:NV
       if rho(i)==0
          mark=[mark;i]; emp=emp+1;
       end
    end
    tr = rand(emp,M); tr = tr./sum(tr,2);
    V(mark,:) = tr;
     
end