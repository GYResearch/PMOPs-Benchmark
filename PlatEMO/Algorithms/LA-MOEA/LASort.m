function [Population] = LASort(Union,pindex,N0,V,alpha)
% using localized a-dominated sorting to sort out N0 solutions
% pindex indicartes the associated reference vector for each solution

% transfer solutions
[pop] = TransferFun(Union,pindex);

% localized Non-a-Dominated Sorting
[pop, F]=LADominatedSorting(pop,alpha);

% Calculate Crowding Distance
pop=CalcCrowdingDistance(pop,F);

% Sort Population
[pop, F]=SortPopulation(pop);

% Truncate
pop=pop(1:N0);
    
% Localized non-a--Dominated Sorting
[pop, F]=LADominatedSorting(pop,alpha);

% Calculate Crowding Distance
pop=CalcCrowdingDistance(pop,F);

% Sort Population
[pop, F]=SortPopulation(pop);

% Transfer pop back to population
Population = TransferBack(pop,Union);

end