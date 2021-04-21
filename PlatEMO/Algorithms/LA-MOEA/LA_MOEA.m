function LA_MOEA(Global)
% <algorithm> <L>
% Localized alpha-MOEA:
% Yu G, Jin Y, Olhofer M. An a priori knee identification multi-objective evolutionary algorithm based on 
% ¦Á-dominance[C]//Proceedings of the Genetic and Evolutionary Computation Conference Companion. 2019: 241-242.

    %% Generate random population
    Population = Global.Initialization();
    % the number of reference vectors may influence the performance of the
    % algorithm, suggest to use small number of reference vectors.
    [V0,~] = UniformPoint(Global.N/10,Global.M);
    V             = V0;
    [Vs,~] = size(V);
    fprintf('Number of reference vectors: %d\n', Vs);
    % when alpha is close to 1.0, then more converged solution set will be
    % obtained. Suggest to use 0.75 or 1.0
    alpha = 0.75;
    index = 0;
    %% Optimization
    while Global.NotTermination(Population)
        MatingPool = randi(length(Population),1,Global.N);
        Offspring  = GA(Population(MatingPool));
        [Population,V] = LAEnvironmentalSelection([Population,Offspring],V,Global.N,alpha);
        index = Global.evaluated/Global.N;
        fprintf('-----The %d th generation------\n',index);
        [FN,~] = NDSort(Population.objs,Population.cons,Global.N);
        Next = FN ==1;
        Pop = Population(Next);
        address1 = ['Pop',num2str(index),'.mat'];
        save(address1,'Pop');
        address2 = ['Reference',num2str(index),'.mat'];
        save(address2,'V');
        fprintf('-----Save Pop and V successfully------\n');
    end
end