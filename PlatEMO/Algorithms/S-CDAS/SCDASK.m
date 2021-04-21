function SCDASK(Global)
% <algorithm> <H-N>
% Self-controlling dominance area of solutions
% Runs --- 20  --- 20 independent runs
% Max_Evas ---1000 --- 1000 Inner_Evas

%------------------------------- Reference --------------------------------
% H. Sato, H. E. Aguirre, and K. Tanaka, Self-controlling dominance area of
% solutions in evolutionary many-objective optimization, Proceedings of the
% Asia-Pacific Conference on Simulated Evolution and Learning, 2010,
% 455-465.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    [Runs,Max_Evas] = Global.ParameterSet(20,10000);
    %% Generate random population
    Population = Global.Initialization();
  
    %% Optimization
    while Global.NotTermination(Population)
        for run = 1:Runs
           fprintf('Runs = %d\n',run);
           Population = Global.Initialization(); 
           [~,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Global.N);
           Eva_ = Global.N;
           while Eva_ < Max_Evas
                MatingPool = TournamentSelection(2,Global.N,FrontNo,-CrowdDis);
                Offspring  = GA(Population(MatingPool));
                [noff,~] = size(Offspring.objs);
                Eva_ = Eva_ + noff;
                fprintf('Evalutions: %d\n',Eva_);
                [Population,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],Global.N);
                if Eva_ >= Max_Evas
                    address = ['C:\Users\Lursonkj\Desktop\KD-MOEA.tex - 4.0\PlatEMO-master\PlatEMO\Data\SCDAS\',num2str(Global.M),'\',num2str(run),'.mat'];
                    save(address, 'Population');
                    fprintf('SAVE Knees \n');
                end
           end
        end
       
    end
end