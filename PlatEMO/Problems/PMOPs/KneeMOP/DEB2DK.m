classdef DEB2DK < PROBLEM
% <problem> <Knee>
% Comparison of Multiobjective Evolutionary Algorithms for knee identification
methods
 
    %% Initialization    
        function obj = DEB2DK()
            obj.Global.Knee = 4; %%4/5 control the number of knee regions
            obj.Global.M = 2;
            if isempty(obj.Global.D)
                obj.Global.D = 7; %% 30
            end
            obj.Global.lower    = zeros(1,obj.Global.D);
            obj.Global.upper    = ones(1,obj.Global.D);
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            K = obj.Global.Knee;
            g = 1+9.*sum(PopDec(:,2:end),2)./(obj.Global.D -1);
            r = 5+10.*(PopDec(:,1)-0.5).*(PopDec(:,1)-0.5)+cos(2.*K.*pi.*PopDec(:,1))./K;
            PopObj(:,1) = g.*r.*sin(PopDec(:,1).*pi./2);
            PopObj(:,2) = g.*r.*cos(PopDec(:,1).*pi./2);
        end
        function P = PF(obj,N) 
            M = obj.Global.M;          
            X1 = M-1;
            K = obj.Global.Knee;
            lower    = zeros(1,X1);
            upper    = ones(1,X1);  
            PopDec   = rand(N,X1).*repmat(upper-lower,N,1) + repmat(lower,N,1);
            r = zeros(N,M-1);
            k = zeros(N,1);
            for i=1:M-1     
                r(:,i) = 5+10.*(PopDec(:,i)-0.5).*(PopDec(:,i)-0.5)+cos(2.*K.*pi.*PopDec(:,i))./K;
            end
            k = sum(r,2)./(M-1); %
            f = [];
            f(:,1) = k.*sin(PopDec(:,1).*pi./2);
            f(:,2) = k.*cos(PopDec(:,1).*pi./2);
            [FrontNo,MaxFNo] = NDSort(f,N);
            t=1;
            for i=1:N   
                if FrontNo(i) ==1
                  P(t,:) = f(i,:);
                  t = t+1;
                end
            end
        end
    end
end