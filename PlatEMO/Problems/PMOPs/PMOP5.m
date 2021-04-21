classdef PMOP5 < PROBLEM
% <problem> <PMOP>
% G. Yu, Y. Jin, and M. Olhofer, ¡°Benchmark problems and performance
% indicators for search of knee points in multi-objective optimization,¡±
% IEEE Transactions on Cybernetics, 2019.
methods
    %% Initialization    
        function obj = PMOP5()
            obj.Global.A = 1; % control the number of knees, and width of knee region.
            obj.Global.B = 1; % control the bias
            obj.Global.S = 2; % control the depth of knee S = 1,
            obj.Global.p = 1; % bias in shape function
            obj.Global.l = 12; % in knee function
            obj.Global.Linkage = 0; %0/1 control the linkage function. 
            if isempty(obj.Global.M)
                obj.Global.M = 3;
            end
            if isempty(obj.Global.D)
                obj.Global.D = obj.Global.M + 9;
            end
            X1 = obj.Global.M-1;
            X2 = obj.Global.D - obj.Global.M + 1;
            obj.Global.lower    = zeros(1,obj.Global.D);
            obj.Global.upper    =  [ones(1,X1),10.*ones(1,X2)]; 
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            % construct g function:
            [N,D] =  size(PopDec);
             M      = obj.Global.M;
             A      = obj.Global.A;  
             B      = obj.Global.B;  
             S      = obj.Global.S;  
             p      = obj.Global.p;  
             l      = obj.Global.l;
             Linkage = obj.Global.Linkage;  
             % construct g function:
             g = zeros(N,1);
            if Linkage == 0    
                for i=1:D-M
                    g =g+100.*(PopDec(:,M-1+i).^2 - PopDec(:,M+i)).^2+(PopDec(:,M-1+i)-1).^2; % g5
                end
            end
            
            temp = zeros(N,D-M+1);
            if  Linkage ~= 0   %% l1 linkage function in g5
                for i=1:D-M+1
                   temp(:,i) =  (1+ i./(D-M+1)).*(PopDec(:,M-1+i) - obj.Global.lower(M-1+i)) - PopDec(:,1).*(obj.Global.upper(M-1+i) - obj.Global.lower(M-1+i));   
                end 
                for i=1:D-M
                    g =g+100.*(temp(:,i).^2 - temp(:,i+1)).^2+(temp(:,i)-1).^2; % g5
                end      
            end
             % construct the knee functions:
            r = zeros(N,M-1);
            k = zeros(N,1);
            for i=1:M-1
                r(:,i) = 2 + min(sin(2.*A.*pi.*power(PopDec(:,i),B)), cos(2.*A.*pi.*power(PopDec(:,i),B)-pi./l))./(power(2,S)); %% A.*power(2,S)
            end
            Tr = r'; %% transposition
             if M>2
              k = prod(Tr(:,1:N))./(M-1); %% multiplicative
             end
             if M==2
              k = Tr(:,1:N)./(M-1); %% multiplicative
             end
             k = (k').^(0.4); %% transposition
           %% construct the shape functions:           
            PopObj = ones(N,M);    
            for i=1:M
                 PopObj(:,i) = PopObj(:,i).*(1.0 + g).* k;
            end
            for i=1:M
               for j=1:M-i
                  PopObj(:,i) =  PopObj(:,i).*power(PopDec(:,j), p);
               end
               if i~=1
                   aux = M - i +1;
                   PopObj(:,i)= PopObj(:,i).*(1 - power(PopDec(:,aux), p));
               end
            end
        end
        function P = PF(obj,N) 
             M      = obj.Global.M;
             A      = obj.Global.A;  
             B      = obj.Global.B;  
             S      = obj.Global.S;  
             p      = obj.Global.p;  
             l      = obj.Global.l;
             Linkage = obj.Global.Linkage;          
             X1 = M-1;
            lower    = zeros(1,X1);
            upper    = ones(1,X1);  
            PopDec   = rand(N,X1).*repmat(upper-lower,N,1) + repmat(lower,N,1);
            r = zeros(N,M-1);
            k = zeros(N,1);
            for i=1:M-1
                r(:,i) = 2 + min(sin(2.*A.*pi.*power(PopDec(:,i),B)), cos(2.*A.*pi.*power(PopDec(:,i),B)-pi./l))./(power(2,S)); %% A.*power(2,S)
            end
           Tr = r'; %% transposition
             if M>2
              k = prod(Tr(:,1:N))./(M-1); %% multiplicative
             end
             if M==2
              k = Tr(:,1:N)./(M-1); %% multiplicative
             end
              k = (k').^(0.4); %% transposition

            PopObj = ones(N,M);    
            for i=1:M
                  PopObj(:,i) = PopObj(:,i).* k;
            end
            for i=1:M
                for j=1:M-i
                    PopObj(:,i) =  PopObj(:,i).*power(PopDec(:,j), p);
                end
                if i~=1
                 aux = M - i +1;
                  PopObj(:,i)= PopObj(:,i).*(1 - power(PopDec(:,aux), p));
                end
            end
            [FrontNo,MaxFNo] = NDSort(PopObj,N);
            
            P = PopObj(FrontNo==1,:);
        end
    end
end