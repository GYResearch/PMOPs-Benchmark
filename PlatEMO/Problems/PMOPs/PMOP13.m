classdef PMOP13 < PROBLEM
% <problem> <PMOP>
% G. Yu, Y. Jin, and M. Olhofer, ¡°Benchmark problems and performance
% indicators for search of knee points in multi-objective optimization,¡±
% IEEE Transactions on Cybernetics, 2019.
methods
    %% Initialization    
        function obj = PMOP13()
            obj.Global.A = 2; % control the number of knees, and width of knee region.
            obj.Global.B = 1; % control the bias
            obj.Global.S = -2; % control the depth of knee S = 1,/ -1
            obj.Global.p = 1; % bias in shape function     
            obj.Global.Linkage = 0; % 0/1 control the linkage function. 
            if obj.Global.M==2
                disp('it is degenerated knee problem, and M must be larger than 2.');
            end
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
             Linkage = obj.Global.Linkage;  
             % construct g function:
             g = zeros(N,1);
            if Linkage == 0   
                [g,~] = max(abs(PopDec(:,M:end)),[],2); %% g1...
            end
            
            temp = zeros(N,D-M+1);
            if  Linkage ~= 0   %% l1 linkage function in g1
                for i=1:D-M+1
                   temp(:,i) =  (1+ i./(D-M+1)).*(PopDec(:,M-1+i) - obj.Global.lower(M-1+i)) - PopDec(:,1).*(obj.Global.upper(M-1+i) - obj.Global.lower(M-1+i));   
                end 
                [g,~] = max(abs(temp(:,1:end)),[],2); %% g1...       
            end
            
           %% construct the -- degenerated -- knee functions: 
            r = zeros(N,M-2);
            k = zeros(N,1);
            for i=1:M-2
                r(:,i) = 5+10.*(PopDec(:,i)-0.5).*(PopDec(:,i)-0.5) + cos(A.*pi.*power(PopDec(:,i),B))./(A.*power(2,S));
            end
             Tr = r'; %% transposition
             if M>3
              k = prod(Tr(:,1:N))./(M-2); %% multiplicative
             end
             if M==3
              k = Tr(:,1:N)./(M-2); %% multiplicative
             end        
             k = sqrt(k'); %% transposition
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
             Linkage = obj.Global.Linkage;          
             X1 = M-1;
            lower    = zeros(1,X1);
            upper    = ones(1,X1);  
            PopDec   = rand(N,X1).*repmat(upper-lower,N,1) + repmat(lower,N,1);
            r = zeros(N,M-2);
            k = zeros(N,1);
            for i=1:M-2
                r(:,i) = 5+10.*(PopDec(:,i)-0.5).*(PopDec(:,i)-0.5) + cos(A.*pi.*power(PopDec(:,i),B))./(A.*power(2,S));
            end
            Tr = r'; %% transposition
             if M>3
              k = prod(Tr(:,1:N))./(M-2); %% multiplicative
             end
             if M==3
              k = Tr(:,1:N)./(M-2); %% multiplicative
             end        
             k = sqrt(k'); %% 

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