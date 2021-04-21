classdef PMOP14 < PROBLEM
% <problem> <PMOP>
% G. Yu, Y. Jin, and M. Olhofer, ¡°Benchmark problems and performance
% indicators for search of knee points in multi-objective optimization,¡±
% IEEE Transactions on Cybernetics, 2019.
methods
    %% Initialization    
        function obj = PMOP14()
            obj.Global.A = 2; % control the number of knees, and width of knee region.
            obj.Global.B = 1; % control the bias
            obj.Global.S = -1; % control the depth of knee S = 1,/ -1
            obj.Global.p = 1; % bias in shape function     
            obj.Global.Linkage = 0; %0/1 control the linkage function. 
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
             g6 = zeros(N,1);
            g8 = zeros(N,1);
            if Linkage == 0  
                %% calc g6
                temp1 = zeros(N,D-M+1);
                temp1 = PopDec(:,M:end).^2 - 10.*cos(2.*pi.*PopDec(:,M:end)) + 10;
                g6 = sum(temp1,2); % g6
                %% calc g8
                temp2 = zeros(N,D-M+1);
                temp3 = zeros(N,D-M+1); 
                temp2 = sum(PopDec(:,M:end).^2,2)./(D-M+1);
                temp3 = sum(cos(2.*pi.*PopDec(:,M:end)),2)./(D-M+1);
                g8 = -20.*exp(-0.2.*sqrt(temp2)) - exp(temp3) + 20 + exp(1); %% g8... 
            end
           if  Linkage ~= 0   %% l1 linkage function in g6 g8
                temp = zeros(N,D-M+1);
                for i=1:D-M+1
                   temp(:,i) =  (1+ i./(D-M+1)).*(PopDec(:,M-1+i) - obj.Global.lower(M-1+i)) - PopDec(:,1).*(obj.Global.upper(M-1+i) - obj.Global.lower(M-1+i));   
                end 
                %% calc g6
                temp1 = zeros(N,D-M+1);
                temp1 = temp(:,1:end).^2 - 10.*cos(2.*pi.*temp(:,1:end)) + 10;
                g6 = sum(temp1,2); % g6
                %% calc g8
                temp2 = zeros(N,D-M+1);
                temp3 = zeros(N,D-M+1); 
                temp2 = sum(temp(:,1:end).^2,2)./(D-M+1);
                temp3 = sum(cos(2.*pi.*temp(:,1:end)),2)./(D-M+1);
                g8 = -20.*exp(-0.2.*sqrt(temp2)) - exp(temp3) + 20 + exp(1); %% g8...
            end
             % construct the knee functions:
            r = zeros(N,M-2);
            k = zeros(N,1);
            for i=1:M-2
                r(:,i) = 1 + exp(sin(A.*power(PopDec(:,i), B).*pi+pi./2))./(power(2, S).* A);
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
                if mod(i,2) == 1
                 PopObj(:,i) = PopObj(:,i).*(1.0 + g6).* k;
                end
                if mod(i,2) == 0 
                 PopObj(:,i) = PopObj(:,i).*(1.0 + g8).* k;
                end
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
                r(:,i) = 1 + exp(sin(A.*power(PopDec(:,i), B).*pi+pi./2))./(power(2, S).* A);
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
                 PopObj(:,i) = PopObj(:,i).*(1.0).* k;
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