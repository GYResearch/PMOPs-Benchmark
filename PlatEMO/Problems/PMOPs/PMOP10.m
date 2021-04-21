classdef PMOP10 < PROBLEM
% <problem> <PMOP>
% G. Yu, Y. Jin, and M. Olhofer, ¡°Benchmark problems and performance
% indicators for search of knee points in multi-objective optimization,¡±
% IEEE Transactions on Cybernetics, 2019.
methods
    %% Initialization    
        function obj = PMOP10()
            obj.Global.A = 1; % control the number of knees, and width of knee region.
            obj.Global.B = 1; % control the bias
            obj.Global.S = 2; % control the depth of knee S = 1,
            obj.Global.p = 1; % bias in shape function   
            obj.Global.l = 12; %% 12
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
             g3 = zeros(N,1);
            g7 = zeros(N,1);
            if Linkage == 0   
                %% calc g3
                g3 = PopDec(:,M:end).^2 - 10.*cos(4.*pi.*PopDec(:,M:end)); %% g3...
                g3 = 1+10.*(D-M+1)+sum(g3,2);
                %% calc g7 
                len = D-M+1;
                temp1 = zeros(N,D);
                temp2 = ones(N,1);
                for i=1:len
                    temp1 = temp1 + PopDec(:,i+M-1).^2./4000;
                    temp2 = temp2.*cos(PopDec(:,i+M-1)./sqrt(i)); 
                end
                     g = sum(temp1,2) + temp2 + 1; % g7
             end
         
            temp = zeros(N,D-M+1);
            if  Linkage ~= 0   %% l1 linkage function in g7
                temp = zeros(N,D-M+1);
                for i=1:D-M+1
                   temp(:,i) =  (1+ cos(0.5.*pi.*i./(D-M+1))).*(PopDec(:,M-1+i) - obj.Global.lower(M-1+i)) - PopDec(:,1).*(obj.Global.upper(M-1+i) - obj.Global.lower(M-1+i));   
                end
               %% calc g3
                g3 =  temp(:,1:end).^2 - 10.*cos(4.*pi.* temp(:,1:end)); %% g3...
                g3 = 1+10.*(D-M+1)+sum(g3,2);
               %% calc g7 
                len = D-M+1;
                temp1 = zeros(N,D);
                temp2 = ones(N,1);
                for i=1:len
                   temp1 = temp1 + temp(:,i).^2./4000;
                   temp2 = temp2.*cos(temp(:,i)./sqrt(i)); 
                end
                  g = sum(temp1,2) + temp2 + 1; % g7
            end
             % construct the knee functions:
            r = zeros(N,M-1);
            k = zeros(N,1);
            for i=1:M-1
                r(:,i) = 2 + min(sin(2.*A.*pi.*power(PopDec(:,i),B)), cos(2.*A.*pi.*power(PopDec(:,i),B)-pi./l))./(power(2,S)); %%A.*power(2,S)
            end
            Tr = r'; %% transposition
             if M>2
              k = prod(Tr(:,1:N))./(M-1); %% multiplicative
             end
             if M==2
              k = Tr(:,1:N)./(M-1); %% multiplicative
             end
             k = (k').^(0.2); %% transposition
           %% construct the shape functions:           
            PopObj = ones(N,M);    
            for i=1:M
                if mod(i,2) == 1
                 PopObj(:,i) = PopObj(:,i).*(1.0 + g3).* k;
                end
                if mod(i,2) == 0 
                 PopObj(:,i) = PopObj(:,i).*(1.0 + g7).* k;
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
             l      = obj.Global.l;
             Linkage = obj.Global.Linkage;          
             X1 = M-1;
            lower    = zeros(1,X1);
            upper    = ones(1,X1);  
            PopDec   = rand(N,X1).*repmat(upper-lower,N,1) + repmat(lower,N,1);
            r = zeros(N,M-1);
            k = zeros(N,1);
            for i=1:M-1
                r(:,i) = 2 + min(sin(2.*A.*pi.*power(PopDec(:,i),B)), cos(2.*A.*pi.*power(PopDec(:,i),B)-pi./l))./(power(2,S)); %%A.*power(2,S)
            end
             Tr = r'; %% transposition
             if M>2
              k = prod(Tr(:,1:N))./(M-1); %% multiplicative
             end
             if M==2
              k = Tr(:,1:N)./(M-1); %% multiplicative
             end
             k = (k').^(0.2); %% transposition
            

            PopObj = ones(N,M);    
            for i=1:M
                if mod(i,2) == 1
                 PopObj(:,i) = PopObj(:,i).*(1.0 + 1.0).* k; %% min(g3) = 1
                end
                if mod(i,2) == 0 
                 PopObj(:,i) = PopObj(:,i).*(1.0 + 0).* k; %% min(g7) = 0
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
            [FrontNo,MaxFNo] = NDSort(PopObj,N);
           
            P = PopObj(FrontNo==1,:);
        end
    end
end