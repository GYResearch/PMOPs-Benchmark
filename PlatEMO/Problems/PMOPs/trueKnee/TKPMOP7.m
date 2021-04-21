P =1;
M = 10; %% objs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=[0:0.001:1];
A=4;
B=1;
S=2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x=[0:0.001:1];
% y=[0:0.001:1];
% [u,v]=meshgrid(x,y);
% m1 = 1 + exp(cos(A.*power(u, B).*pi+pi./2))./(power(2, S).* A);
% m2 = 1 + exp(cos(A.*power(v, B).*pi+pi./2))./(power(2, S).* A);
% ma=(m1.*m2)./2;
% m=3.^(ma)
% f1=m.*u.*v;
% f2=m.*u.*(1-v);
% f3=m.*(1-u);
% plot3(f1,f2,f3,'k');
% hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% parameters in knee function;
m1 =  1 + exp(cos(A.*power(x, B).*pi+pi./2))./(power(2, S).* A);
%% find the index of the minima  
[~,N] = size(x); %% the size of x
index = round(N/2); %% there are two minima of the knee function. 
[mm,mxid] = min(m1(1,1:index));
[mm2,mxid2] = min(m1(1,index:N));
%% corresponding x 
xknee = [];
minx = x(mxid);
minx2 = x(mxid2+index-1);
xknee =[xknee;minx;minx2];
fxknee = repmat(xknee',M-1,1); %% the xknee set.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% parameters in knee function;
m1 = 1 + exp(cos(A.*power(x, B).*pi+pi./2))./(power(2, S).* A);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find the index of the minima  
[~,N] = size(x); %% the size of x
index = round(N/2); %% there are two minima of the knee function. 
[mm,mxid] = min(m1(1,1:index));
[mm2,mxid2] = min(m1(1,index:N));
%% corresponding x 
xknee = [];
minx = x(mxid);
minx2 = x(mxid2+index-1);
xknee =[xknee;minx;minx2];
fxknee = repmat(xknee',M-1,1); %% the xknee set.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% do combinations among the fxknee. [1 2;4 5;7 8;10 11;13 14....] 
%% if there are 1 interval, then there are in the same line
com = zeros(M-1,2); %% 2 --> two minima
for i=1:M-1
   com(i,1) = 3*i-2;
   com(i,2) = 3*i-1;
end
com2 = reshape(com,[1,2*(M-1)]);
com3 = combnk(com2,M-1);
fcom = [];
[nr,nc] = size(com3);
for i=1:nr
    com4 = combnk(com3(i,:),2);
    [row,col] = size(com4);
    com5 = com4(:,2)-com4(:,1);
    if ~any((com5+1)==0) & ~any((com5-1)==0)
        fcom = [fcom;com3(i,:)];
    end
end
fcom = sort(fcom,2); %% sort array and save the position of the xknee.
[nrf,ncf] = size(fcom);
oxknee = xknee';
vark = []; %% save the true x knee point
ffcom = mod(fcom,3);
for i=1:nrf
    for j =1: ncf
        if ffcom(i,j) == 1
            tmp =  oxknee(1,1);
            vark(i,j) = tmp;
        end
        if  ffcom(i,j) == 2
            tmp =  oxknee(1,2);
            vark(i,j) = tmp;
        end
    end
end
mk = zeros(nrf,M-1);
mknee = zeros(nrf,1);
for i=1:M-1
    mk(:,i)  =  1 + exp(cos(A.*power(vark(:,i), B).*pi+pi./2))./(power(2, S).* A);
end
Tr = mk'; %% transposition
if M>2
   mknee = prod(Tr(:,1:nrf))./(M-1); %% multiplicative
end
if M==2
   mknee = Tr(:,1:nrf)./(M-1); %% multiplicative
end
mknee = 3.^(mknee'); %% transposition
%% construct the shape functions:           
PopObj = ones(nrf,M);     
for i=1:M
    PopObj(:,i) = PopObj(:,i).*mknee;  
 end
 for i=1:M
               for j=1:M-i
                  PopObj(:,i) =  PopObj(:,i).*power(vark(:,j), P);
               end
               if i~=1
                   aux = M - i +1;
                   PopObj(:,i)= PopObj(:,i).*(1 - power(vark(:,aux), P));
               end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot3(PopObj(:,1),PopObj(:,2),PopObj(:,3),'r*');
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
saveKnee('PMOP7-',PopObj,M);
% view(124,21)
