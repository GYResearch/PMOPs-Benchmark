%% PMOP13 and 14 are degenrated problems.
P =1;
M = 10; %% objs
x=[0:0.001:1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=2;
B=1;
S=-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x=[0:0.01:1];
% y=[0:0.01:1];
% [u,v]=meshgrid(x,y);
% m1 = 1 + exp(sin(A.*power(u, B).*pi+pi./2))./(power(2, S).* A);
% mkk=sqrt(m1);
% f1=mkk.*u.*v;
% f2=mkk.*u.*(1-v);
% f3=mkk.*(1-u);
% plot3(f1,f2,f3,'k');
% hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% parameters in knee function;
m1 = 1 + exp(sin(A.*power(x, B).*pi+pi./2))./(power(2, S).* A);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find the index of the minima  
[~,N] = size(x); %% the size of x
[mm,mxid] = min(m1(1,1:N));
%% [mm2,mxid2] = min(m1(1,index:N));
%% corresponding x 
xknee = [];
minx = x(mxid);
%% minx2 = x(mxid2+index-1);
xknee =[xknee;minx];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nx,ny] = size(x);
vark = zeros(ny,M-1);
tmmp = repmat(xknee,ny,M-2);
vark(:,1:M-2) = tmmp;
vark(:,M-1) = x';
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kneef = 1 + exp(sin(A.*power(xknee, B).*pi+pi./2))./(power(2, S).* A);
kneef = sqrt(kneef);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PopObj = ones(N,M);  
for i=1:M
      PopObj(:,i) = PopObj(:,i).*(1.0).* kneef;
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
plot3(PopObj(:,1),PopObj(:,2),PopObj(:,3),'r*');
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
saveKnee('PMOP14-',PopObj,M);
% view(124,21)
