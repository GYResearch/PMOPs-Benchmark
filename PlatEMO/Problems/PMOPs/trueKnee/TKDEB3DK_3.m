clc;
clear all;
K=3;
M=3;
x=[0:0.001:1];
y=[0:0.001:1];
[u,v]=meshgrid(x,y);
m1 = 5+10.*(u-0.5).*(u-0.5)+2.*cos(2.*K.*pi.*u)./K;
m2 = 5+10.*(v-0.5).*(v-0.5)+2.*cos(2.*K.*pi.*v)./K;
r = (m1+m2)./2;
f1 = r.*sin(u.*pi./2).*sin(v.*pi./2);
f2 = r.*sin(u.*pi./2).*cos(v.*pi./2);
f3 = r.*cos(u.*pi./2);
plot3(f1,f2,f3);
m1 = 5+10.*(x-0.5).*(x-0.5)+2.*cos(2.*K.*pi.*x)./K;
% plot(x,m1);
hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,N] = size(x); %% the size of x
index = round(N/K); %% there are two minima of the knee function. 
[mm,mxid] = min(m1(1,1:index));
[mm2,mxid2] = min(m1(1,index:2*index));
[mm3,mxid3] = min(m1(1,2*index:N));
%%%%%%% corresponding x %%%%%%%%%%%%
xknee = [];
minx = x(mxid);
minx2 = x(mxid2+index-1);
minx3 = x(mxid3+2*index-1);
xknee =[xknee;minx;minx2;minx3];
%%%%%%%% verify the minia %%%%%%%%%%%
% for i=1:3
%     m1 = 5+10.*(xknee(i)-0.5).*(xknee(i)-0.5)+2.*cos(2.*K.*pi.*xknee(i))./K;
%     plot(xknee(i),m1,'r*');
%     hold on;
% end
%%%%%%%%%%%%%%% find the x -- knee %%
A = [1 2 3   6 7 8]
B = combnk(A,2);
X = [];
[nr,nc] = size(B);
for i=1:nr
   com = abs(B(i,2)-B(i,1));
   if ~any(com<3)
        X = [X;B(i,:)];
   end 
end
Y = mod(X,5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nr,~] = size(Y);
vark = zeros(nr,M-1);
for i=1:nr
   for j=1:M-1
    Y(i,j)
    vark(i,j) = xknee(Y(i,j),1);
   end
end
%%%%%%%%%%%%% find the knee points %%
mk = zeros(nr,M-1);
mknee = zeros(nr,1);
for i=1:M-1
  mk(:,i) = 5+10.*(vark(:,i)-0.5).*(vark(:,i)-0.5) + 2*cos(2*K*pi.*vark(:,i))./K
end 
mknee = sum(mk,2)/(M-1); %% transposition
PopObj = ones(nr,M);    
for i=1:M
    PopObj(:,i) = PopObj(:,i).* mknee;
end
 PopObj(:,1) = PopObj(:,1).*sin(vark(:,1).*pi./2).*sin(vark(:,2).*pi./2);
 PopObj(:,2) = PopObj(:,2).*sin(vark(:,1).*pi./2).*cos(vark(:,2).*pi./2);
 PopObj(:,3) = PopObj(:,3).*cos(vark(:,1).*pi./2);
 plot3(PopObj(:,1),PopObj(:,2), PopObj(:,3),'r*');
% saveKnee('DEB3DK-',PopObj,K);
view(115,16);