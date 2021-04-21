clc;
clear all;
K=5;
x=0:0.001:1;
x = [0:0.001:1];
m1 = 5+10.*(x-0.5).*(x-0.5)+cos(2.*K.*pi.*x)./K;
f1 = m1.*sin(x.*pi./2);
f2 = m1.*cos(x.*pi./2);
plot(f1,f2);
% plot(x,m1);
hold on;
%% find the index of the minima  
[~,N] = size(x); %% the size of x
index = round(N/K); %% there are two minima of the knee function. 
[mm,mxid] = min(m1(1,1:0.6*index));
[mm2,mxid2] = min(m1(1,index:2*index));
[mm3,mxid3] = min(m1(1,2*index:3*index));
[mm4,mxid4] = min(m1(1,3.409*index:4*index));
[mm5,mxid5] = min(m1(1,4.41*index:N));

xknee = [];
minx = x(mxid);
minx2 = x(mxid2+index-1);
minx3 = x(mxid3+2*index-1);
minx4 = x(round(3.409*index));
minx5 = x(round(4.41*index));
xknee =[xknee;minx;minx2;minx3;minx4;minx5];
%% verify the minima
PopObj = zeros(5,2);
for i=1:5
    m1 = 5+10*(xknee(i,1)-0.5).*(xknee(i,1)-0.5)+cos(2*K.*pi.*xknee(i,1))./K;
    f1 = m1.*sin(xknee(i,1).*pi./2);
    f2 = m1.*cos(xknee(i,1).*pi./2);
    PopObj(i,1) = f1;
    PopObj(i,2) = f2;
end
plot(PopObj(:,1),PopObj(:,2),'r*');
saveKnee('DEB2DK-',PopObj,K);
