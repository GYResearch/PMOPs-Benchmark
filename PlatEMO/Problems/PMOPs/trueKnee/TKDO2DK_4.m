clc;
clear all;
K=4;
S=1;
x=0:0.001:1;
m1 = 5+10.*(x-0.5).*(x-0.5)+cos(2.*K.*pi.*x).*power(2,S./2)./K;
f1 = m1.*(sin(x.*pi./power(2,S+1)+(1+(power(2,S)-1)./(power(2,S+2))).*pi)+1);
f2 = m1.*(cos(pi+x.*pi./2)+1);
plot(f1,f2);
hold on;
%% find the index of the minima  
[~,N] = size(x); %% the size of x
index = round(N/K); %% there are two minima of the knee function. 
[mm,mxid] = min(m1(1,1:index));
[mm2,mxid2] = min(m1(1,index:2*index));
[mm3,mxid3] = min(m1(1,2*index:3*index));
[mm4,mxid4] = min(m1(1,3*index:N));
xknee = [];
minx = x(mxid);
minx2 = x(mxid2+index-1);
minx3 = x(mxid3+2*index-1);
minx4 = x(mxid4+3*index-1);
xknee =[xknee;minx;minx2;minx3;minx4];
PopObj = zeros(4,2);
for i=1:4
    m1 = 5+10.*(xknee(i,1)-0.5).*(xknee(i,1)-0.5)+cos(2.*K.*pi.*xknee(i,1)).*power(2,S./2)./K
    f1 = m1.*(sin(xknee(i,1).*pi./power(2,S+1)+(1+(power(2,S)-1)./(power(2,S+2))).*pi)+1);
    f2 = m1.*(cos(pi+xknee(i,1).*pi./2)+1);
    PopObj(i,1) = f1;
    PopObj(i,2) = f2;
% %     plot(f1,f2,'*');
% %     plot(xknee(i,1),m1,'*');
    hold on;
end
plot(PopObj(:,1),PopObj(:,2),'r*');
saveKnee('DO2DK-',PopObj,K);



