clc;
clear all;
K=3;
S=1;
x=0:0.001:1;
m1 = 5+10.*(x-0.5).*(x-0.5)+cos(2.*K.*pi.*x).*power(2,S./2)./K;
f1 = m1.*(sin(x.*pi./power(2,S+1)+(1+(power(2,S)-1)./(power(2,S+2))).*pi)+1);
f2 = m1.*(cos(pi+x.*pi./2)+1);
plot(f1,f2);
% % % % plot(x,m1);
hold on;
%% find the index of the minima  
[~,N] = size(x); %% the size of x
index = round(N/K); %% there are two minima of the knee function. 
[mm,mxid] = min(m1(1,1:index));
[mm2,mxid2] = min(m1(1,index:2*index));
[mm3,mxid3] = min(m1(1,2*index:N));
xknee = [];
minx = x(mxid);
minx2 = x(mxid2+index-1);
minx3 = x(mxid3+2*index-1);
xknee =[xknee;minx;minx2;minx3];
%% verify the minima
%% save the true knees
PopObj = zeros(3,2);
for i=1:3
    m1 = 5+10.*(xknee(i,1)-0.5).*(xknee(i,1)-0.5)+cos(2.*K.*pi.*xknee(i,1)).*power(2,S./2)./K;
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


