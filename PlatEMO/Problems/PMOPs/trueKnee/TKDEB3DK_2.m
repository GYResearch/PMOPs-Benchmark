clc;
clear all;
K=2;
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
hold on;
m1 = 5+10.*(x-0.5).*(x-0.5)+2.*cos(2.*K.*pi.*x)./K;
% plot(x,m1);
[~,N] = size(x); %% the size of x
index = round(N/2); %% there are two minima of the knee function. 
[mm,mxid] = min(m1(1,1:index));
[mm2,mxid2] = min(m1(1,index:N));
%% corresponding x 
xknee = [];
minx = x(mxid);
minx2 = x(mxid2+index-1);
xknee =[xknee;minx;minx2];
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
  mk(:,i) = 5+10.*(vark(:,i)-0.5).*(vark(:,i)-0.5) + 2*cos(2*K*pi.*vark(:,i))./K
end 
mknee = sum(mk,2)/(M-1); %% transposition
PopObj = ones(nrf,M);    
for i=1:M
    PopObj(:,i) = PopObj(:,i).* mknee;
end
 PopObj(:,1) = PopObj(:,1).*sin(vark(:,1).*pi./2).*sin(vark(:,2).*pi./2);
 PopObj(:,2) = PopObj(:,2).*sin(vark(:,1).*pi./2).*cos(vark(:,2).*pi./2);
 PopObj(:,3) = PopObj(:,3).*cos(vark(:,1).*pi./2);
 plot3(PopObj(:,1),PopObj(:,2), PopObj(:,3),'r*');
% saveKnee('DEB3DK-',PopObj,K);
view(115,16);