A1 = load('C:\Users\lurso\Desktop\PlatEMO\PlatEMO v1.5 (2017-12)\Problems\XPMOP\trueKnee\PMOP12-3.mat');
A2 = load('C:\Users\lurso\Desktop\PlatEMO\PlatEMO v1.5 (2017-12)\Problems\XPMOP\POF\XPMOP12-3.mat');
Ak1 = A1.Objs; %% true knee
Ak2 = A2.Objs; %% PoF
[n1,c1] = size(Ak1);
[n2,c2] = size(Ak2);
mam = max(Ak2(:,1:c2));
mim = min(Ak2(:,1:c2));
cm0 = [mam;mim];
pra = 2*(1/2)^(c1-2);
delta = sqrt(pdist(cm0,'euclidean'))/(4*pra*n1);
% delta = min(range)/(2*n1); %% how large the radius of the knee regions.
RefPoF = [];
for i=1:n1
    for j = 1:n2
        cm = [Ak1(i,:);Ak2(j,:)]
        dis = pdist(cm,'euclidean');
        if dis<delta
            RefPoF = [RefPoF;Ak2(j,:)];
        end
    end
end
plot3(Ak2(:,1),Ak2(:,2),Ak2(:,3),'o');
hold on
RefPoF = [RefPoF;Ak1];
RefPoF = unique(RefPoF,'rows');
saveKneeRef('RefPMOP12_3_',RefPoF,c1);
plot3(RefPoF(:,1),RefPoF(:,2),RefPoF(:,3),'*');

