A1 = load('C:\Users\lurso\Desktop\PlatEMO\PlatEMO v1.5 (2017-12)\Problems\XPMOP\trueKnee\PMOP12-8.mat');
A2 = load('C:\Users\lurso\Desktop\PlatEMO\PlatEMO v1.5 (2017-12)\Problems\XPMOP\POF\XPMOP12-8.mat');
Ak1 = A1.Objs; %% true knee
Ak2 = A2.Objs; %% PoF
[n1,c1] = size(Ak1);
[n2,c2] = size(Ak2);
mam = max(Ak2(:,1:c2));
mim = min(Ak2(:,1:c2));
% range = mam - mim;
cm0 = [mam;mim];
pra = 1.5*(1/2)^(c1-2);
delta = sqrt(pdist(cm0,'euclidean'))/(pra*n1);
% delta = min(range)/(n1);
RefPoF = [];
for i=1:n1
    for j = 1:n2
        cm = [Ak1(i,:);Ak2(j,:)]
        dis = sqrt(pdist(cm,'euclidean'));
        if dis<delta
            RefPoF = [RefPoF;Ak2(j,:)];
        end
    end
end
% plot(Ak2(:,1),Ak2(:,2),'o');
hold on
RefPoF = [RefPoF;Ak1];
RefPoF = unique(RefPoF,'rows');
saveKneeRef('RefPMOP12_8_',RefPoF,c1);
% plot(RefPoF(:,1),RefPoF(:,2),'o');

