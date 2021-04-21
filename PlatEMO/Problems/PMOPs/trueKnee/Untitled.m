clc;
clear all;
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
