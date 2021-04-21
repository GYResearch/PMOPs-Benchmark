function diff=LADominates(x,y,alpha)
    
    xobjs = x.Cost;
    yobjs = y.Cost;
    
    [~,M] = size(xobjs);
    g = zeros(1,M);
    
    for i=1:M
       g(i) =  xobjs(1,i)-yobjs(1,i);
       for j=1:M
           if j~=i
               g(i) = g(i)+ alpha*(xobjs(1,j)-yobjs(1,j));
           end
       end
    end
    
    if x.Associate==y.Associate    
        diff=all(g<=0) && any(g<0);
    end
    if x.Associate~=y.Associate  
        diff = false;
    end

end