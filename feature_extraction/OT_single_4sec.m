function [OT,allOT1] = OT_single_4sec(TF,f,pad,targetlength)
    allConceFT = TF;
    m = length(TF(:,1));
    n = length(f);
    l_i = zeros(m,n);
    v_i = zeros(m,n);
    
    k = find(~allConceFT);
    allConceFT(k) = 1.0000e-05;
   for i = 1:length(f)
        allConceFT(1:f(i)-10,i)=1.0000e-05;
        allConceFT(f(i)+15:end,i) = 1.0000e-05;
    end
    %%
    for yy = 1:n
        
        l_i(f(yy),yy) = sum(allConceFT(f(yy)-2:f(yy)+2,yy));
        l_i(:,yy) =  l_i(:,yy)./sum(l_i(:,yy));
        
        v_i(f(yy)-15:f(yy)+15,yy)= abs(allConceFT(f(yy)-15:f(yy)+15,yy));
        v_i(:,yy) = v_i(:,yy)./sum(v_i(:,yy));% power spectrum of the real TF result, normalized 
    end
    %%
     allOT1 = zeros(1,length(f));
    for i = 1:length(f)
        allOT1(i) = slicedOT(v_i(:,i),l_i(:,i));
    end
    OT = median(allOT1(pad*5+1:(pad+targetlength)*5));