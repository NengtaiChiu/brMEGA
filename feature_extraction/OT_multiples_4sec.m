function [OT,allOT1,power] = OT_multiples_4sec(TF,f,c1,c2,c3,pad,targetlength)
    allConceFT = TF;
    m = length(TF(:,1));
    n = length(f);
    l_i = zeros(m,n);
    v_i = zeros(m,n);
    
    k = find(~allConceFT);
    allConceFT(k) = 1.0000e-05;
   for i = 1:length(f)
        allConceFT(1:f(i)-10,i)=1.0000e-05;
        allConceFT(c3(i)+15:end,i) = 1.0000e-05;
    end
    %%
    for yy = 1:n
        
        l_i(f(yy),yy) = sum(allConceFT(f(yy)-2:f(yy)+2,yy));
        l_i(c1(yy),yy) = sum(allConceFT(c1(yy)-2:c1(yy)+2,yy));%keyboard
        l_i(c2(yy),yy) = sum(allConceFT(c2(yy)-2:c2(yy)+2,yy));
        l_i(c3(yy),yy) = sum(allConceFT(c3(yy)-2:c3(yy)+2,yy));
        l_i(:,yy) =  l_i(:,yy)./sum(l_i(:,yy));
        
        v_i(f(yy)-15:f(yy)+15,yy)= abs(allConceFT(f(yy)-15:f(yy)+15,yy));
        v_i(c1(yy)-15:c1(yy)+15,yy)= abs(allConceFT(c1(yy)-15:c1(yy)+15,yy));
        v_i(c2(yy)-15:c2(yy)+15,yy)= abs(allConceFT(c2(yy)-15:c2(yy)+15,yy));
        v_i(c3(yy)-15:c3(yy)+15,yy)= abs(allConceFT(c3(yy)-15:c3(yy)+15,yy));
        v_i(:,yy) = v_i(:,yy)./sum(v_i(:,yy));% power spectrum of the real TF result, normalized 
    end
    %%
     allOT1 = zeros(1,length(f));
    for i = 1:length(f)
        allOT1(i) = slicedOT(v_i(:,i),l_i(:,i));
    end
    OT = median(allOT1(pad*5+1:(pad+targetlength)*5));
    power = sum(sum(abs(allConceFT(:,pad*5+1:(pad+targetlength)*5))));