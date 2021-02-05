function [f,c1,c2,c3] = IHR(TF,tic,sampling_rate,resolution,pad)
%% curve extraction for fundamental frequency and multiples
% TF  = the time frequency anlysis for the epoch
%%
    tfrsqtic = tic;
    basicTF.fs = sampling_rate;
    basicTF.fr = resolution;
    allConceFT = TF;
    %%
    idx0 = find(tfrsqtic*basicTF.fs>0.75&tfrsqtic*basicTF.fs<1.5) ;
    [f] = CurveExt_M(abs(allConceFT(idx0,:))',1) ;
    f = f + idx0(1) - 1 ;
    rg = quantile(f(5*pad+1:5*(pad+4)),0.75)/2*resolution;
    rg = round(rg,1);
    ConceFT_c1 = zeros(1+2*rg/resolution,length(f));
    for i = 1:length(f)
        ConceFT_c1(:,i) = allConceFT(2*f(i)-rg/basicTF.fr:2*f(i)+rg/basicTF.fr,i);
    end
    idx1 = find(tfrsqtic*basicTF.fs<2*rg) ;
    [c] = CurveExt_M(abs(ConceFT_c1(:,:))',1) ;
    c = c + idx1(1) - 1 ;
    c1 = c + 2*f - rg/basicTF.fr; 


    %%
    f_median = round(median(f(5*pad+1:5*(pad+4))));
    ConceFT_c2 = zeros(1+2*rg/basicTF.fr,length(f));
    for i = 1:length(f)
        ConceFT_c2(:,i) = allConceFT(3*f_median-rg/basicTF.fr:3*f_median+rg/basicTF.fr,i);
    end
    
    idx2 = find(tfrsqtic*basicTF.fs<2*rg) ;
    [c] = CurveExt_M(abs(ConceFT_c2(:,:))',1) ;
    c = c + idx2(1) - 1 ;
    c2 = c + 3*f_median - rg/basicTF.fr; 

    %%
    ConceFT_c3 = zeros(1+2*rg/basicTF.fr,length(f));
    for i = 1:length(f)
        ConceFT_c3(:,i) = allConceFT(4*f_median-rg/basicTF.fr:4*f_median+rg/basicTF.fr ,i);
    end
    idx3 = find(tfrsqtic*basicTF.fs<2*rg) ;
    [c] = CurveExt_M(abs(ConceFT_c3(:,:))',1) ;
    c = c + idx3(1) - 1 ;
    c3 = c + 4*f_median - rg/basicTF.fr; 