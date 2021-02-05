function [all_pk] = NLEM_new(I, Loc, N, fs, BFORDER, CUTOFF, SHRINKAGE, DIFFUSION)

xd = zeros(size(I));


L_1 = 0.2*fs;
L_2 = 0.5*fs;
BT = 5 ;
if Loc(1) < L_1 ; Loc = Loc(2:end) ; end
if Loc(end)+ L_2 > length(length(I)) ; Loc = Loc(1:end-1) ; end


[b_hp,a_hp] = butter(BFORDER, CUTOFF/(fs/2),'high');
IHP = filtfilt(b_hp,a_hp, I);
% window = 0.05;%seconds
%     K = floor(window/2*fs);
%     [trend]=median_filter(IHP,K);
%     trend = smooth(trend,3*K,'loess');%smooth the trend
%     IHPks = IHP -trend;
IHPks = smooth(IHP, 'loess', fs*0.005) ;  % ks = kernel smoothing


% You can well replace this high pass filter step by subtracting the original signal
% by the ICA reconstructed EEG signal. The goal is to get a "clean underlying
% artifact" to determine the neighbors so that the EEG signal is better preserved
% N = numbers of KNN to search

Nbeat = length(Loc) ;
fprintf(['\tOnly ',num2str(Nbeat),' beats\n']) ;
beats = zeros(L_1+L_2+1, Nbeat) ;

beatsHP = zeros(L_1+L_2+1, Nbeat) ;
beatsHPks = zeros(L_1+L_2+1, Nbeat) ;

beatsidx = zeros(L_1+L_2+1, Nbeat) ;

beatHeight = zeros(1, Nbeat) ;
beatWidth = zeros(1, Nbeat) ;

    % to evaluate the noise level
rawEEG = nan(size(I)) ;

fprintf('\tPreparing dataset...\n') ;
for ii = 1:Nbeat
    idx = Loc(ii)-L_1: Loc(ii) + L_2  ;
    beats(:,ii) = I(idx) ;
    beatsHP(:,ii) = IHP(idx) ;
    beatsHPks(:,ii) = IHPks(idx) ;
    beatsidx(:,ii) = idx ;
    beatHeight(:, ii) = I(Loc(ii)) ; 
    if ii > 1
        RRI = Loc(ii) - Loc(ii-1);
        if RRI < 1.5 *fs
        rawEEG(Loc(ii-1)+L_2+1:Loc(ii)-L_1-1) = I(Loc(ii-1)+L_2+1:Loc(ii)-L_1-1) ; 
        end
    end
end


%%
fprintf('\tFinding neighbors...\n') ; t0=tic;
clear X ;
%%
ttt = 500;
Idnx_1 = randperm(Nbeat);
group_len = floor(Nbeat/ttt);
%%
all_pk = zeros(size(beatsidx));
X_ops = zeros(size(beatsidx));
for rr = 1:group_len
    if rr ~= group_len
        nn = ttt;
        beats_seg = zeros(L_1+L_2+1, nn) ;
        beatsidx_seg =  zeros(L_1+L_2+1, nn) ;
        for pp = 1:ttt
            tt = (rr-1)*nn+pp;
            idx_s = Loc(Idnx_1(tt))-L_1: Loc(Idnx_1(tt)) + L_2  ;
            beatsidx_seg(:,pp) = idx_s ;
            beats_seg(:,pp) = IHPks(idx_s);
        end
        
    else
        nn = length(Idnx_1) - (group_len-1)*ttt ;
        beats_seg = zeros(L_1+L_2+1, nn) ;
        beatsidx_seg =  zeros(L_1+L_2+1, nn) ;
        
        for pp = 1:nn
            tt = (rr-1)*ttt+pp;
            idx_s = Loc(Idnx_1(tt))-L_1: Loc(Idnx_1(tt)) + L_2  ;
            beatsidx_seg(:,pp) = idx_s ;
            beats_seg(:,pp) = IHPks(idx_s);
        end
    end
    %===============================================
    % local optimal shrinkage
   sigma = std(rawEEG(~isnan(rawEEG))) ;
    beta = (L_1+L_2+1) ./ nn ;
    [u,l,v] = svd(beats_seg./(sqrt(nn)*sigma)); J = zeros(size(l));
    y = diag(l) ; eta = sqrt( (y.^2-beta-1).^2 - 4*beta) ./ y ;
    
    tmp = find(y<=1+sqrt(beta)) ; eta(tmp) = 0 ;
    tmp = min(size(J)) ; J(1:tmp,1:tmp) = diag(eta) ;
    Xc_s = (sigma*sqrt(nn))*u*J*v' ;
    %% put the Optimal shrinkage of tje beats back to its place in time
    if rr ~= group_len
        nn = ttt;
        for pp = 1:ttt
            tt = (rr-1)*nn+pp;
            X_ops(:,Idnx_1(tt)) = Xc_s(:,pp);
        end
        
    else
        nn = length(Idnx_1) - (group_len-1)*ttt ;
        for pp = 1:nn
            tt = (rr-1)*ttt+pp;
            X_ops(:,Idnx_1(tt)) = Xc_s(:,pp);
        end
    end
end 
    %%
    yyy = 250;
    ovlp = 250;
    group_len = floor(Nbeat/yyy);
    X = beatsHPks';
    beatHeight_OS = max(X(:,L_1-0.02*fs:L_1+0.02*fs),[],2);
    for ppp = 1:group_len
     if ppp ~= group_len 
        nn = yyy+ovlp;
        beatWidth_seg = zeros(1,nn);
        beatHeight_OS_seg = zeros(1,nn);
        beats_seg = zeros(L_1+L_2+1, nn) ;
        beats_OS_seg = zeros(L_1+L_2+1, nn) ;
        beatsidx_seg =  zeros(L_1+L_2+1, nn) ;
        try
        for pp = 1:yyy+ovlp
            tt = (ppp-1)*yyy+pp;
            idx_s = Loc(tt)-L_1: Loc(tt) + L_2  ;
            beatWidth_seg(:,pp)  = beatWidth(tt);
            beatHeight_OS_seg(:,pp) = beatHeight_OS(tt);
            beatsidx_seg(:,pp) = idx_s ;
            beats_seg(:,pp) = I(idx_s);
            beats_OS_seg(:,pp) = X_ops(:,tt);
        end
        catch
            keyboard
        end
     else
        nn = Nbeat - (group_len-1)*yyy ;
        beatWidth_seg = zeros(1,nn);
        beatHeight_OS_seg = zeros(1,nn);
        beats_seg = zeros(L_1+L_2+1, nn) ;
        beats_OS_seg = zeros(L_1+L_2+1, nn) ;
        beatsidx_seg =  zeros(L_1+L_2+1, nn) ;
        for pp = 1:nn
            tt = (ppp-1)*yyy+pp;
            idx_s = Loc(tt)-L_1: Loc(tt) + L_2  ;
            beatWidth_seg(:,pp)  = beatWidth(tt);
            beatHeight_OS_seg(:,pp) = beatHeight_OS(tt);
            beatsidx_seg(:,pp) = idx_s ;
            beats_seg(:,pp) = I(idx_s);
            beats_OS_seg(:,pp) = X_ops(:,tt);
        end
     end
     X_s = beats_OS_seg';
     %%
    [index, distance] = knnsearch(X_s,X_s, 'k', N) ;
    J1 = zeros(N,nn);
    
    for k = 1 : N; J1(k,:) = 1:nn; end
    J1 = J1(:);

    J2 = index.'; J2 = J2(:); DD = distance' ;
    Z = sparse(J1, J2, DD(:), nn, nn, nn*N) ;
    %===============================================
    % prepare data for diffusion distance
    tmp = quantile(distance(:), .5) ;
    QQ = 1./(1+(distance/tmp).^2) ; QQ = QQ' ;
    WW = sparse(J1, J2, QQ(:), nn, nn, nn*N) ;
    %===============================================
    % run diffusion distance
    if DIFFUSION
        fprintf('\t*** Usediffusion..\n') ;
 
        WW = sparse(WW'*WW) ;
        DD = sum(WW,2);

        Dinv = sparse(1:length(DD), 1:length(DD), 1./sqrt(DD)) ;
        AA = Dinv*WW*Dinv ;

        [UU,LL] = eigs(AA, 50) ;
        UU = Dinv * UU ;
        XX = UU(:, 2:end) * LL(2:end, 2:end) ;

        [indexD, distanceD] = knnsearch(XX, XX, 'k', N) ;
        J2 = indexD.'; J2 = J2(:); DDD = distanceD' ;
        Z = sparse(J1, J2, DDD(:), nn, nn, nn*N) ;
        toc(t0)
    end

    %===============================================
    % prepare weight information for NLEM
    ZW = zeros(N,nn) ; ZH = zeros(N,nn) ;
    for i = 1:nn
        ZW(:,i) = abs(beatWidth_seg(index(i,:)) - beatWidth_seg(i)) ;
        ZH(:,i) = abs(beats_seg(index(i,:)) - beats_seg(i)) ;
    end
    Zwidth = sparse(J1, J2, ZW, nn, nn, nn*N) ;
    Zheight = sparse(J1, J2, ZH, nn,nn, nn*N) ;
    clear J1; clear J2;
    %===============================================
    fprintf('\tRun NLEM...\n') ; t0=tic;
    fprintf(['Total: ',num2str(Nbeat), ' beats: ']) ;
    fprintf('%05d',0) ;
    
    xd0 = zeros(size(beatsidx_seg)) ;
    if ppp ~= group_len ; nn = yyy;end
    for i = 1:nn
        if ~mod(i,100) ; fprintf('\b\b\b\b\b') ; fprintf(['%05d'],i) ; end
        
        % don't take the beat itself into account
        tmp = Z(i,:) ; tmp(find(tmp==0)) = inf ;
        z = sort(tmp,'ascend');
        th = z(N-1) ;
        thTemp(i) = th;
        if isinf(th)
            Nidx = find(~isinf(tmp)) ;
        else
            % don't count the beat itself to avoid oversmoothing
            Nidx = find(tmp<=th) ;
        end
        Vx = beats_seg(:,Nidx);
        m = median(Vx,2);
        xd0(:, i) = m ;
    end
    for i = 1:nn
        tt = (ppp-1)*yyy+i;
        all_pk(:,tt) = xd0(:,i);
    end
    clear nn
    end