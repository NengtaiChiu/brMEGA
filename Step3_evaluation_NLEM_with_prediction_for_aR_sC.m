clear all;close all;
addpath('feature_extraction');
addpath('libsvm');
addpath('preprocess_method');
addpath('TF_anaylsis');

foldername = 'database';

ggg = dir([foldername,'/*.mat']) ;

%%
sampling_rate = 250 ;
Dnotchfilter = designfilt('bandstopiir','FilterOrder',2, ...
    'HalfPowerFrequency1',49,'HalfPowerFrequency2',51, ...
    'DesignMethod','butter','SampleRate',sampling_rate);
cfolder = 'output_NLEM_aR_sC';
mkdir(cfolder);

for jj = 1:length(ggg)
    load([foldername,'/',ggg(jj).name]) ;
    fprintf(['Analyze ',foldername,'/',ggg(jj).name,'\n']) ;
    load(['output_prediction_result/','prediction_',ggg(jj).name(1:end)]);
    %%
    ep = prediction_result(:,1);% number of 20 sec epochs with prediction
    pred = prediction_result(:,2:end);
    ep_20 = [];
    for i = 1:length(ep)
        a = find(pred(i,:) >= 2);
        if length(a) >= 3
            ep_20 = [ep_20;ep(i)];
        end
    end
    if length(ep_20)<= 30
        continue
    end
    %% basic parameters for STFT
    basicTF.win = sampling_rate*10+1;
    basicTF.hop = 50; %note that sampling_rate is 250, so the hopping frequency is 5hz
    basicTF.fs = sampling_rate ;
    basicTF.fr = 0.02; % frequency resolution
    basicTF.feat = 'SST11'; % two option: STFT or SST11 (one window rejection)
    
    % advanced parameters for STFT
    advTF.num_tap = 1; % Number of tap in ConceFT
    advTF.win_type = 'Gauss'; % Only 2-tap basis here
    advTF.Smo = 1; % alpha smoothing; 1 = no smoothing
    advTF.Rej = 0; % The bandwidth of window rejection;
    advTF.ths = 1E-9; % Global threshold of STFT
    advTF.HighFreq = 8/basicTF.fs; % highest frequency/sampling freq
    advTF.LowFreq = 0.02/basicTF.fs; % lowest frequency/sampling freq
    advTF.lpc = 0;
 
    %% create signal for peak removal
    len = length(ep_20);
    pad = 10;
    EEG_pd2 = zeros(len*(20+2*pad)*sampling_rate,1);
    IHR_pd2 = zeros(len*(20+2*pad)*sampling_rate,1);
    count  = 0;
    
    locEEG = [];% peak location in the whole signal
    locEEG_pd2 = []; % peak location in EEG_pd2
    for i = 1:len %the i-th 20 sec epoch
        t = ep_20(i)-1;
        count = count +1;
        ori = EEG_FpzA2((20*t-pad)*sampling_rate+1:(20*(t+1)+pad)*sampling_rate);
        %% 
        [dEEG] = EEG_preprocess_relu(ori,sampling_rate,pad,20,Dnotchfilter);
        
        [~, ~, ~, ConceFT, tfrsqtic] = ConceFT_sqSTFT_C(dEEG,...
            0,advTF.HighFreq,0.02/basicTF.fs, basicTF.hop, basicTF.win, 1, 6, 1, 1, 0) ;
        
        [f,~,~,~] = IHR(ConceFT,tfrsqtic,sampling_rate,basicTF.fr,pad);
        G = repmat(f, 1, basicTF.hop) ; G = G' ;
        f_250 = G(:)' ;
        
        %% locate beats
        locEEG1 = beat_simple(dEEG',sampling_rate,(f_250.*0.02),30);
        
        a = find(locEEG1<(pad+20)*sampling_rate&locEEG1>pad*sampling_rate);
        a(dEEG(locEEG1(a))< 3) = [];
        
        locEEG_i = locEEG1(a)- pad*sampling_rate+(20*t)*sampling_rate;
        locEEG = [locEEG;locEEG_i']; % peak location in the whole signal
        
        locEEG_j = locEEG1(a) + (2*pad+20)*(count-1)*sampling_rate;
        locEEG_pd2 = [locEEG_pd2;locEEG_j']; % peak location in EEG_pd2
        
        %% taper the padding signal with window function
        [dEEG] = EEG_preprocess_detrend(-ori,sampling_rate,pad,20,Dnotchfilter);
        
        W = ones(length(dEEG),1);
        left_overlap = 2*sampling_rate;
        right_overlap = 2*sampling_rate;
        W(1:left_overlap) = sin(linspace(0,pi/2,left_overlap)).^2;
        W(end:-1:end-right_overlap+1) = sin(linspace(0,pi/2,right_overlap)).^2;
        dEEG = dEEG.*W;
        %%
        EEG_pd2((count-1)*sampling_rate*(20+2*pad)+1:(count)*sampling_rate*(20+2*pad)) = dEEG;
    end
    %% NLEM
        [all_pk] = NLEM_new(EEG_pd2,locEEG_pd2, 100, sampling_rate, 5, 10 ,1, 1);
%%
    %locEEG = ori_loc;
  StimArt = zeros(size(EEG_FpzA2)) ;
    
    pL = 0.2* sampling_rate ;
    pR = 0.5* sampling_rate ;
    
    BDRY = 8 ;
    WIN0 = ones(pL+pR+1, 1) ;
    WIN0(1: BDRY) = sin(pi*[0:BDRY-1]'/2/BDRY).^2 ;
    WIN0(end-BDRY+1: end) = cos(pi*[1:BDRY]'/2/BDRY).^2 ;
    
    % truncate to avoid the boundary effect
    all_pk2 = diag(WIN0) * all_pk ;
     

    beats2 = locEEG;
    %ori_loc = locEEG;
    rm = [];
    for i = 1:length(locEEG)-1
        RRI = locEEG(i+1)-locEEG(i);
        if RRI < 0.4*sampling_rate
            [~,a] = min([max(all_pk2(:,i)),max(all_pk2(:,i+1))]);
            rm = [rm;(i+a-1)];
        end
    end
   beats2(rm) = [];
   all_pk2(:,rm) = [];

    % ONE is for the debugging purpose.
    ONE = zeros(size(EEG_FpzA2)) ;

    for qq = 2: length(beats2)-1 %beat2 = peakloc
        
        
        BDRYL = max(1, pL+pR+2 - (beats2(qq) - beats2(qq-1) + 1) - 1) ;
        BDRYR = max(1, pL+pR+2 - (beats2(qq+1) - beats2(qq) + 1) - 1) ;
        
        
        WIN = ones(1, pL+pR+1) ;
        
        
        WIN(1: BDRYL) = sin(pi*[0:BDRYL-1]'/2/BDRYL).^2 ;
        
        
        WIN(end-BDRYR+1: end) = cos(pi*[1:BDRYR]'/2/BDRYR).^2 ;
        
        
        tmp = all_pk2(:, qq)' .* WIN ;
        
        
        idx = beats2(qq) - pL : beats2(qq) + pR ;
        StimArt(idx) = StimArt(idx) + tmp' ;
        ONE(idx) = ONE(idx) + WIN' ;
    end
    estimated_artifact = -StimArt;
    EEG_all_clean =   EEG_FpzA2 - estimated_artifact;
    locEEG = beats2;
    %%
    f1 = [cfolder,'/',ggg(jj).name(1:end-4),'_NLEM_prediction_evaluation.mat'];
%   save(f1,'ori_loc','EEG_FpzA2','sampling_rate','EEG_all_clean','estimated_artifact','ep_20','prediction_result','locEEG','EEG_pd2','locEEG_pd2','all_pk');
  save(f1,'EEG_FpzA2','sampling_rate','EEG_all_clean','estimated_artifact','ep_20','prediction_result','locEEG');
end
