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
cfolder = 'output_EAS_aR_sC';
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
        %%  detrend 
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
    %% EAS method
      [xd] = artifact_remove_Hae_Jeong(EEG_pd2,locEEG_pd2,sampling_rate);% neighbors=50; (with diffusion = 1 ; No diffusion = 0)
    %% remove artifact from original signal
    len = length(ep_20);
    estimated_artifact = zeros(length(EEG_FpzA2),1);
    count  = 0;
    pad=10;
    for i = 1:len
        t = ep_20(i)-1;
        count = count+1;
        xd_i = xd((count-1)*sampling_rate*(20+2*pad)+sampling_rate*pad+1:(count)*sampling_rate*(20+2*pad)-sampling_rate*pad);
        estimated_artifact(20*t*sampling_rate+1:(20*t+20)*sampling_rate) = -xd_i;
    end
    EEG_all_clean =   EEG_FpzA2 - estimated_artifact;
    %%
    f1 = [cfolder,'/',ggg(jj).name(1:end-4),'_EAS_prediction_evaluation.mat'];
   save(f1,'EEG_FpzA2','sampling_rate','EEG_all_clean','estimated_artifact','ep_20','prediction_result','locEEG');
end
