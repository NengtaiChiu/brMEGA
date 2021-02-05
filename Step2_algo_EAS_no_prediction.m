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
cfolder ='output_EAS_no_prediction';
mkdir(cfolder);   
%%
for jj = 1:length(ggg)
     load([foldername,'/',ggg(jj).name]) ;
    fprintf(['Analyze ',foldername,'/',ggg(jj).name,'\n']) ;
   %% prediction result: 6 colums 
    ep = 1:floor(length(EEG_FpzA2)/sampling_rate/20);% index of 20 sec epochs 
    pred = zeros(length(ep),5); % prediction outcome for each 4s subepoch
    pred (:,:)= 2; % predict all epochs
    %% basic parameters for STFT
    basicTF.win = sampling_rate*10+1;
    basicTF.hop = sampling_rate/5; %note that sampling_rate is 250, so the hopping frequency is 5hz
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
    len = length(ep);
    pd2_num = length(find(pred == 2)); % number of subepochs with predction value == 2
    pad = 2;
    EEG_pd2 = zeros(pd2_num*(4+2*pad)*sampling_rate,1);
    IHR_pd2 = zeros(pd2_num*(4+2*pad)*sampling_rate,1);
    correct_polarity = zeros(length(ep),5);
    count  = 0;
    
    locEEG = [];% peak location in the whole signal
    locEEG_pd2 = []; % peak location in EEG_pd2
    for i = 2:len-1 %the i-th 20 sec epoch
        for j = 1:5 % the j-th sub epoch in the i-th 20 sec epoch
            t = ep(i)-1;
            if pred(i,j) == 2
                count = count+1;
                st = 4*(j-1);
                ori = EEG_FpzA2((20*t+st-pad)*sampling_rate+1:(20*t+st+4+pad)*sampling_rate);
                %% detrend
                [dEEG] = EEG_preprocess_relu(ori,sampling_rate,pad,4,Dnotchfilter);
                
                [~, ~, ~, ConceFT, tfrsqtic] = ConceFT_sqSTFT_C(dEEG,...
                    0,advTF.HighFreq,0.02/basicTF.fs, basicTF.hop, basicTF.win, 1, 6, 1, 1, 0) ;
                [f,~,~,~] = IHR(ConceFT,tfrsqtic,sampling_rate,basicTF.fr,pad);
                % recover IHR at the original sampling rate
                G = repmat(f, 1, basicTF.hop) ; G = G' ;
                f_250 = G(:)' ;
                 % handle polarity
                [dEEG_r,correct] = EEG_preprocess_relu_correct_polarity(ori,f_250,sampling_rate,pad,4,Dnotchfilter);
                correct_polarity(i,j) = correct;
                if correct ==1
                    [dEEG] = EEG_preprocess_detrend(-ori,sampling_rate,pad,4,Dnotchfilter);
                else
                    continue
                end
                %% locate beats
                
                locEEG1 = beat_simple(dEEG_r',sampling_rate,(f_250.*0.02),30);
                
                a = find(locEEG1<(pad+4)*sampling_rate&locEEG1>pad*sampling_rate); 
                a(dEEG(locEEG1(a))< 3) = [];
                
                locEEG_i = locEEG1(a)- pad*sampling_rate+(20*t+st)*sampling_rate;
                locEEG = [locEEG;locEEG_i']; % peak location in the whole signal
                
                locEEG_j = locEEG1(a) + (2*pad+4)*(count-1)*sampling_rate;
                locEEG_pd2 = [locEEG_pd2;locEEG_j']; % peak location in EEG_pd2
                
                %% taper the padding signal with window function
                W = ones(length(dEEG),1);
                left_overlap = pad*0.5*sampling_rate;
                right_overlap = pad*0.5*sampling_rate;
                W(1:left_overlap) = sin(linspace(0,pi/2,left_overlap)).^2;
                W(end:-1:end-right_overlap+1) = sin(linspace(0,pi/2,right_overlap)).^2;
                dEEG = dEEG.*W;
                %%
                EEG_pd2((count-1)*sampling_rate*(4+2*pad)+1:(count)*sampling_rate*(4+2*pad)) = dEEG;
            end
        end
    end
    
    %% PJP method of removing artifact
     [xd] = artifact_remove_Hae_Jeong(EEG_pd2,locEEG_pd2,sampling_rate);
    %% remove artifact from original signal
    len = length(ep);
    estimated_artifact = zeros(length(EEG_FpzA2),1);
    count  = 0;
    pad=2;
    for i = 2:len-1
        for j = 1:5
            t = ep(i)-1;
            st = 4*(j-1);
            if pred(i,j) == 2
                count = count+1;
                xd_ij = xd((count-1)*sampling_rate*(4+2*pad)+sampling_rate*pad+1:(count)*sampling_rate*(4+2*pad)-sampling_rate*pad);
                if correct_polarity(i,j) == 1
                    xd_ij = -xd_ij;
                end
                estimated_artifact((20*t+st)*sampling_rate+1:(20*t+st+4)*sampling_rate) = xd_ij;
            end
        end
    end
    
    EEG_all_clean =   EEG_FpzA2 - estimated_artifact;
    %%
    f1 = [cfolder,'/',ggg(jj).name(1:end-4),'_EAS_no_prediction.mat'];
%   save(f1,'EEG_FpzA2','sampling_rate','ECG_scoring','EEG_all_clean','estimated_artifact','ep','pred','correct_polarity','locEEG','EEG_pd2','locEEG_pd2');
   save(f1,'EEG_FpzA2','sampling_rate','ECG_scoring','EEG_all_clean','estimated_artifact','ep','pred','correct_polarity','locEEG');
end
