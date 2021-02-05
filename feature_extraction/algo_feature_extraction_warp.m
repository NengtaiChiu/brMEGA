function [peak_height,ot_s,ot_m]=algo_feature_extraction_warp(X,epochs,sampling_rate)
%%
% X = input signal;
% epochs = the 20s epochs that needs feature extraction for every 4s
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
    P.num_s = 1; % (harmonic 的根數, 只做簡單探討的話設1就好)
    P.num_c = 1;

    % parameters for cepstral representation
    cepR.g = 0.3; % for generalized cepstrum: 
    cepR.Tc=0; % Global threshold of cepstrum

    lam_curve = 10 ;
    lam_beat = 20;
    num_nonlocal = 10;
   %% feature extraction
    epochs_num = length(epochs);%numbers of 20sec epochs
    pad = 10 ; % padding length
    peak_height = zeros(epochs_num,5);
    ot_s = zeros(epochs_num,5);
    ot_m = zeros(epochs_num,5);
    parfor i = 1:epochs_num
        t1 = 20*(epochs(i)-1);
        pk_i = zeros(1,5);
        ots_i = zeros(1,5);
        otm_i = zeros(1,5);
        for j = 1:5
            st = 4*(j-1);
            ori = X((t1+st-pad)*sampling_rate+1:(t1+st+4+pad)*sampling_rate);
            [dEEG] = EEG_preprocess_relu(ori,sampling_rate,pad,4);
            %% correct polarity
            [~, ~, ~, ConceFT, tfrsqtic] = ConceFT_sqSTFT_C(dEEG,...
        0,advTF.HighFreq,0.02/basicTF.fs, basicTF.hop, basicTF.win, 1, 6, 1, 1, 0) ;
            [f,~,~,~] = IHR(ConceFT,tfrsqtic,sampling_rate,basicTF.fr,pad);
            f_250 = repmat(f',1,basicTF.hop);
            [dEEG,~] = EEG_preprocess_relu_correct_polarity(ori,f_250,sampling_rate,pad,4);
            %% warping
        loc = beat_simple(dEEG',sampling_rate,(f_250.*0.02),30);
        T = loc - loc(1) ;
        T(1) = [];
        dEEG(1:loc(1)) = [];
        
        V = 2*pi*[1:length(T)];
        phi = interp1(T,V,[1:length(dEEG)],'pchip','extrap');
        dEEG_unwrap = interp1(phi*length(dEEG)/max(phi),dEEG,[1:length(dEEG)],'pchip','extrap');
        dEEG_unwrap(1:100) = 0;
        dEEG_unwrap = [zeros(1,loc(1)), dEEG_unwrap];
            %% oridinary SST
            [~, ~, ~, ConceFT, tfrsqtic] = ConceFT_sqSTFT_C(dEEG_unwrap',...
        0,advTF.HighFreq,0.02/basicTF.fs, basicTF.hop, basicTF.win, 1, 6, 1, 0, 0) ;
            [f,c1,c2,c3] = IHR(ConceFT,tfrsqtic,sampling_rate,basicTF.fr,pad);
            %% feature extraction
            f_250 = repmat(f',1,basicTF.hop);
            [pk_h_j] = peak_height_SVM(dEEG_unwrap',sampling_rate,f_250,pad);
            [OT_s_j,~] = OT_single_4sec(ConceFT,f,pad,4);            
            [OT_m_j,~] = OT_multiples_4sec(ConceFT,f,c1,c2,c3,pad,4);
            %%
            pk_i(j)= pk_h_j;
            ots_i(j) = OT_s_j;
            otm_i(j) = OT_m_j;
        end
        peak_height(i,:) = pk_i;
        ot_s(i,:) = ots_i;
        ot_m(i,:) = otm_i;
    end