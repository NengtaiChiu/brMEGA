function [X,correct] = EEG_preprocess_detrend_correct_polarity(rawEEG,f,sampling_rate,pad,length,d)
    mn = mean(rawEEG(sampling_rate*pad+1:sampling_rate*(pad+length)));
    EEG = rawEEG-mn;
    %% notch filter to remove powerline interference
    EEG = filtfilt(d,EEG); 
     %% median filter
    window = 0.1; % seconds
    K = floor(window/2*sampling_rate);
    [trend]=median_filter(EEG,K);
    trend = smooth(trend,3*K,'loess');% smooth the trend
    dEEG = EEG -trend;
    %%    
    locEEG1 = beat_simple(dEEG',sampling_rate,(f.*0.02),30);
    locEEG2 = beat_simple(-dEEG',sampling_rate,(f.*0.02),30);
       if mean(dEEG(locEEG1)) > mean(-dEEG(locEEG2))
        X= dEEG ;
        correct = 0;
        else
        X= -dEEG ;
        correct = 1;
        end
