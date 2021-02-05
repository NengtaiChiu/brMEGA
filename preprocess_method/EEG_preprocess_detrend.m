function [X] = EEG_preprocess_detrend(rawEEG,sampling_rate,pad,length,d)
    mn = mean(rawEEG(sampling_rate*pad+1:sampling_rate*(pad+length)));
    EEG = rawEEG-mn;
    %% notch filter to remove powerline interference
    EEG = filtfilt(d,EEG); 
     %% median filter
    window = 0.1;%seconds
    K = floor(window/2*sampling_rate);
    [trend]=median_filter(EEG,K);
    trend = smooth(trend,3*K,'loess');%smooth the trend
    X = EEG -trend;
