function [X] = EEG_preprocess_detrend_hp(rawEEG,sampling_rate,cut,pad,length)
mn = mean(rawEEG(sampling_rate*pad+1:sampling_rate*(pad+length)));
EEG = rawEEG-mn;
%% notch filter to remove powerline interference
d = designfilt('bandstopiir','FilterOrder',2, ...
    'HalfPowerFrequency1',49,'HalfPowerFrequency2',51, ...
    'DesignMethod','butter','SampleRate',sampling_rate);
EEG = filtfilt(d,EEG);
%% median filter
[b_hp,a_hp] = butter(5, cut/(sampling_rate/2),'high');
EEG = filtfilt(b_hp,a_hp, EEG);
X = smooth(EEG, 'loess', sampling_rate*0.005) ;
