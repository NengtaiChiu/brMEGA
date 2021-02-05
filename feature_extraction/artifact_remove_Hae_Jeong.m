function [xd] = artifact_remove_Hae_Jeong(X,peakloc,sampling_rate)
% X = input signal
% IHR = corresponding instantaneous heart rate
% N = number of neighbors for NLEM
%%  process signal for every 2 hours
len  = floor(length(X)/sampling_rate/60/60);
if len <= 1
    [xd] = EAS_method(I, peakloc, sampling_rate,0.3);
else
    xd = [];
    for i = 1:len-1
        X_i = X((i-1)*sampling_rate*60*60+1:i*sampling_rate*60*60);
        a = find(peakloc< i*sampling_rate*60*60 & peakloc>(i-1)*sampling_rate*60*60);
        locEEG = peakloc(a) - (i-1)*sampling_rate*60*60;
        %keyboard
        [xd_i] = EAS_method(X_i', locEEG, sampling_rate,0.3);
        xd = [xd;xd_i'];
    end
    %keyboard
    X_i = X((len-1)*sampling_rate*60*60+1:end);
    a = find( peakloc>(len-1)*sampling_rate*60*60);
    locEEG = peakloc(a) - (len-1)*sampling_rate*60*60;
    [xd_i] = EAS_method(X_i', locEEG, sampling_rate,0.3);
    xd = [xd;xd_i'];
end
end