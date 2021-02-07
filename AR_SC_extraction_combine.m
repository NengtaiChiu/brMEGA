clear all;close all;
addpath('feature_extraction');
addpath('preprocess_method');
addpath('TF_anaylsis');

foldername = 'output_brMEGA_aR_sC';
ggg = dir([foldername,'/*.mat']) ;
%%
sampling_rate = 250;
d = designfilt('bandstopiir','FilterOrder',2, ...
    'HalfPowerFrequency1',49,'HalfPowerFrequency2',51, ...
    'DesignMethod','butter','SampleRate',sampling_rate);
mkdir 'output_aR_sC' ;
%%
aR_brMEGA = [];sC_brMEGA = [];patient_name =[];
for jj = 1:length(ggg)
    load([foldername,'/',ggg(jj).name]) ;
    fprintf(['Analyze ',foldername,'/',ggg(jj).name,'\n']) ;%keyboard
    %%
    EEG_FpzA2 = filtfilt(d,EEG_FpzA2);
    window = 0.1;%seconds
    K = floor(window/2*sampling_rate);
    [trend]=median_filter(EEG_FpzA2,K);
    trend = smooth(trend,3*K,'loess');%smooth the trend
    EEG = EEG_FpzA2 -trend;
    EEG_c = EEG - estimated_artifact;
    sC_all = [];aR_all = [];
    %%
    for i = 1:length(ep_20)-1
        s = (ep_20(i)-1)*20*sampling_rate; % starting time
        t = ep_20(i)*20*sampling_rate; % ending time 
        
        EEG_r_seg = EEG(s+1:t); % contaminated EEG
        EEG_c_seg = EEG_c(s+1:t); % cleaned EEG
        a = find(locEEG<t & locEEG>s); % get peak location
        loc = locEEG(a)-s;
        
        T = loc - loc(1) ;
        T(1) = [];
        EEG_r_seg(1:loc(1)) = [];
        EEG_c_seg(1:loc(1)) = [];
        
        V = 2*pi*[1:length(T)];
        phi = interp1(T,V,[1:length(EEG_r_seg)],'pchip','extrap');
        raw_unwrap = interp1(phi*length(EEG_r_seg)/max(phi),EEG_r_seg,[1:length(EEG_r_seg)],'pchip','extrap');
        raw_unwrap(1:100) = 0;
        
        clean_unwrap = interp1(phi*length(EEG_c_seg)/max(phi),EEG_c_seg,[1:length(EEG_c_seg)],'pchip','extrap');
        clean_unwrap(1:100) = 0;
        %%
        warp_peak_loc = round(2*pi*[1:length(T)]*length(EEG_r_seg)/max(phi));
        [aR,sC,~] = nlem_performance(clean_unwrap, raw_unwrap, warp_peak_loc);
        aR_all = [aR_all;aR];
        sC_all = [sC_all;sC];
    end
    %%
     % k = {'S025','S027','S032','S033','S036','S039','S041','S008','S011','S015','S016'};
     aR_all = log10(aR_all);

     aR_all_m =rmoutliers(aR_all,'percentiles',[1 99]);
     sC_all_m =rmoutliers(sC_all,'percentiles',[1 99]);
     rr = length(aR_all_m );
      patient_name_i = repelem(string(ggg(jj).name(1:4)), rr)';
      patient_name = [patient_name ; patient_name_i];
     aR_brMEGA = [aR_brMEGA ;aR_all_m];sC_brMEGA = [sC_brMEGA ;sC_all_m];
     
end

%%
foldername = 'output_EAS_aR_sC';
fff = dir([foldername,'/*.mat']) ;
aR_EAS = [];sC_EAS = [];
%%
for jj = 1:length(fff)
    load([foldername,'/',fff(jj).name]) ;
    fprintf(['Analyze ',foldername,'/',fff(jj).name,'\n']) ;%keyboard
    %%
    EEG_FpzA2 = filtfilt(d,EEG_FpzA2);
    window = 0.1;%seconds
    K = floor(window/2*sampling_rate);
    [trend]=median_filter(EEG_FpzA2,K);
    trend = smooth(trend,3*K,'loess');%smooth the trend
    EEG = EEG_FpzA2 -trend;
    EEG_c = EEG - estimated_artifact;
    sC_all = [];aR_all = [];
    %%
    for i = 1:length(ep_20)-1
        s = (ep_20(i)-1)*20*sampling_rate; % starting time
        t = ep_20(i)*20*sampling_rate; % ending time 
        
        EEG_r_seg = EEG(s+1:t); % contaminated EEG
        EEG_c_seg = EEG_c(s+1:t); % cleaned EEG
        a = find(locEEG<t & locEEG>s); % get peak location
        loc = locEEG(a)-s;
        
        T = loc - loc(1) ;
        T(1) = [];
        EEG_r_seg(1:loc(1)) = [];
        EEG_c_seg(1:loc(1)) = [];
        
        V = 2*pi*[1:length(T)];
        phi = interp1(T,V,[1:length(EEG_r_seg)],'pchip','extrap');
        raw_unwrap = interp1(phi*length(EEG_r_seg)/max(phi),EEG_r_seg,[1:length(EEG_r_seg)],'pchip','extrap');
        raw_unwrap(1:100) = 0;
        
        clean_unwrap = interp1(phi*length(EEG_c_seg)/max(phi),EEG_c_seg,[1:length(EEG_c_seg)],'pchip','extrap');
        clean_unwrap(1:100) = 0;
        %%
        warp_peak_loc = round(2*pi*[1:length(T)]*length(EEG_r_seg)/max(phi));
        [aR,sC,~] = nlem_performance(clean_unwrap, raw_unwrap, warp_peak_loc);
        aR_all = [aR_all;aR];
        sC_all = [sC_all;sC];
    end
    %%
     aR_all = log10(aR_all);

     aR_all_m =rmoutliers(aR_all,'percentiles',[1 99]);
     sC_all_m =rmoutliers(sC_all,'percentiles',[1 99]);
   
     aR_EAS = [aR_EAS ;aR_all_m];sC_EAS = [sC_EAS ;sC_all_m];
end
     %%
     k = {'patient';'brMEGA';'EAS'};
     T = table(patient_name,aR_brMEGA,aR_EAS,'VariableNames',k);
     fname = ['output_aR_sC/','aR_extract_all','.xlsx'];
     writetable( T,fname);
     
     T = table(patient_name,sC_brMEGA,sC_EAS,'VariableNames',k);
     fname = ['output_aR_sC/','sC_extract_all','.xlsx'];
     writetable( T,fname);