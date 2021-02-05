clear all;close all;
addpath('feature_extraction');
addpath('libsvm');
addpath('libsvm/matlab');
addpath('preprocess_method');
addpath('TF_anaylsis');

foldername = 'database';
ggg = dir([foldername,'/*.mat']) ;
%%
sampling_rate = 250 ;
% notch filter to remover the powerline interference 
Dnotchfilter = designfilt('bandstopiir','FilterOrder',2, ...
    'HalfPowerFrequency1',49,'HalfPowerFrequency2',51, ...
    'DesignMethod','butter','SampleRate',sampling_rate);
mkdir 'output_prediction_result';
%%
for jj = 1:length(ggg)
    load([foldername,'/',ggg(jj).name]) ;
    load('ALGO_train.mat');% load the training model from step 1
    fprintf(['Analyze ',foldername,'/',ggg(jj).name,'\n']) ;
    %% feature extraction
    epochs = 2:floor(length(EEG_FpzA2)/sampling_rate/20)-1;
    [peak_height,ot_s,ot_m,ot_power]=algo_feature_extraction(EEG_FpzA2,epochs,sampling_rate,Dnotchfilter);
    
    mean = Zscore.mean;sigma = Zscore.sigma;% apply the zscore from the training data
    peak_height = (peak_height-mean(1))/sigma(1);
    ot_s = (ot_s-mean(2))/sigma(2);
    ot_m = (ot_m-mean(3))/sigma(3);
    ot_power = (ot_power-mean(4))/sigma(4);
    
    peak_height = reshape(peak_height,[],1);
    ot_s = reshape(ot_s,[],1);
    ot_m = reshape(ot_m,[],1);
    ot_power = reshape(ot_power,[],1);
    %%
    [prediction]=svmpredict(zeros(length(peak_height),1),[peak_height,ot_s,ot_m,ot_power],train_model,'-q');
    prediction = reshape(prediction,[],5);
    %%
    prediction_result = [epochs',prediction];
    filename = ['output_prediction_result/','prediction_',ggg(jj).name];
    save(filename,'prediction_result');
end
