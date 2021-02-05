clear all;close all;
addpath('feature_extraction');
addpath('libsvm');
addpath('libsvm/matlab');
addpath('preprocess_method');
addpath('TF_anaylsis');

foldername = 'database'; % folder name of the EEG signal...
labels_mat = 'data_label'; % folder name of the labels provided by the neuroscientist

ggg = dir([foldername,'/*.mat']) ;
%%
sampling_rate = 250 ;
% notch filter to remover the powerline interference 
Dnotchfilter = designfilt('bandstopiir','FilterOrder',2, ...
    'HalfPowerFrequency1',49,'HalfPowerFrequency2',51, ...
    'DesignMethod','butter','SampleRate',sampling_rate);

SVM_train_label = [];
SVM_train_feature = [];
for jj = 1:length(ggg)
    load([foldername,'/',ggg(jj).name]) ;
    fprintf(['Analyze ',foldername,'/',ggg(jj).name,'\n']) ;
    %% read label from CL filled mat file    
    ff = strcat(labels_mat,'/',string(ggg(jj).name(1:end-4)),'_labels_CL');
    load([ff]) ;
    % variable 'train_epochs'& 'train_label'is in the mat file loaded above
    %% feature extraction
    [peak_height,ot_s,ot_m,ot_power]=algo_feature_extraction(EEG_FpzA2,train_epochs,sampling_rate,Dnotchfilter);
    
    peak_height = reshape(peak_height',[],1);
    ot_s = reshape(ot_s',[],1);
    ot_m = reshape(ot_m',[],1);
    ot_power  = reshape(ot_power',[],1);
    train_label = reshape(train_label',[],1);
    
    SVM_train_feature =[SVM_train_feature;peak_height,ot_s,ot_m,ot_power];
    SVM_train_label = [SVM_train_label;train_label];
end
%%
    [SVM_train_feature,mean,sigma] = zscore(SVM_train_feature);
    Zscore.mean = mean;
    Zscore.sigma = sigma;
%% let label(1) be label (2), so SVM trains label(0) vs label(1)
zz = find(SVM_train_label ==  1);
SVM_train_label(zz) = 2;
%%
c = 1;g = (1/5);w1 = 1; w2 = 3;
cmd = ['-s 0 -q -h 1 -t 2 -m 16384 -c ' num2str(c)...
        ' -g ' num2str(g)...
        ' -b 0 -w1 ' num2str(w1) ' -w2 ' num2str(w2)  ];
train_model=svmtrain(SVM_train_label,SVM_train_feature,cmd);
save('ALGO_train','train_model','Zscore');