% Script to calculate performance measures for ICA and NLM
clear
clc

load('C:\Users\FrohlichLab\Dropbox (Frohlich Lab)\Codebase\CodeSankar\ECoG_Project\Non Local Euclidean Median\Data\Subject_017_StimElec1_Post60Clean.mat')
load('Subject_017_StimSegments.mat', 'stimElec1Clean60')
w0 = 60/500; bw = w0/35;
[b,a] = iirnotch(w0,bw);
%%
rawData = stimElec1Clean60;
cleanDataNLM = cleanData1;
% rawData = rawData - repmat(nanmean(rawData,2),[1,size(rawData,2)]);
% cleanDataNLM = cleanDataNLM - repmat(nanmean(cleanDataNLM,2),[1,size(cleanDataNLM,2)]);
load('C:\Users\FrohlichLab\Dropbox (Frohlich Lab)\Codebase\CodeSankar\ECoG_Project\Non Local Euclidean Median\Data\S017_ICA_CleanData.mat','stimElec1CleanICA')
cleanDataICA = stimElec1CleanICA;
% cleanDataICA = cleanDataICA - repmat(nanmean(cleanDataICA,2),[1,size(cleanDataICA,2)]);
stimChannel = stimElec1(4,:);
[~,stimLoc] = findpeaks(stimChannel,'MinPeakHeight',1000, 'MinPeakDistance',75) ;

for i = 1:985
    stimInd(i,:) = stimLoc(i)-10:stimLoc(i)+10;
end
stimAmp = nanmean(rawData(:,stimInd),2);
stimAmp(4:5) = nan;
% Plot all data
rD = rawData;
rD(4:5,:) = NaN;
t = linspace(0,size(rawData,2)/1000,size(rawData,2));
plot_timeseries(rD,t,1:112,[],300)
%%
aR_ICA = zeros(size(rawData,1),1); aR_NLM = zeros(size(rawData,1),1);
sC_ICA = zeros(size(rawData,1),1); sC_NLM = zeros(size(rawData,1),1);

for iChannel = 1:size(rawData,1)
    [aR_ICA(iChannel),sC_ICA(iChannel),atICA(iChannel,:)] = nlem_performance(cleanDataICA(iChannel,:), rawData(iChannel,:), stimLoc);
    [aR_NLM(iChannel),sC_NLM(iChannel),atNLM(iChannel,:)] = nlem_performance(cleanDataNLM(iChannel,:), rawData(iChannel,:), stimLoc);
end
%% Spectra Calculation
for iChannel = 1:size(rawData,1)
    [powerRaw(iChannel,:) ,f] = pwelch(rawData(iChannel,:), 5000, [], 5000, 1000);
    powerCleanICA(iChannel,:) = pwelch(cleanDataICA(iChannel,:), 5000, [], 5000, 1000);
    powerCleanNLM(iChannel,:) = pwelch(cleanDataNLM(iChannel,:), 5000, [], 5000, 1000);
    
end

%% Plotting
% Channels 1,11,21 for high, medium, low
channelNumber = 21;
figure(1)
plot(t,rawData(channelNumber,:),t,cleanDataICA(channelNumber,:),t,cleanDataNLM(channelNumber,:))
set(gca,'fontsize',20)
xlabel('Time [s]')
ylabel('Amplitude [\muV]')
legend('Raw','ICA','NLM')
xlim([36.75 37.15])

figure(2)
plot(f,10*log10(powerRaw(channelNumber,:)),f,10*log10(powerCleanICA(channelNumber,:)),f,10*log10(powerCleanNLM(channelNumber,:)))
set(gca,'fontsize',20)
xlabel('Frequency [Hz]')
ylabel('Power [\muV^2]')
legend('Raw','ICA','NLM')
ylim([-15 30])
%%
x1 = 0.35*rand(114,1);
x2 = 0.35*rand(114,1);

aR_ICA(4:5) = NaN;
sC_ICA(4:5) = NaN;

figure(3)
plot(x1-0.2,log10(aR_ICA),'ko','markersize',10,'markerfacecolor','b')
hold on
plot(x2+1-0.2,log10(aR_NLM),'ko','markersize',10,'markerfacecolor','r')
set(gca,'fontsize',20)
set(gca,'Xtick',0:1,'Xticklabel',{'ICA','NLEM'})
xlim([-0.5 1.5])
ylim([0 1.5])

figure(4)
plot(x1-0.15,sC_ICA,'ko','markersize',10,'markerfacecolor','b')
hold on
plot(x2+1-0.15,sC_NLM,'ko','markersize',10,'markerfacecolor','r')
set(gca,'fontsize',20)
set(gca,'Xtick',0:1,'Xticklabel',{'ICA','NLEM'})
xlim([-0.5 1.5])
% ylim([0 1.5])
%%
mn_aR = [nanmean(log10(aR_ICA)) nanmean(log10(aR_NLM))];
sem_aR = [nanstd(log10(aR_ICA))./sqrt(112) nanstd(log10(aR_NLM))./sqrt(112)];
d_aR = (nanmean(log10(aR_NLM)) - nanmean(log10(aR_ICA)))./sqrt(nanstd(log10(aR_NLM)).^2 + nanstd(log10(aR_ICA)).^2 );

mn_sC = [nanmean(sC_ICA) nanmean(sC_NLM)];
sem_sC = [nanstd(sC_ICA)./sqrt(112) nanstd(sC_NLM)./sqrt(112)];
d_sC = (nanmean(sC_ICA) - nanmean(sC_NLM))./sqrt(nanstd(sC_NLM).^2 + nanstd(sC_ICA).^2 );
