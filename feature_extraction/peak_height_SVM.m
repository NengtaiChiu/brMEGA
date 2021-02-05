function [pk_h] = peak_height_SVM(EEG,sampling_rate,IHR,pad)
if length(EEG)~= length(IHR)
    %keyboard;
end
f = IHR;
  %% find peak
         locEEG1 = beat_simple((EEG)',sampling_rate,(f.*0.02),20);
         peak = zeros(1,length(f));
         for k = locEEG1
             peak(k) = 1;
         end
         %%
        clean = EEG;
        for i  = 2:length(locEEG1)-1
            clean(locEEG1(i)-0.06*sampling_rate:locEEG1(i)+0.06*sampling_rate) = 0;
        end
         %% remove padding
            EEG =EEG(sampling_rate*pad+1:end-sampling_rate*pad);
            p =peak(sampling_rate*pad+1:end-sampling_rate*pad);
            clean  = clean(sampling_rate*pad+1:end-sampling_rate*pad);
            %% find relative peak height
            [~,mean,sigma] = zscore(clean);
             qq = [];rr = find(p == 1);
        for ww = 1:length(rr)
            relative_height = (EEG(rr(ww))-mean)./sigma;
            qq = [qq;relative_height];
        end
        pk_h = median(qq);
        %keyboard
