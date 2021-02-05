function [xd] = EAS_method(I, Loc, sampling_rate,alpha)
xd = zeros(size(I));
Nbeat = length(Loc) ;
L = 1.5*sampling_rate;

fprintf(['\tOnly ',num2str(Nbeat),' beats\n']) ;%keyboard
beats = zeros(L+1,Nbeat);
beats_mod = zeros(L+1,Nbeat);
for ii = 2:Nbeat-1
    loc_s = Loc(ii) - 0.2*sampling_rate;
    loc_t = Loc(ii+1) - 0.2*sampling_rate;
    beat_len  = loc_t - loc_s;
    if beat_len < 1.2 * sampling_rate
        beats(1:beat_len,ii) = I(loc_s:loc_t-1);
    else
        beats(1:0.6*sampling_rate+1,ii) = I(loc_s:loc_s+0.6*sampling_rate);
    end
end

%% split beats into groups
N_100 = floor(Nbeat/100);
for i = 1:N_100
    if i == 1
        seg_beats = beats(1:end,1:100);
        avg_beats = mean(seg_beats,2);
        beats_mod(1:length(avg_beats),1:100) = repmat(avg_beats,1,100);
    else
        seg_beats = beats(1:end,100*(i-1)+1:100*i); m = mean(seg_beats,2);
        avg_beats = avg_beats.*alpha+m.*(1-alpha);
        beats_mod(1:length(avg_beats),100*(i-1)+1:100*i) = repmat(avg_beats,1,100);
    end
end
%% determine Time delay
for ii = 2:Nbeat-1
    loc_s = Loc(ii) - 0.2*sampling_rate;
    loc_t = Loc(ii+1) - 0.2*sampling_rate;
    beat_len  = loc_t - loc_s;
    if beat_len < 1.2 * sampling_rate
        xd_i = I(loc_s:loc_t-1);
        r = xcorr(xd_i,beats_mod(1:beat_len,ii));
        [~,k] = max(r);
        delay = k - beat_len;
        delay = 0;
        xd(loc_s+delay:loc_t+delay-1) = beats_mod(1:beat_len,ii);
    else
        beat_len = 0.6*sampling_rate;
        xd_i = I(loc_s:loc_s+beat_len);
        r = xcorr(xd_i,beats_mod(1:beat_len,ii));
        [~,k] = max(r);
        delay = k - beat_len;
        delay = 0;
        xd(loc_s+delay:loc_s+beat_len+delay-1) = beats_mod(1:beat_len,ii);
    end
end
end