function power_iteration(me,X)

% Refine the feature and filter estimates obtained with HOSD through a power iteration method.
% The ability of HOSD to recover a moment-maximizing filter is only
% approximate, but it yields a starting estimate which can be quickly
% refined using a power iteration. This is implemented by updating the detection
% filter with the partial delay filter computed for the current feature estimate, 
% then recomputing the feature by aligning on the outputs of the new
% feature.


pdfilt =  me.partial_delay_filt(me.feature);

me.filterfft =pdfilt;
% delt = me.delay;

[~,FXsh] = me.apply_filter(X,false);

me.wavefft = nanmean(FXsh,2);

% me.delay = delt;
