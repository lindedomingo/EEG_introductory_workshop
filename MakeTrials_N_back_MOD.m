function [trl, event] = MakeTrials_ANT_Embodiment(cfg)

event = ft_read_event(cfg.dataset);


% searchg for stimulus events
value1 = [event(find(strcmp(char(cfg.trigval(1)),{event.value}))).value];
sample1 = [event(find(strcmp(char(cfg.trigval(1)),{event.value}))).sample];
% number of samples before and after onset
pretrig = -round(cfg.trialdef.pre * cfg.fsample);%need to find it from the data
posttrig = round(cfg.trialdef.post * cfg.fsample);

% look for all s2 triggers
trl1 = [];
for j = 1:length(sample1)
    trlbegin1 = sample1(j)+pretrig;
    trlend1 = sample1(j)+posttrig;
    offset1 = pretrig;
    newtrl1 = [trlbegin1 trlend1 offset1];
    trl1 = [trl1; newtrl1];
end

trl_all = [];
trl_all = [trl1];

[B,I] = sort(trl_all(:,1),1);

trl =[];
trl = trl_all(I,:);
