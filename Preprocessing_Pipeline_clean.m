% PREPROCESSING PIPELINE

% Juan Linde Domingo, March 2024, CIMCYC
% EEG preprocessing - Course -- N-back task

%% PART I. Visual inspection and cleaning

%% Adding Fieldtrip path and preparing the environment

clear all
clc
close all

addpath /Users/juanlindedomingo/toolboxes/fieldtrip-20190224 %modify to your path
ft_defaults; %initialize fieldtrip with default values

%% Load EEG files
% Let us load the EEG data of the participant 32 and the layout of the EEG cap

load('EEG_dataset/raw_data/lay.mat') %layout file


cfg=[];
cfg.dataset= 'EEG_dataset/raw_data/p_32_nback.eeg'; %data file
cfg.headerfile   = 'EEG_dataset/raw_data/p_32_nback.vhdr'; %header file

% The data is loaded in a raw format. We need to convert it to a format that Fieldtrip can handle.

cfg.continuous = 'yes'; %continuous data
data_format = ft_preprocessing(cfg); %preprocess the data as continuous

%% Downsampling
% The data is sampled at 1000 Hz. We will downsample it to 250 Hz to speed up the analysis.

cfg=[];
cfg.resamplefs      = 250; %new sampling rate
[data_resampled] = ft_resampledata(cfg, data_format)

%% First visual check. 
% We will use the data browser to inspect the data, see how electrodes look like
% and identify bad channels. Also, let take a look to some general artifacts.
% We will write down the bad channels.

% It is very important to do this before any filtering or interpolation, as
% these steps can introduce changes that are not present in the original data.

% The point here is to have a general idea of the data quality.

cfg=[]
cfg.viewmode = 'vertical';
artf = ft_databrowser(cfg,data_resampled);


%% Inspect raw data. Kicking out bad channels
% Let check the statistics of the data and kick out bad channels.

% Summary
cfg          = [];
cfg.method   = 'summary';
dummy        = ft_rejectvisual(cfg,data_resampled);


%% Interpolate bad channels
% We will interpolate the bad channels using the average of the surrounding channels.
% We will use the function ft_channelrepair for this.
% The function requires that the data is in a specific format. We will convert the data to this format.

% Interpolation is a method that allows us to estimate the value of a variable at a point based on 
% the values of the surrounding points. In the case of EEG data, we can estimate the value of 
%   the signal at a bad channel based on the values of the surrounding channels.

data = dummy;

allLabels = data_resampled.label; %
allPositions = lay.pos; % get the positions of all channels
misschans=find(~ismember(allLabels , data.label)); %indices of  the channels that are not in the data
nmiss=length(misschans); %the amount of missing channels
cleanchans=find(ismember(allLabels,  data.label)); %indices of labels of the clean channels
nclen=size(allLabels,1) -nmiss; %the amount of clean channels


for t=1:length(data.trial) %for every trial
    tmp0=data.trial{t};
    tmp1=tmp0(1:end,:);  % we write that trial into tmp1
    tmp2=[tmp1;zeros(nmiss,length(tmp1))]; % then we add zeros for every channel that is missing at the end
    tmp3=[tmp2,[cleanchans;misschans]]; % then we write the indices of the channels in the last column...
    tmp4=sortrows(tmp3,length(tmp3)); %... and sort rows along that column
    data.trial{t}=tmp4(:,1:length(tmp4)-1); % finally we overwrite the trial
end

data.label = allLabels; % and we overwrite the labels

cfg=[];
cfg.method = 'triangulation';
%alternatively you can specify cfg.elec
cfg.elecfile = 'standard_1005.elc'; %note if you get a warning: Your data and configuration allow for multiple sensor definitions, you have .elec in the data so you don't need to specify cfg.elecfile
cfg.channel = 'all';
cfg.feedback  = 'yes';
neighbours = ft_prepare_neighbours(cfg, data);


cfg=[];
cfg.badchannel = allLabels(misschans);
cfg.elecfile = 'standard_1005.elc'; %note if you get a warning: Your data and configuration allow for multiple sensor definitions, you have .elec in the data so you don't need to specify cfg.elecfile
cfg.neighbours = neighbours;
data_interp = ft_channelrepair(cfg, data);


%% Re-reference to average EEG channels
% We will re-reference the data to the average of all EEG channels.

cfg=[];
cfg.reref         = 'yes' %(default = 'no')
cfg.refchannel    = data_interp.label(1:128);%cell-array with new EEG reference channel(s), this can be 'all' for a common average reference
cfg.refmethod     = 'avg'% 'median', or 'bipolar' for bipolar derivation of sequential channels (default = 'avg')

data_reref = ft_preprocessing(cfg,data_interp);

%% Filtering

% We will apply a band-stop filter to remove the 50 Hz noise and its harmonics.
% We will also apply a high-pass filter to remove slow drifts in the data.
% Finally, we will apply a low-pass filter to remove high-frequency noise.

% Keep in mind that, depending of your analyses of interest, it is likely that you 
% will need to change the filter settings. For example, if you are interested in
% high-frequency activity, you might want to use a higher low-pass filter.

cfg=[];
cfg.channel= 'all';
cfg.continuous = 'yes';

% Band-stop filter
cfg.bsfilter = 'yes';
cfg.bsfreq = [48 52];
% High-pass filter
cfg.hpfilter = 'yes';
cfg.hpfreq = 0.1;
cfg.hpfilttype ='firws';
% Low-pass filter
cfg.lpfilter = 'yes';
cfg.lpfreq = 30;

data_filter = ft_preprocessing(cfg,data_reref); %preprocess the data

%% Visual inspection (cont. data)
% Now, we are gonna visually inspect the data and mark the artifacts that we want to remove.
% In particular, muscular artifacts and other artifacts that are not related to the brain activity.

% Since we will apply later the ICA for blinks and eye movements, we will not remove these artifacts.
% However, if there a clear sequence of blinks or eye movements together that cannot be isolated, 
% we can remove them.

cfg=[];
cfg.viewmode = 'vertical';
artfx = ft_databrowser(cfg,data_filter) % Here we will have a timestamp of the artifacts.

% Having a timestamp of the artifacts, we can now mark them as "nan" in the data.
% Also, we can save this information for further analysis, so we do not need to do it again.

% IMPORTANT: Keep in mind that any downsampling of the data will change the timestamps of the artifacts.
% We will take care of this later.

%% Marking artifacts selected as "nan"

% We will mark the artifacts selected as "nan" in the data.

cfg=[];
artf2.artfctdef.reject='nan';
data_vis_inspected= ft_rejectartifact(artf2, data_filter);


%% Saving the data as "visual inspected"

% I would save now the data as "visual inspected" to avoid doing this again.
fold_tosave = strcat('xxxx') %WARNING Change to name of the partipants
 save (fold_tosave,'artf2','data_vis_inspected','-v7.3')

 

 
%% PART II. Defining trials and redefine trials

% This part will allow us to load the data after visual inspection and define the trials.
% So, in case we need to change the epochs or the trials, we can do it without repeating the previous steps.

%% Defining trials

%IMPORTANT: Check and recheck files names!

% Define trials use the original signal (only way)
% We will define the trials based on the triggers in the data.
% We will cut the data from 2 seconds before the trigger to 2 seconds after the trigger.
% We will use a function that we have created to define the trials.

cfg=[];
cfg.dataset= 'EEG_dataset/raw_data/p_32_nback.eeg'; % This is the original data, not visual inspected
cfg.headerfile   = 'EEG_dataset/raw_data/p_32_nback.vhdr' % This is the original data, not visual inspected
cfg.fsample = 1000; %sampling rate of the original data
cfg.trigval = [{'S  2'}]; %Object onset
cfg.trialfun = 'MakeTrials_N_back_MOD';
cfg.trialdef.pre= 2; %cut prestimulus in seconds
cfg.trialdef.post= 2; %cut poststimulus in seconds

original_data_epoch=ft_definetrial(cfg);

%% Redifine the trial in the visual inspected dataset (not kicking out trials)

% Now, we will redefine the trials in the visual inspected data.
% Important: We will redefine the trials based on the timestamps of the original data.
% This is important because the timestamps of the artifacts are based on the original data. 

original_samplingrate=1000; %in hz
actual_data_samplingrate=250; %in hz

%First, changing timepoints due to a downsampling
neotrial=round(original_data_epoch.trl/(original_samplingrate/actual_data_samplingrate));

% Redefine trials
cfg=[];
cfg.trl=neotrial%neotrial
data_epoch=ft_redefinetrial(cfg,data_vis_inspected)


%% Checking trials with a visual artifact and select only clean trials
% Since we have marked the artifacts in the data, we can now check the trials that have artifacts and select only the clean trials.
% Data with artifacts will be removed.
% In our case, artifacts are marked as "nan".

counter=0; %counter for bad trials
bad_trial=[]; %bad trials
alltrials=(1:length(data_epoch.trial)); %all trials

% Check for bad trials
for triales=1:length(data_epoch.trial);
ttt=data_epoch.trial{triales};

nodata=isnan(ttt(1,:)); %check if there is any nan in the trial

if sum(nodata)>0
    counter=counter+1;
    bad_trial(counter)=triales;
end

end

cleantrial=setdiff(alltrials,bad_trial); 

% Select clean

cfg=[];
cfg.trials=cleantrial;
data_goodtrials=ft_selectdata(cfg,data_epoch)


%% ICA 
% Now, we will apply the ICA to remove eye movements and blinks.
% We should write down the components that are related to eye movements and blinks.
% Remember, if you don't know or you are not sure, don't remove them!

cfg        = [];
cfg.method = 'runica';
cfg.numcomponent = 20; %number of components to be calculated (default = all). I used this "low" number just for the example
cfg.channel = data_goodtrials.label(1:129);
data_ica= ft_componentanalysis(cfg,data_goodtrials);

data_comp = data_ica; %save the data with ICA components

%% Ensure the ICA components are ordered correctly and set dimord fields

layout_file = 'standard_1005.elc';
elec = ft_read_sens(layout_file);

% Ensure the ICA components are ordered correctly
if size(data_comp.topo, 2) == 1  % If topo is 129x1, reshape it
    data_comp.topo = reshape(data_comp.topo, [size(data_comp.label, 1), size(data_comp.unmixing, 1)]);
end

% Fix unmixing matrix format if needed
if size(data_comp.unmixing, 1) == 1
    data_comp.unmixing = data_comp.unmixing';  % Ensure [components × channels]
end

% Ensure `dimord` is set correctly
data_comp.dimord = 'chan_comp';
data_comp.topodimord = 'chan_comp';
data_comp.unmixingdimord = 'chan_chan';

% Restore correct electrode information
data_comp.elec = ft_datatype_sens(data_comp.elec);

% Ensure layout matches ICA components
cfg = [];
cfg.layout = lay; % Use appropriate layout
cfg.elec = elec;
layout = ft_prepare_layout(cfg, data_comp);

%% Plotting the first components to remove eye movs. and blinks.
% This way plotting only works when using a version of Fieldtrip that is not the latest one.
% In this case, we are using the version 20190224.

layout_file = 'standard_1005.elc';
elec = ft_read_sens(layout_file);

figure;
ax_ = nan(10,1);
ind_ = 0;
for comp = 1:10
    ind_ = ind_+1;
    pspctrm = 0;
    for tr = 1 : numel(data_comp.trial)
        signal  = data_comp.trial{1,tr}(comp,:);
        N       = length(signal);
        nyquist = data_comp.fsample/2;
        fourierCoefs = zeros(size(signal));
        frequencies = linspace(0,nyquist,floor(N/2)+1);
        fourierCoefsF = fft(signal) / N;
        pspctrm = pspctrm + abs(fourierCoefsF(1:length(frequencies)))*2;
    end
    pspctrm = pspctrm/numel(data_comp.trial);
    subplot(5,4,2*(ind_-1)+1),hold on
    
    cfg           = [];
    cfg.component = [comp:comp];       % specify the component(s) that should be plotted
    cfg.elec    = elec; % specify the layout file that should be used for plotting
    cfg.comment   = 'no';
    cfg.layout = lay;
    ft_topoplotIC(cfg, data_comp)
    ax_(ind_)=subplot(5,4,2*ind_);
    set(ax_(ind_),'xlim', [0 60]); %store the handles to the axes of all subplots in an array
    hold on;
    plot(frequencies,pspctrm)
    
end
hax = @(src, ax_) set(ax_, 'xlim', [0 get( src, 'Value' )] );
uicontrol('Style', 'slider', 'Units','normalized', 'Min',20,'Max',300,'Value',60, 'Position', [0 0 1 0.05], 'Callback', @(src,evt) hax( src, ax_ ) ); % allows for zooming on the x axis
%% Reject components
% After visual inspection, we will reject the components that are related to 
% eye movements and blinks.
% Since we wrote them down, we can now reject them.

cfg = [];
cfg.component = []; % an array with the bad components
data = ft_rejectcomponent(cfg, data_comp);
data.reject=cfg.component; %save info of rejected components

data_ica_rejected = data;


%% Recheck signal before saving
% Let see how the signal looks like after removing the components.
% We should expect a cleaner signal, without eye movements and blinks.
% However, ICA is not perfect and some artifacts might still be present.

% Data browser
cfg = [];
cfg.viewmode = 'vertical';
cfg.channel = 'all'
cfg.allowoverlap = 'yes';
ft_databrowser(cfg, data_ica_rejected);


%% I also like to see the data in a multiplot and explore it, to make sure that everything looks fine

cfg = [];
cfg.xlim = [-0.1 0.6]; %time window
cfg.layout=lay; %layout
figure;
ft_multiplotER(cfg,data_ica_rejected)

%% Done. Save it.

% Save the data
fold_tosave = strcat('xxxx') %WARNING Change to name of the partipants
save (fold_tosave,'data_ica_rejected','-v7.3')



