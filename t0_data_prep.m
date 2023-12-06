% Script for wave segmentation tutorial
%
% https://jp.mathworks.com/help/signal/ug/waveform-segmentation-using-deep-learning.html?lang=en

load('mat/wave_table.mat');

% Initialize settings
layers = {'CA1'; 'Interlayer'; 'DG'};
waveforms = {'w1'; 'w2'; 'w3';'bf'};

timeseries = [];
label = [];
index = [];

idx_start = [];
idx_end = [];
curr_idx = 0;

i = 1;

% time_labels = strings(50000,0);
time_labels = {};

% For each waveform type
for w = 1:length(waveforms)
  wf = waveforms{w};
  
  % For each layer type
  for l = 1:length(layers)
    ly = layers{l};
    ts = wave.(ly).(strcat('shf_',wf));
    ts = cellfun(@transpose,ts,'UniformOutput',false);
%     ts_ = vertcat(ts{:});
%     timeseries = [timeseries; ts_];
    timeseries = vertcat(timeseries, ts);
%     timeseries = vertcat(timeseries, ts);

    ln = cellfun('length', ts);
    
    [idx_start, idx_end] = make_index(ln);
    idx_start = idx_start + curr_idx;
    idx_end = idx_end + curr_idx;

    idx = [idx_start, idx_end];
    index = [index; idx];

    label = [label; repmat(wf, length(ln), 1)];    

    label_idx = [idx_start(1):idx_end(end)];
%     time_labels(label_idx) = wf;
    time_labels = vertcat(time_labels, repmat({wf}, length(ts),1));

    curr_idx = idx_end(end);
  end
end

% time_labels(1) = [];

wave_index_labels = table(index, string(label));
data = {timeseries; wave_index_labels};

% M = signalMask(data{2});
% plotsigroi(M,data{1});

% proj_path = fileparts(mfilename('fullpath'));
% datasetFolder = fullfile(proj_path, 'ML_dataset');
% fname = 'Mouse_191127';
% save_data(timeseries, wave_index_labels, fps, fname, datasetFolder);

%%

% Create start and end index based on segment lengths
function [idx_start, idx_end] = make_index(ln)
n = length(ln);
idx_start = ones(5,1);
idx_end = ln;

for i = 2:n
%     idx_end(i) = sum(idx_end(i-1:i))+1;
    idx_end(i) = sum(idx_end(i-1:i));
end

idx_start(2:n) = idx_end(1:n-1)+1;


end

%%

% Save datastore

function dat = save_data(signal, labels, Fs, fname, datasetFolder)
timeseries = signal;
wave_index_labels = labels;

save(fullfile(datasetFolder,fname), 'Fs', 'timeseries', 'wave_index_labels')

end
