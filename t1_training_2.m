
rng default %Default chooses same indices for repeatability
% [trainIdx,validIdx,testIdx] = dividerand(length(timeseries),0.60,0.10,0.30);
[trainIdx,~,testIdx] = dividerand(length(timeseries),0.60,0,0.30);

%% Training data
XTrain = timeseries(trainIdx);
YTrain = categorical(time_labels(trainIdx));

for k = 1:length(XTrain)
XTrain{k,1} = vertcat(XTrain{k,1},fsst(XTrain{k,1},10,kaiser(100,10)));
end

DataParts = cell(size(XTrain,1), size(XTrain,2), 2);

for k = 1:length(XTrain)
DataParts{k,:,1} = real(XTrain{k,1});
%DataParts{k,:,2} = imag(XTrain{k,1});

XTrain{k,:} = DataParts{k,:,1};
XTrain{k,:} = vertcat(XTrain{k,:},DataParts{k,:,2});
end


 %% Plot sample fsst spectogram
% w1_idx = 1:6;   % wf 1, layer 1
% w2_idx = 19:24; % wf 2, layer 1; with one extra
% w3_idx = 34:39; % wf 3, layer 1; with one extra
% wall_idx = [w1_idx, w2_idx, w3_idx];
% 
% figure;
% n = length(wall_idx);
% for p = 1:n
%     subplot(3,6,p)
%     i = wall_idx(p);
% 
%     if (i==24 || i==39)
%         continue
%     end
% 
%     fsst(timeseries{i},10,kaiser(100,10),'yaxis');
%     ylim([0 0.05])
%     colorbar
%     caxis([0 250])
% end
% 
% 
% 
% 
% 
% subplot(1,n,k)
% fsst(timeseries{29},10,kaiser(100,10),'yaxis');
% ylim([0 0.05])
% colorbar
% caxis([0 1000])
% subplot(1,n,k)
% fsst(timeseries{46},10,kaiser(100,10),'yaxis');
% ylim([0 0.05])
% colorbar
% caxis([0 1000])
% 

% 
% figure;
% subplot(2,1,1)
% plot(wave.CA1.shf(1:2500))
% subplot(2,1,2)
% plot(Y_Pred_all(1:2500))
% linkaxes

%% Validation data

% XValid = timeseries(validIdx);
% YValid = categorical(time_labels(validIdx));
% 
% for k = 1:length(XValid)
% XValid{k,1} = vertcat(XValid{k,1},fsst(XValid{k,1},10,kaiser(100,10)));
% end
% 
% DataParts = cell(size(XValid,1), size(XValid,2), 2);
% 
% for k = 1:length(XValid)
% DataParts{k,:,1} = real(XValid{k,1});
% %DataParts{k,:,2} = imag(XTrain{k,1});
% 
% XValid{k,:} = DataParts{k,:,1};
% XValid{k,:} = vertcat(XValid{k,:},DataParts{k,:,2});
% end

%% Testing data

XTest = timeseries(testIdx);
YTest = categorical(time_labels(testIdx));

for k = 1:length(XTest)
XTest{k,1} = vertcat(XTest{k,1},fsst(XTest{k,1},10,kaiser(100,10)));
end

DataParts = cell(size(XTest,1), size(XTest,2), 2);

for k = 1:length(XTest)
DataParts{k,:,1} = real(XTest{k,1});
%DataParts{k,:,2} = imag(XTrain{k,1});
end

for k = 1:length(XTest)
XTest{k,:} = DataParts{k,:,1};
XTest{k,:} = vertcat(XTest{k,:},DataParts{k,:,2});
end

%% Pre-visualiztion and sorting

figure
plot(XTrain{1}')
xlabel("Time Step")
title("Training Observation 1")
numFeatures = size(XTrain{1},1);
legend("Feature " + string(1:numFeatures),'Location','northeastoutside')

numObservations = numel(XTrain);
for i=1:numObservations
    sequence = XTrain{i};
    sequenceLengths(i) = size(sequence,2);
end

[sequenceLengths,idx] = sort(sequenceLengths);
XTrain = XTrain(idx);
YTrain = YTrain(idx);

figure
bar(sequenceLengths)
% ylim([0 30])
xlabel("Sequence")
ylabel("Length")
title("Sorted Data")

% for k = 1:length(XTrain)
% XTrain{k,1} = vertcat(XTrain{k,1},fsst(XTrain{k,1},10));
% end

%% Set-up neural net

numClasses = 4; %Number of label outputs
miniBatchSize = 13; %Divide dataset to 13 batches so that lengths are equal
inputSize = numFeatures; %52 features
numHiddenUnits = 187; %Minimum number of timesteps?
maxEpochs = 30;

layers = [ ...
    sequenceInputLayer(numFeatures,'Normalization','zscore') 
    bilstmLayer(numHiddenUnits,'OutputMode','last')
    dropoutLayer(0.2)
    fullyConnectedLayer(numClasses)
    softmaxLayer
    classificationLayer]

options = trainingOptions('adam', ...
    'ExecutionEnvironment','cpu', ...
    'GradientThreshold',1, ...
    'MaxEpochs',maxEpochs, ...
    'MiniBatchSize',miniBatchSize, ...
    'SequenceLength','longest', ...
    'Shuffle','never', ...
    'Verbose',0, ...
    'Plots','training-progress')
    %'ValidationData',{XValid,YValid});

% tallTrainSet = tall(trainDs);
% tallTestSet = tall(testDs);
% trainData = gather(tallTrainSet);
% testData = gather(tallTestSet);

% fsstTrainDs = transform(sds,@(x)extractFSSTFeatures(x,10));
% fsstTallTrainSet = tall(fsstTrainDs);
% fsstTrainData = gather(fsstTallTrainSet);
return
net = trainNetwork(XTrain,YTrain,layers,options);

%% Check with testing data

% predTest = classify(rawNet,testData(:,1),'MiniBatchSize',50);
YPred = classify(net,XTest,'MiniBatchSize',50);

checklist = [YTest, YPred]
confusionchart([YTest],[YPred],'RowSummary','row-normalized','ColumnSummary','column-normalized','Normalization','absolute');
