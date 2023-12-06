
%% Access calcium imaging data

filename = '/Users/joshwoah/Documents/NAIST/Graduate Research/Results/Kainic Acid Results/Normalize of Cleaned All KA Mice.csv'; %Changed file to processed data from GraphPad 
mytable = readtable(filename, 'Delimiter', ',','ReadRowNames', true, 'ReadVariableNames', true,'VariableNamingRule', 'preserve');
mynewtable = rows2vars(mytable, 'VariableNamingRule', 'preserve');
layer_assignments = mynewtable{:, 1};
my_array = table2array(mynewtable(:,2:end));

%Remove NaN values (does not work for beginning values)

for i = 1:size(my_array, 1)
    row = my_array(i,:);
    idx = ~isnan(row);

    xq = cumsum(idx-diff([1,idx])/2); %makes a vector where the index distance/difference between two NaN values is 1 which will be used as
    row_interp = interp1(1:nnz(idx),row(idx),xq);
    new_array(i,:) = row_interp;
  
end

%Replace beginning NaN values with median

for i = 1:size(new_array, 1)
    row = new_array(i, :);
    
    % Calculate median of the current row
    rowMedian = median(row, 'omitnan');
    
    % Replace NaN values with row median
    row(isnan(row)) = rowMedian;
    
    % Update the matrix with the modified row
    new_array(i, :) = row;
end

%% Access behavioral data

filename = '/Users/joshwoah/Documents/NAIST/Graduate Research/Results/Kainic Acid Results/Racine scale behavior.csv';
mytableb = readtable(filename, 'Delimiter', ',','ReadRowNames', true, 'ReadVariableNames', true,'VariableNamingRule', 'preserve');
mynewtableb = rows2vars(mytableb, 'VariableNamingRule', 'preserve');
my_arrayb = table2array(mynewtableb(:,2:end));

% Repeat each state 30 times to match time resolution of imaging data

% Initialize the output matrix and padding
new_arrayb = zeros(size(my_arrayb, 1), size(my_arrayb, 2) * 30 + 900);

% Iterate over each row of the input matrix
for i = 1:size(my_arrayb, 1)
    % Repeat each element in the row 30 times
    row = repelem(my_arrayb(i, :), 30);
    % Pad 900 zeros at the beginning of the row (Stage 0 before KA)
    row = [zeros(1, 900), row];
    % Append the modified row to the output matrix
    new_arrayb(i, :) = row;
end

%Add 1 to each behavioral state, where 1 is Stage 0 and 6 is Stage 5 (coz HMM doesn't work with states numbered as 0)
new_arrayb = new_arrayb + 1;

%Repeat rows per number of ROIs per layer

n_mouse1_CA1 = 4;
n_mouse1_IL = 7;
n_mouse1_DG = 4;

n_mouse3_CA1 = 4;
n_mouse3_IL = 4;
n_mouse3_DG = 10;

n_mouse4_CA1 = 6;
n_mouse4_IL = 5;
n_mouse4_DG = 1;

n_mouse5_CA1 = 7;
n_mouse5_IL = 8;
n_mouse5_DG = 5;

n_mouse6_CA1 = 7;
n_mouse6_IL = 6;
n_mouse6_DG = 9;

n_mouse7_CA1 = 7;
n_mouse7_IL = 13;
n_mouse7_DG = 9;

CA1repcounts = [n_mouse1_CA1, n_mouse3_CA1, n_mouse4_CA1, n_mouse5_CA1, n_mouse6_CA1, n_mouse7_CA1];
ILrepcounts = [n_mouse1_IL,n_mouse3_IL, n_mouse4_IL, n_mouse5_IL, n_mouse6_IL, n_mouse7_IL];
DGrepcounts = [n_mouse1_DG,n_mouse3_DG, n_mouse4_DG, n_mouse5_DG, n_mouse6_DG, n_mouse7_DG];

% CA1repcounts = [n_mouse1_CA1, n_mouse3_CA1, n_mouse4_CA1];
% ILrepcounts = [n_mouse1_IL,n_mouse3_IL, n_mouse4_IL];
% DGrepcounts = [n_mouse1_DG,n_mouse3_DG, n_mouse4_DG];


CA1_behav = repelem(new_arrayb, CA1repcounts, 1);
IL_behav = repelem(new_arrayb, ILrepcounts, 1);
DG_behav = repelem(new_arrayb, DGrepcounts, 1);

final_behav_array = [CA1_behav;IL_behav;DG_behav];

%% Initialize values - replace with actual behavioral data later

% seq is my calcium imaging data
% states is my Racine behavioral scale (behavioral state)

test_seq = new_array;
test_states = final_behav_array; 

% Convert calcium data and behavioral states to sequence and state cells
seq_cell = mat2cell(test_seq, ones(1,size(test_seq,1)), size(test_seq,2));
seq_cell = seq_cell';
states_cell = mat2cell(test_states, ones(1,size(test_states,1)), size(test_states,2));
states_cell = states_cell';

% Combine all sequences and states into one array
all_seq = [seq_cell{:}];
all_states = [states_cell{:}];
uniqueSymbols = unique(all_seq);
uniqueStates = unique(all_states);

% %% Determine transition and emission probability matrices (below is from MATLAB example)
% 
% TRANS = [.9 .1; .05 .95];
% 
% EMIS = [1/6, 1/6, 1/6, 1/6, 1/6, 1/6;...
% 7/12, 1/12, 1/12, 1/12, 1/12, 1/12];
% 
% %Generate random emission and states
% 
% [seq,states] = hmmgenerate(1000,TRANS,EMIS); %'seq' and  'states' are a random sequence of emissions (# from dice or CA IMG data) and states (red or green or behavior), respectively
% % output: # of states based on rows of EMIS (2) and emissions based on columns (6)
% 


%% Estimate the transition and emission matrices

%Using hmmestimate - requires that you know the sequence of states that the model went through to generate seq

[TRANS_EST, EMIS_EST] = hmmestimate(all_seq,all_states,'Symbols',uniqueSymbols,'Statenames',uniqueStates);

%Using hmmtrain if you do not know the sequence of states, but you have initial guesses for TRANS and EMIS
% I will use TRANS_EST and EMIS_EST from hmmestimate as initial guesses and feed it to hmmtrain

%% Split the train, validation, and test data

% Set the random seed for reproducibility (optional)
rng(42);


% Define the split ratios
trainRatio = 0.7;
testRatio = 0.3;


% Compute the number of samples for each set
numSamples = numel(seq_cell);
numTrain = floor(trainRatio * numSamples);
numTest = numSamples - numTrain;


% Create a random partition for train and test
partitionTrainTest = cvpartition(numSamples, 'HoldOut', testRatio);

% Training data
trainIndicesAll = training(partitionTrainTest);
trainDataSeqAll = seq_cell(trainIndicesAll);
trainDataStateAll = states_cell(trainIndicesAll);
trainLayerAssignAll = layer_assignments(trainIndicesAll);

% Test data
testIndices = test(partitionTrainTest);
testDataSeq = seq_cell(testIndices);
testDataState = states_cell(testIndices);

% Set the number of folds (K) for cross-validation
K = 5;

% Create a random stratified partition for cross-validation
partitionTrainVal = cvpartition(numTrain, 'KFold', K);

% Initialize variables to store evaluation metrics
logLikelihoods = zeros(K, 1);
predictions = cell(K, 1);
trans_mat = cell(K,1);
emis_mat = cell(K,1);
accuracy_mat = zeros(K, 1);

% Start the timer
tic;

% Perform K-fold cross-validation
for k = 1:K
    % Get the indices of the training and validation sets for the current fold
    trainIndices = partitionTrainVal.training(k);
    valIndices = partitionTrainVal.test(k);

    % Split the data into training and validation sets
    trainDataSeq = trainDataSeqAll(trainIndices);
    valDataSeq = trainDataSeqAll(valIndices);
    trainDataState = trainDataStateAll(trainIndices);
    valDataState = trainDataStateAll(valIndices);

    % Train the HMM using hmmtrain with the training data
    [TRANS_EST2, EMIS_EST2] = hmmtrain(trainDataSeq, TRANS_EST, EMIS_EST,'Symbols',uniqueSymbols,'Verbose',true,'Tolerance',1e-6,'Maxiterations',200);


    % Evaluate the HMM on the validation data
    numSequences = numel(valDataSeq); % Number of sequences
    decodedStates = cell(1,numSequences); % Cell array to store decoded states
    decodedlogP = zeros(1,numSequences); % Cell array to store log likelihood
    decodedAccuracy = zeros(1,numSequences);% Cell array to store decoding accuracies (across the cell array since hmmviterbi doesn't take cells)

    % Iterate over each sequence in the cell array
    for i = 1:numSequences
        sequence = valDataSeq{i}; % Get the sequence from the cell array

        % Add pseudo-transition
        minNonZeroValue = min(TRANS_EST2(TRANS_EST2 > 0), [], 'omitnan');
        pseudoTransition = minNonZeroValue;  % Small non-zero probability for pseudo-transitions
        TRANS_EST2(TRANS_EST2 == 0) = pseudoTransition;

        % Add pseudo-emissions
        epsilon = 1e-6;  % Small non-zero value for smoothing
        zeroIndices = (EMIS_EST2 == 0);
        EMIS_EST2(zeroIndices) = epsilon;

        % Normalize the emission probabilities
        EMIS_EST2 = EMIS_EST2 ./ sum(EMIS_EST2, 2);


        % Perform decoding using hmmviterbi
        [likelystates, logP] = hmmviterbi(sequence, TRANS_EST2, EMIS_EST2, 'Symbols', uniqueSymbols, 'Statenames', uniqueStates);
        decoding_accuracy = sum(valDataState{i}==likelystates)/numel(valDataState{i});

        % Store the decoded sequence in the output cell array
        decodedStates{i} = likelystates;
        decodedlogP(i) = logP;
        decodedAccuracy(i) = decoding_accuracy;
        
    end

   fold_likelihood = mean(decodedlogP);
   fold_accuracy = mean(decodedAccuracy);
    disp(['Iteration ', num2str(k), ': Fold Accuracy = ', num2str(fold_accuracy)]);


    % Store the results for this fold
    logLikelihoods(k) = fold_likelihood;
    predictions{k} = decodedStates;
    trans_mat{k} = TRANS_EST2;
    emis_mat{k} = EMIS_EST2;
    accuracy_mat(k) = fold_accuracy;


    % Update the initialModel with the trained Model for the next iteration
    TRANS_EST = TRANS_EST2;
    EMIS_EST = EMIS_EST2;
end

% Compute the average log likelihood across all folds
averageLogLikelihood = mean(logLikelihoods);

% Select the best model based on the highest average log likelihood
bestModelIndex = find(logLikelihoods == max(logLikelihoods));
[TRANS_EST_BEST, EMIS_EST_BEST] = hmmtrain(trainDataSeqAll, trans_mat{bestModelIndex},emis_mat{bestModelIndex},'Symbols',uniqueSymbols,'Verbose',true,'Tolerance',1e-6,'Maxiterations',200);  % Train the best model using all available data

% % Save the matrices
filename = 'hmm_results.mat';
save(filename, 'TRANS_EST_BEST','EMIS_EST_BEST');


% Stop the timer
elapsedTime = toc;

% Display the elapsed time
fprintf('Elapsed time: %.2f seconds\n', elapsedTime);

return

% %% Visualize
% 
% state_names = uniqueStates;
% figure;
% imagesc(TRANS_EST);
% caxis([0, 0.1]); 
% colorbar;
% colormap('parula');
% xlabel('To State');
% ylabel('From State');
% set(gca, 'XTick', 1:numel(state_names), 'XTickLabel', state_names, ...
%     'YTick', 1:numel(state_names), 'YTickLabel', state_names);
% title('Transition Probability Matrix');



%% Decoding - getting the likely states based on seq

% Evaluate the HMM on the test data
numSequences = numel(testDataSeq); % Number of sequences
test_decodedStates = cell(1,numSequences); % Cell array to store decoded states
test_decodedlogP = zeros(1,numSequences); % Cell array to store log likelihood
test_decodedAccuracy = zeros(1,numSequences);% Cell array to store decoding accuracies (across the cell array since hmmviterbi doesn't take cells)

for i = 1:numel(testDataSeq)

    % Add pseudo-transition
    minNonZeroValue = min(TRANS_EST_BEST(TRANS_EST_BEST > 0), [], 'omitnan');
    pseudoTransition = minNonZeroValue;  % Small non-zero probability for pseudo-transitions
    TRANS_EST_BEST(TRANS_EST_BEST == 0) = pseudoTransition;

    % Add pseudo-emissions
    epsilon = 1e-6;  % Small non-zero value for smoothing
    zeroIndices = (EMIS_EST_BEST == 0);
    EMIS_EST_BEST(zeroIndices) = epsilon;

    % Normalize the emission probabilities
    EMIS_EST_BEST = EMIS_EST_BEST ./ sum(EMIS_EST_BEST, 2);

    %
     [likelystates, logP] = hmmviterbi(testDataSeq{i}, TRANS_EST_BEST, EMIS_EST_BEST, 'Symbols', uniqueSymbols, 'Statenames', uniqueStates);

    % Test accuracy of hmmviterbi (given the actual states)
    test_decoding_accuracy = sum(testDataState{i}==likelystates)/numel(testDataState{i});

    % Store the decoded sequence in the output cell array
    test_decodedStates{i} = likelystates;
    test_decodedlogP(i) = logP;
    test_decodedAccuracy(i) = test_decoding_accuracy;

end

% bestTestIndex = find(test_decodedlogP == max(test_decodedlogP));
bestTestIndex = find(test_decodedAccuracy == max(test_decodedAccuracy));
test_fold_likelihood = mean(test_decodedlogP);
test_fold_accuracy = mean(test_decodedAccuracy);
max_fold_accuracy = max(test_decodedAccuracy);
 disp(['Mean Accuracy = ', num2str(test_fold_accuracy)]);
 disp(['Max Accuracy = ', num2str(max_fold_accuracy)]);
bestLogLikelihood = test_decodedlogP(bestTestIndex);


 %% Visualize

state_names = uniqueStates;
figure;
imagesc(TRANS_EST_BEST);
caxis([0, 0.1]); 
colorbar;
colormap('parula');
xlabel('To State');
ylabel('From State');
set(gca, 'XTick', 1:numel(state_names), 'XTickLabel', state_names, ...
    'YTick', 1:numel(state_names), 'YTickLabel', state_names);
title('Transition Probability Matrix');
saveas(gcf, 'transition_matrix.png');

%% Sample visualization
time_series = testDataSeq{bestTestIndex};  % Your continuous time series cell array
classification = testDataState{bestTestIndex};  % Your categorical classification cell array

% Define color map
color_map = parula(numel(uniqueStates));  % RGB values for each category

% Define the labels
categories = {'Stage 0', 'Stage 0.5', 'Stage 1', 'Stage 1.5', 'Stage 2', 'Stage 2.5', 'Stage 3', 'Stage 3.5', 'Stage 4'};

% Map classification values to integers
classification_int = zeros(size(classification));
for i = 1:numel(uniqueStates)
    classification_int(classification == uniqueStates(i)) = i;
end

% Plot the time series data
fig_actual = figure('Position', get(0, 'ScreenSize'));
hold on;
for i = 1:numel(time_series)
    x = i;
    y = time_series(i);
    c = classification_int(i);
    scatter(x, y, 10, c, 'filled');
end
hold off;

% Add color bar with color-coded labels
colormap(color_map);
c = colorbar('Ticks', 1:numel(categories), 'TickLabels', categories, 'Location', 'eastoutside');
c.Label.String = 'Racine Stage';
c.Label.FontSize = 14;
c.FontSize = 12;

% Add axis labels and title
ylabel('Fluorescence Intensity', 'FontSize', 16);
xlabel('Frames', 'FontSize', 16);
xlim([1000, max(x)]);
title('Calcium Imaging Data and Racine Behavior', 'FontSize', 18);

% Adjust font size and style of axis ticks
set(gca, 'FontSize', 15);  % Adjust the font size here

% Adjust font size and style of tick labels
set(gca, 'FontName', 'Arial');  % Adjust the font style here

% Adjust figure size
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);



%% Discretize the test sequence for hmmdecode

numBins = 81270; % Specify the number of bins

% Compute the histogram and assign data points to bins
[counts, edges] = histcounts(testDataSeq{bestTestIndex}, numBins);

% Map the data points to discrete symbols based on the bin assignments
discreteSeq = discretize(testDataSeq{bestTestIndex}, edges);

% Assuming you have the sequence stored in the variable 'seq'

% Calculate the histogram
binEdges = min(discreteSeq):1:max(discreteSeq);
binCounts = histcounts(discreteSeq, binEdges);

% Plot the histogram
figure;
bar(binEdges(1:end-1), binCounts);
xlabel('States');
ylabel('Frequency');
title('State Histogram');



%% Decoding - estimating the posterior probability of seq given the current state

% The posterior state probabilities are the probabilities that the model is in a particular state when it generates a symbol in seq, given that seq is emitted. 
% You compute the posterior state probabilities with hmmdecode


[PSTATES,logpseq] = hmmdecode(discreteSeq,TRANS_EST_BEST,EMIS_EST_BEST)

% Display the state probabilities and log likelihood
disp('State Probabilities:');
disp(PSTATES);
disp('Log Likelihood:');
disp(logpseq);

% Output: PSTATES is an M-by-L matrix, where M is the number of states and L is the length of seq
% PSTATES(i,j) is the probability that the model is in state i when it generates the jth symbol of seq
% Example: hmmdecode begins with the model in state 1 at step 0, prior to the first emission. PSTATES(i,1) is the probability that the model is in state i at the following step 1
%logpseq will get the log of the probability because the probability of a sequence tends to 0 as the length of the sequence increases, and the probability of a sufficiently long sequence becomes less than the smallest positive number your computer can represent. hmmdecode returns the logarithm of the probability to avoid this problem.


%% Visualize PSTATES

% Calculate the minimum and maximum values from the matrix
min_value = min(PSTATES(:));
max_value = max(PSTATES(:));

% Define the stage names in inverted order
categories = {'Stage 4', 'Stage 3.5', 'Stage 3', 'Stage 2.5', 'Stage 2', 'Stage 1.5', 'Stage 1', 'Stage 0.5', 'Stage 0'};

% Create the grayscale heatmap
fig_pred = figure;
imagesc(flipud(PSTATES));  % Invert the rows using flipud
colormap(gray);
caxis([min_value, max_value]);

% Add colorbar and labels
colorbar;
xlabel('Time');
xlim([1000, size(PSTATES, 2)]);
ylabel('Racine Stage');
title('Predicted seizure behavior');

% Set the y-axis tick labels with inverted stage names
yticks(1:size(PSTATES, 1));
yticklabels(categories);

% Invert the y-axis limits
ylim([0.5, size(PSTATES, 1) + 0.5]);

% Adjust font size and style of axis ticks
set(gca, 'FontSize', 15);  % Adjust the font size here

% Adjust font size and style of tick labels
set(gca, 'FontName', 'Arial');  % Adjust the font style here

%%
% Create a subplot with two rows and one column
subplot(2, 1, 1);
copyobj(allchild(fig_actual), gca);

subplot(2, 1, 2);
copyobj(allchild(fig_pred), gca);


%% Comparing the cross-entropy and accuracy with random data

% Number of unique states (classes)
ran_numStates = 9;

% Number of samples
ran_numSamples = 81370;  % Update with your actual number of samples

% Generate random probabilities for each sample
randomProbabilities = rand(ran_numSamples, ran_numStates);
randomProbabilities = randomProbabilities ./ sum(randomProbabilities, 2); % Normalize to sum to 1

% Choose the true class for each sample (example, replace with your own true classes)
trueClasses = randi(ran_numStates, ran_numSamples, 1);  % Assuming classes are represented as integers

% Extract the true class probabilities from randomProbabilities
prob_true_class = randomProbabilities(sub2ind(size(randomProbabilities), 1:ran_numSamples, trueClasses'));

% Calculate cross-entropy for the random baseline
crossEntropy_baseline = -log(prob_true_class);

fprintf('Cross-Entropy for Random Baseline: %f\n', mean(crossEntropy_baseline));

% Choose the class with the highest random probability as the predicted label
[~, predictedLabels] = max(randomProbabilities, [], 2);

% Generate random true labels (example, replace with your own true labels)
trueLabels = randi(ran_numStates, ran_numSamples, 1);  % Assuming labels are represented as integers

% Calculate accuracy
accuracy = sum(predictedLabels == trueLabels) / ran_numSamples;

fprintf('Accuracy of Random HMM Model: %.2f%%\n', accuracy * 100);

