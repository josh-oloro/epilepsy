# Decoding Epileptic Seizure States from Hippocampal Calcium Imaging Using Machine Learning

This repository contains the MATLAB code accompanying the research paper published in the *Journal of Neural Engineering*:

> **[Brain-implantable needle-type CMOS imaging device enables multi-layer dissection of seizure calcium dynamics in the hippocampus](https://iopscience.iop.org/article/10.1088/1741-2552/ad5c03)**  
> *Journal of Neural Engineering*, 2024. DOI: 10.1088/1741-2552/ad5c03

---

## Overview

Epilepsy is characterized by recurrent, unprovoked seizures. Understanding how neural activity in the hippocampus evolves across seizure stages is critical for developing better diagnostic and therapeutic tools. This project combines two-photon calcium imaging with machine learning to:

1. **Classify hippocampal calcium waveform types** using a Bidirectional Long Short-Term Memory (BiLSTM) deep neural network.
2. **Decode seizure behavioral states** (Racine scale) from hippocampal calcium imaging signals using a Hidden Markov Model (HMM).

Data were collected from the hippocampus of kainic acid (KA)-treated mice — a well-established model of temporal lobe epilepsy — across three anatomical layers:
- **CA1** (Stratum oriens/pyramidale region, OPR)
- **Interlayer / LM** (Lacunosum moleculare, LM)
- **Dentate Gyrus / DG** (Granule/hilus region, GH)

---

## Repository Structure

```
.
├── t0_data_prep.m                          # Data preparation for BiLSTM waveform classification
├── t1_training_2.m                         # BiLSTM model training and hyperparameter setup
├── t2_testing.m                            # Sliding-window inference using the trained BiLSTM
├── hmm.m                                   # HMM training, cross-validation, and decoding
├── Normalize of Cleaned All KA Mice.csv    # Preprocessed calcium imaging data (6 mice)
├── wave_table.mat                          # Waveform segments by type and layer
└── hmm_results.mat                         # Saved best HMM transition/emission matrices
```

---

## Dataset

### Calcium Imaging Data
- **File:** `Normalize of Cleaned All KA Mice.csv`
- **Subjects:** 6 mice (Mouse 1, 3, 4, 5, 6, 7) injected with kainic acid to induce status epilepticus
- **Layers recorded:** CA1 (OPR), Interlayer (LM), Dentate Gyrus (GH)
- **Imaging rate:** 10 Hz (1 frame per 100 ms; behavioral data is repeated at 30× to match)
- **Signal:** Normalized fluorescence intensity (% ΔF/F or equivalent)
- **Total ROIs per layer per mouse:**

| Mouse | CA1 ROIs | Interlayer ROIs | DG ROIs |
|-------|----------|-----------------|---------|
| 1     | 4        | 7               | 4       |
| 3     | 4        | 4               | 10      |
| 4     | 6        | 5               | 1       |
| 5     | 7        | 8               | 5       |
| 6     | 7        | 6               | 9       |
| 7     | 7        | 13              | 9       |

### Behavioral Data (Racine Scale)
Seizure severity was scored using the **modified Racine scale** with 9 states:

| State | Description                          |
|-------|--------------------------------------|
| 0     | No behavioral change (baseline)      |
| 0.5   | Mouth/facial movements               |
| 1     | Head nodding                         |
| 1.5   | Forelimb clonus                      |
| 2     | Bilateral forelimb clonus            |
| 2.5   | Bilateral clonus with rearing        |
| 3     | Rearing and falling                  |
| 3.5   | Continuous rearing and falling       |
| 4     | Tonic-clonic seizure                 |

- Behavioral scoring was done at 1-minute intervals and up-sampled to match the 10 Hz imaging frame rate.
- 900 frames (15 minutes) of baseline activity before KA injection are prepended as Stage 0.

---

## Pipeline 1: Waveform Classification (BiLSTM)

This pipeline classifies individual calcium imaging waveform segments into one of four categories.

### Waveform Types
| Label | Description                  |
|-------|------------------------------|
| `w1`  | Waveform type 1              |
| `w2`  | Waveform type 2              |
| `w3`  | Waveform type 3              |
| `bf`  | Background fluorescence      |

### Scripts

#### `t0_data_prep.m` — Data Preparation
- Loads pre-segmented waveform data from `wave_table.mat`.
- Iterates over all waveform types (`w1`, `w2`, `w3`, `bf`) and all three hippocampal layers.
- Assembles a time-series cell array (`timeseries`) paired with segment labels and index ranges.
- Produces variables `timeseries`, `time_labels`, and `wave_index_labels` for downstream use.

#### `t1_training_2.m` — Feature Extraction and Model Training
- Splits data into training (60%) and test (30%) sets using `dividerand` with a fixed random seed.
- **Feature extraction:** Applies the **Fourier Synchrosqueezed Transform (FSST)** to each waveform segment using a Kaiser window (`kaiser(100, 10)`) at 10 Hz. The real part of the FSST coefficients is stacked with the raw signal to form a multi-channel input.
- **Sequence sorting:** Sequences are sorted by length for efficient mini-batch processing.
- **Model architecture:**

  | Layer                    | Configuration                    |
  |--------------------------|----------------------------------|
  | Sequence Input           | `numFeatures` features, z-score normalization |
  | Bidirectional LSTM       | 187 hidden units, output last state |
  | Dropout                  | 20%                              |
  | Fully Connected          | 4 outputs (one per class)        |
  | Softmax                  | —                                |
  | Classification           | —                                |

- **Training options (Adam optimizer):**
  - Max epochs: 30
  - Mini-batch size: 13
  - Gradient threshold: 1
  - Sequence padding: longest
  - Execution: CPU

#### `t2_testing.m` — Sliding-Window Inference
- Applies the trained BiLSTM model to continuous hippocampal recordings using a **sliding window** approach.
- Window size: 500 frames (~50 seconds at 10 Hz).
- Window step: 1 frame (stride of 1).
- For each window, FSST features are extracted and the model predicts a waveform class.
- Predictions are concatenated and visualized against the raw CA1 signal.

---

## Pipeline 2: Seizure State Decoding (HMM)

This pipeline models the relationship between calcium imaging fluorescence and Racine behavioral states using a **discrete Hidden Markov Model**.

### Script: `hmm.m`

#### Data Loading and Preprocessing
1. **Calcium imaging:** Loads `Normalize of Cleaned All KA Mice.csv`, transposes, interpolates internal NaN values (linear interpolation), and replaces leading NaN values with the row median.
2. **Behavioral data:** Loads `Racine scale behavior.csv`, up-samples each 1-minute behavioral score by repeating it 30 times to match the imaging frame rate, and prepends 900 frames of Stage 0 (baseline).
3. **State encoding:** Behavioral stages (0–4) are shifted by +1 (so Stage 0 → State 1, Stage 4 → State 5) to be compatible with MATLAB's HMM functions.
4. **ROI replication:** Behavioral sequences are replicated per-ROI so each imaging channel has a corresponding state sequence.

#### HMM Training (K-Fold Cross-Validation)
- **Train/test split:** 70% training, 30% test (fixed seed `rng(42)`).
- **Initial parameter estimation:** `hmmestimate` is used to compute initial transition (`TRANS_EST`) and emission (`EMIS_EST`) matrices from all sequences and known states.
- **Cross-validation:** 5-fold cross-validation on the training set.
  - Each fold refines the model using `hmmtrain` (Baum-Welch / EM algorithm, tolerance 1e-6, max 200 iterations).
  - Viterbi decoding (`hmmviterbi`) is applied to each validation sequence.
  - **Pseudo-count smoothing** is applied to prevent zero-probability transitions/emissions.
- **Best model selection:** The fold with the highest log-likelihood is used to retrain on all training data (`hmmtrain` with full training set).
- **Results saved:** Best transition (`TRANS_EST_BEST`) and emission (`EMIS_EST_BEST`) matrices are saved to `hmm_results.mat`.

#### Decoding and Evaluation
- **Test set decoding:** Viterbi algorithm decodes the most likely seizure state sequence for each test ROI.
- **Metrics reported:**
  - Mean decoding accuracy across test sequences
  - Maximum decoding accuracy
  - Log-likelihood of decoded sequences
- **Random baseline comparison:** Cross-entropy and accuracy of a uniformly random classifier are computed for reference.

#### Posterior State Probabilities
- `hmmdecode` computes the posterior probability of being in each Racine stage at every time step, given the best HMM and a discretized test sequence.
- Results are visualized as a grayscale heatmap (time × stage) showing the predicted seizure trajectory.

#### Visualization Outputs
| Figure | Content |
|--------|---------|
| Transition probability matrix | Heatmap of state-to-state transition probabilities (parula colormap) |
| Calcium + Racine overlay | Scatter plot of fluorescence intensity color-coded by Racine stage |
| Posterior state heatmap | Grayscale heatmap of `PSTATES` (probability of each seizure stage over time) |

---

## Requirements

- **MATLAB** R2021a or later (recommended)
- **Signal Processing Toolbox** (for `fsst`)
- **Deep Learning Toolbox** (for `bilstmLayer`, `trainNetwork`, `classify`)
- **Statistics and Machine Learning Toolbox** (for `cvpartition`, `hmmestimate`, `hmmtrain`, `hmmviterbi`, `hmmdecode`)

---

## Usage

### Step 1: Waveform Classification

```matlab
% 1. Prepare data (loads wave_table.mat)
run('t0_data_prep.m')

% 2. Train BiLSTM model
run('t1_training_2.m')

% 3. Classify continuous recordings with sliding window
run('t2_testing.m')
```

> **Note:** `t1_training_2.m` contains a `return` statement before `trainNetwork`. Remove it to execute training.

### Step 2: Seizure State Decoding (HMM)

```matlab
% Update file paths at the top of hmm.m to point to your data files, then:
run('hmm.m')
```

> **Note:** Update the `filename` variables in `hmm.m` to point to the local paths of your CSV data files.

---

## Key Methods Summary

| Component | Method | Tool |
|-----------|--------|------|
| Feature extraction | Fourier Synchrosqueezed Transform (FSST) | `fsst` (Signal Processing Toolbox) |
| Waveform classification | Bidirectional LSTM | Deep Learning Toolbox |
| Seizure state modeling | Hidden Markov Model | Statistics and Machine Learning Toolbox |
| HMM parameter estimation | Baum-Welch (EM) algorithm | `hmmtrain` |
| Seizure state decoding | Viterbi algorithm | `hmmviterbi` |
| Posterior decoding | Forward-Backward algorithm | `hmmdecode` |
| Validation | 5-fold cross-validation | `cvpartition` |

---

## Citation

If you use this code or data in your research, please cite:

```bibtex
@article{oloro2024epilepsy,
  title   = {Brain-implantable needle-type CMOS imaging device enables multi-layer dissection of seizure calcium dynamics in the hippocampus},
  journal = {Journal of Neural Engineering},
  year    = {2024},
  doi     = {10.1088/1741-2552/ad5c03},
  url     = {https://iopscience.iop.org/article/10.1088/1741-2552/ad5c03}
}
```

---

## License

Please refer to the journal article and institutional guidelines for data and code licensing terms.
