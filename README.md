
#  Machine Learning Prediction of Reduction Potentials Using DFT-Calculated Energy Gaps

**Random Forest Regression (RFR)** to predict ground-state reduction potentials of molecular structures using **DFT-calculated HOMO-LUMO energy gaps** and **molecular fingerprints**.

---

##  Table of Contents

- [Overview](#overview)
- [Pipeline](#pipeline)
- [Requirements](#requirements)
- [Usage](#usage)
  - [Step 1: Prepare ORCA Input Files](#step-1-prepare-orca-input-files)
  - [Step 2: Run DFT Calculations](#step-2-run-dft-calculations)
  - [Step 3: Extract Energy Gaps](#step-3-extract-energy-gaps)
  - [Step 4: Calculate Reduction Potentials](#step-4-calculate-reduction-potentials)
  - [Step 5: Train and Evaluate Models](#step-5-train-and-evaluate-models)
- [Results](#results)
  - [Model Comparison](#model-comparison)
  - [Fingerprint Comparison](#fingerprint-comparison)
  - [Cross-Validation](#cross-validation)
- [Key Findings](#key-findings)

---

## Overview

This project predicts **reduction potentials** of molecular structures by:

1. Computing **HOMO-LUMO energy gaps** via Density Functional Theory (DFT) using [ORCA](https://www.orcasoftware.de/)
2. Generating **molecular fingerprints** from SMILES representations
3. Training a **Random Forest Regressor** on the combined feature set (fingerprints + energy gaps)

The dataset contains approximately **477 molecular structures**. Multiple models and fingerprint types were benchmarked to identify the optimal configuration.

---

## Pipeline

```
XYZ Files
   │
   ▼
┌──────────────────────┐
│  Prepare ORCA Inputs  │  ← Convert XYZ to ORCA input format
└──────────┬───────────┘
           ▼
┌──────────────────────┐
│  Run DFT Calculations │  ← ORCA on server
└──────────┬───────────┘
           ▼
┌──────────────────────┐
│  Extract HOMO-LUMO    │  ← Parse ORCA output files
│  Energy Gaps          │  ← Store in Excel
└──────────┬───────────┘
           ▼
┌──────────────────────┐
│  Calculate Reduction  │  ← Derive from energy data
│  Potentials           │  ← Store in Excel
└──────────┬───────────┘
           ▼
┌──────────────────────┐
│  ML Model Training    │  ← RFR + fingerprints + energy gaps
│  & Evaluation         │  ← Cross-validation & hyperparameter tuning
└───────────────────────┘
```

---

## Requirements

```
Python >= 3.8
ORCA (for DFT calculations)
```

### Python Dependencies

```bash
pip install pandas numpy scikit-learn matplotlib rdkit openpyxl
```

| Package | Purpose |
|---|---|
| `pandas` | Data manipulation and Excel I/O |
| `numpy` | Numerical operations |
| `scikit-learn` | ML models, train/test split, cross-validation |
| `matplotlib` | Data visualisation and performance plots |
| `rdkit` | Molecular fingerprint generation (Morgan & Atom Pair) |

---

## Usage

### Step 1: Prepare ORCA Input Files

Converts `.xyz` coordinate files into ORCA `.inp` input files with appropriate charge and coordinate formatting.

```bash
python generate_DFT_input.py
```

- Reads all `*_state1.xyz` files from the `data/xyz_files/` directory
- Extracts charge information from the file name
- Generates `.inp` files with DFT calculation parameters

### Step 2: Run DFT Calculations

Batch executes ORCA calculations on a server.

```bash
./run_orca.sh
```

- Runs ORCA on each `.inp` file
- Outputs `.out` result files
- 60-second timeout between jobs

### Step 3: Extract Energy Gaps

Parses ORCA output files to extract HOMO and LUMO orbital energies and calculates the energy gap.

```bash
python database-HOMO-LUMO.py
```

- Reads orbital energies from ORCA output
- Identifies **HOMO** (last occupied orbital) and **LUMO** (first virtual orbital)
- Calculates energy gap: `E_gap = E_LUMO - E_HOMO`
- Appends results to `dataset.xlsx`

### Step 4: Calculate Reduction Potentials

Derives reduction potentials from the energy data and stores them in the dataset.

```bash
python database-redox.py
```

- Reads the HOMO-LUMO Excel file
- Splits file names to identify base molecules and states (ground, oxidised, reduced)
- Calculates reduction potentials using the defined formula
- Cleans data by removing intermediate oxidised/reduced state rows
- Saves final dataset to `dataset.xlsx`

### Step 5: Train and Evaluate Models

Trains ML models, compares performance, and performs cross-validation with hyperparameter tuning.

```bash
python learning_model.py
```

- Generates data distribution plots
- Generates molecular fingerprints (Morgan & Atom Pair)
- Trains and compares multiple models (Linear, SVR, SGD, KNN, Random Forest)
- Performs 5-fold cross-validation with grid search
- Outputs performance plots and `result.txt`

---

## Results

### Model Comparison

| Model | R² Score |
|---|---|
| Linear Regression | < 0.94 |
| Linear SVR | < 0.94 |
| SGD Regressor | < 0.94 |
| KNN Regressor | < 0.94 |
| **Random Forest Regressor** | **0.987** |

### Fingerprint Comparison

| Fingerprint | nBits | Notes |
|---|---|---|
| Morgan | 2048 | Captures local atomic environment; prone to overfitting at higher nBits |
| **Atom Pair** | **2048** | Captures interatomic distances; outperforms Morgan at optimal nBits; less overfitting |

**Atom Pair fingerprints with 2048 bits** were selected as the optimal configuration.

### Cross-Validation

5-fold cross-validation with 90/10 train/test split:

**Optimal Hyperparameters:**

| Parameter | Value |
|---|---|
| `max_depth` | None |
| `max_features` | None |
| `n_estimators` | 150 |

**Results:**

| Metric | Training (avg) | Testing (avg) |
|---|---|---|
| R² Score | Stable (~0.99) | **0.942** |

> Training R² remained stable across all 5 random sets; testing R² showed more variance but averaged at 0.942.

---

## Key Findings

- **Random Forest Regressor** significantly outperforms linear models for this task due to the **non-linear** nature of the data
- **Deep learning** was not suitable due to the small dataset size (~100 samples), high risk of overfitting, and the need for heavy regularisation
- **KNN** was too sensitive to noise in small datasets
- **Linear models** (LSVR, SGD, LR) underfit as they cannot capture non-linear relationships
- **Atom Pair fingerprints** outperform Morgan fingerprints at optimal nBits and are less prone to overfitting
- **2048 bits** is sufficient — increasing beyond this causes overfitting (Morgan) or instability (Atom Pair) due to sparse feature space
- Energy gap and reduction potential data show a **non-linear relationship**, making tree-based models ideal

---
