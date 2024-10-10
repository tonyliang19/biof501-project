# Biof 501 Special Topics in Bioinformatics Project

Benchmark and comparison of ML/DL methods that uses multiomics data for classification

Author: Tony Liang

## Project Outline


**Question**: Identify strengths and pitfalls of ML, DL methods from different languages (R, Python) that utilizes multiomics data (genomics, transcriptomics, proteomics, and others) from binary classification setting.

**Data input**: Various CSV files, where each file correspond to one bulk omics measurement, where n is number of observations, and p is number of features. 

**Output**: Summary report that contains visualization of results and patterns observed from classification problem, like metric performances, feature identified, and etc.

### Workflow stages/steps

1. Read in csvs and bundle them into MultiAssayExperiment (R) and MuData (Python) 
2. Preprocess these data input in a common way like removing NAs, dropping low variance columns
3. Take a common splitting out for each multiomics dataset using the MuData for cross validation usage
    - Done by stratification, likely using sklearn to the splitting
4. Passing in these splits to each method's preprocessor to prepare the data input required by method itself and do the splitting there
5. For each fold on each method, perform an inner loop of cross validation to tune for hyperparameters
    - I.e. Method 1, fold 1 do CV
    - I.e. Method 2, fold 1 do CV
    - I.e. Method 1, fold k do CV
    - I.e. Method 2, fold k do cv
    - These should be done all in parallel (from channels of nextflow)
6. Collect results for each method's folds into 1 result dataframe, repeat for all methods. Then this is wrapped more as 1 single data.
7. Visualization and summary report

