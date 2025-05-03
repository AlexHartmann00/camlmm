# camlmm
Empirical validation of the contrast avoidance model using social media sentiment analysis

# Data acquisition

Before running the analyses, acquire the data and data_with_text folder containing the raw tweets from alex.hartmann@student.maastrichtuniversity.nl.
The obtained data will also contain shapefiles used for mapping.

# Data pre-processing

## Creating base data

The notebook *create_base_data.ipynb* will clean and filter the raw data sources. It will also create new files for the aggregated (6-hour window and daily window) datasets.

## RoBERTa

The notebook *transformer_sentiment.ipynb* will run the RoBERTa sentiment analysis for the sub-sample of users and create all necessary data files, just like *create_base_data.ipynb*

# Analyses

## Descriptives

Descriptive statistics are contained in *descriptives.ipynb*.

The mapping graphic for time zones is created in *time_zones_plot.ipynb*.

## Simple mean and standard deviation comparisons

*simple_comparisons.ipynb* and *simple_comparisons_roberta.ipynb* run all necessary analyses for this.

- raw tweet-level data comparisons
- 6-hour windowed comparisons
- daily comparisons
- plots for all analyses
- mean vs median comparisons

All of the above analyses exist for both VADER and RoBERTa-based sentiment scores in the respective notebooks.

## Analyses with constant means

The analyses of standard deviation differences when holding the mean NA levels constant can be found in *test_sd_dependence.ipynb* and *test_sd_dependence_roberta.ipynb* respectively.

## LMM analyses

The R scripts are the ones prefixed with "LMM_*".