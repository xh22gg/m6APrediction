# m6APrediction

Predict N6-methyladenosine (m6A) RNA modification sites from genomic and
sequence features using a fitted random forest model. The package provides two
exported functions: `prediction_multiple()` for batch prediction over a feature
table, and `prediction_single()` for predicting a single genomic position.

## Installation

```r
# install.packages("remotes")   # run once if you do not have remotes
remotes::install_github("xh22gg/m6APrediction")
```

## Usage

```r
library(m6APrediction)

ml_fit <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))
example_df <- read.csv(system.file("extdata", "m6A_input_example.csv", package = "m6APrediction"))

# Batch prediction
result <- prediction_multiple(ml_fit, head(example_df, 5))
head(result[, c("predicted_m6A_prob", "predicted_m6A_status")])

# Single-site prediction
prediction_single(ml_fit, gc_content = 0.6, RNA_type = "mRNA", RNA_region = "CDS",
                  exon_length = 12, distance_to_junction = 5,
                  evolutionary_conservation = 0.8, DNA_5mer = "GGACA")
```

## Model performance

The random forest model achieves strong discrimination of m6A sites, shown by
the ROC and precision-recall curves below.

![ROC curve](roc.png)

![PRC curve](prc.png)
