# m6APrediction

Predict N6-methyladenosine (m6A) RNA modification sites from genomic and
sequence features using a fitted random forest model. The package provides two
exported functions: `prediction_multiple()` for batch prediction over a feature
table, and `prediction_single()` for predicting a single genomic position.

## Installation

```r
# install.packages("remotes")   # run once if you do not have remotes
remotes::install_github("xh22gg/m6APrediction")
