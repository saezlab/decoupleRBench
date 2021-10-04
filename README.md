# decoupleRBench
`decoupleRBench` allows to evaluate the performance of biological activity 
inference methods using perturbation experiments. It builds on `decoupleR`, and 
more specifically the `decouple` wrapper function. As such, it requires the 
decoupleR package to be installed and it is recommended that the user is familiar 
with the [basics of decoupleR](https://saezlab.github.io/decoupleR/articles/decoupleR.html#basics-1).
The benchmark pipeline requires an input tibble with user-specified settings,
benchmark data in the form of a count table and a corresponding metadata table.

For more information, please check:

- `decoupleRBench` vignette: https://github.com/saezlab/decoupleRBench/blob/main/vignettes/how-to-bench.Rmd
- `decoupleR` repository: https://github.com/saezlab/decoupleR
- Manuscript repository: https://github.com/saezlab/decoupleR_manuscript

## Install
To install `deocupleRBench` please run:
```
devtools::install_github('saezlab/decoupleRBench')
```

## Evaluation
For a given `decoupleR` method, activities are inferred for each regulator and 
experiment. To evaluate their performance, all experiments are concatenated 
together to generate a response vector (whether a regulator is perturbed or not)
and a predictor vector (the regulator activities). Then, using different 
thresholds we can calculate AUROC and AUPRC for each method. Given that the true 
positive classes are limited by the regulators covered in the perturbation 
experiments, we use a downsampling strategy, where for each permutation an 
equal number of  negative classes are randomly sampled.


