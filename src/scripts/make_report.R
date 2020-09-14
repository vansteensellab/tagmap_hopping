#!/usr/bin/env Rscript
library(knitr)
library(rmarkdown)
opt = commandArgs(trailingOnly = TRUE)
root_dir = getwd()

render(opt[[1]], output_format='html_document', output_file=opt[[2]],
       knit_root_dir=root_dir)
