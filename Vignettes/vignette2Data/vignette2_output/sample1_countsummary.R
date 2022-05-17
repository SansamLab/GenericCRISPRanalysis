Sweave("sample1_countsummary.Rnw");
library(tools);

texi2dvi("sample1_countsummary.tex",pdf=TRUE);

