# Exploring Seed Effects in Causal Inference with Machine Learning

This repository contains code, data, and other material for a project exploring the effects of setting seeds when using double-robust methods for effect estimation with machine learning algorithms. The code in this repo should generate the following 
set of results with the NuMoM2b data:

1) Distribution of risk differences under a wide range of seeds

2) The "pseudo-bias" of extreme results: take largest and smallest 5 or 3 results, and compare them to mean/median of all results

3) Impact on standard errors: how much larger do they get if we incorporate variability due to seeds?

4) Distribution of p-values for all seeds. Impact on error rates (type 1 and type 2)?