# London-disease-Spatial-data-analysis-and-modelling-in-R

Two spatial models, Possion regression and CAR model with Poisson regression were applied to estimate the effects of air pollution and socio-economic variables on the risk of admission to hospital caused by respiratory disease. Both models were performed on R with OpenBUGS for fitting a Bayesian model using Markov chain Monte Carlo (MCMC) methods. .

1. Summary
	1.1  Spatial plot was created to explore the number of respiratory incident on each region in London.
	1.2  Poisson regression model using OpenBugs for MCMC method was built. Poisson distribution was selected with the 	prior of beta0, beta1, beta2 and beta3 were dnorm (0, 0.001). 
	1.3 Spatial random effects  with a Conditional-autoregressive prior distribution for Poisson regression was fitted 	using OpenBugs.
	1.4 Gelman-Rubin diagnostic (Rhat) and Geweke diagnostic were plotted to assess the converge for both models.
        1.5 Pearson residuals was calculated and used to check the model assumptions for both models.

2. Tool: R and OpenBugs.


3. Library: sped, sp, R2OpenBUGS, CARBayes, coda.


4. Dataset: Annual hospital admissions of 624 electoral wards of Greater London for  7year period spanning 2003 to 2009.
The dataset consists of 6 variables: Observed, Expected, JSA, Price, PM2.5 and year, with a total of 624 rows.