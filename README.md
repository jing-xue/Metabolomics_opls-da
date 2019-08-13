# Metabolomics_opls-da
Dataset generated from genetically diverse collaborative cross strains

"LVDRegression_MainEff_G03D8Str_OPLSDA.R"
- Oplsda was performed on strain-adjusted metabolite dataset. Adjustment was performed using linear regression.

"LinearRegression3diets8strains_G0Metabolomics.R" validate high VIP metabolites derived from strain-adjusted oplsda analysis.
- streamlined data transformation (box-cox for dealing with non-normal residual, tukey's ladder of powers for dealing with heterscedastic residual) and 3 different linear models in reponse to different assumption complience (OLS, OLS using heteroscedasticity consistent residuals and robust linear regression).
