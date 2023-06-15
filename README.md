# Zeroing_out_emerging_contagions
The files in dataset provide epidemiological information, records of public health interventions, and a subset of auxiliary data. The auxiliary data comprises control variables obtained from GEE, namely temperature and humidity, and is available in the form of zip files: Temp***.zip and Humidity***.zip.

The ajusted reported cases and estimated Rt are stored in Estimated_Rt.csv.

The code for estimating Rt is located in the files 1 merge_data.R, 2 NPI_pre-processed.R, and 4 Rt_estimation.R. The 5 auxdata_addition.R, 6 NPI_lag_effect.R, 7 NPI_check.R, 8.1 NPI_effective_estimation_withcontrol_pooled.R and 8.2 NPI_effective_estimation_withcontrol_group.R in R_code were all used to run the analysis based on Bayesian inference model. Additionally, please download stanmodel.zip. The Sensitivity_analysis_bayesian.R and Sensitivity_analysis_bayesian_plot.R have the code for sensitivity analysis of Bayesian inference model and its plotting. To validate the Bayesian model and reproduce relevant figures, please download Leave_one_out_validation_Bayesian.R and Leave-one-out-validation-plot.R. 

The main code to run analysis for each outbreak using the Intervention-SEIR-Vaccination model (ISEIRV) can be found in ISEIRV_code. ISEIRV_Simulation_timing_intensity.m contains the code for simulating the impact of NPI timing and intensity. ISEIRV_Simulation_NPIcombined.m provides code for multiple scenario simulations of optimal intervention strategies aimed at eliminating emerging infections. To validate the ISEIRV model, you can use the code available in ISEIRV_CrossValidation.m. Additionally, the code for sensitivity analysis is provided in ISEIRV_sensitivity_analysis.m. 

To reproduce all the figures, you can use Fig***.R and ***plot.R. Note that the output figures here are meta data without embellishment.
