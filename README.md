# Zeroing_out_emerging_contagions
The Dataset.zip provides epidemiological information, records of public health interventions, and a subset of auxiliary data. The auxiliary data comprises control variables obtained from GEE, namely temperature and humidity, and is available in the form of zip files: Temp***.zip and Humidity***.zip.

The ajusted reported cases and estimated Rt are stored in Estimated_Rt.csv.

The code for estimating Rt is located in the files merge_data.R, 2_NPI_group, and 3_Rt_estimation. The 4 auxdata_addition, 5 NPI_smoothed and 6 Bayesian_estimated_variant were used to run the analysis based on Bayesian inference model. Additionally, please download stanmodel.zip. The Sensitivity_analysis_bayesian.R and Sensitivity_analysis_bayesian_plot.R have the code for sensitivity analysis of Bayesian inference model and its plotting. To validate the Bayesian model and reproduce relevant figures, please download Leave_one_out_validation_Bayesian.R and Leave-one-out-validation-plot.R. 

The main code to run analysis for each outbreak using the Intervention-SEIR-Vaccination model (ISEIRV) can be found in ISEIRV_main.m. ISEIRV_Simulation_timing&intensity.m contains the code for simulating the impact of NPI timing and intensity. ISEIRV_Simulation_NPIcombined.m provides code for multiple scenario simulations of optimal intervention strategies aimed at eliminating emerging infections. To validate the ISEIRV model, you can use the code available in ISEIRV_CrossValidation.m. Additionally, the code for sensitivity analysis is provided in ISEIRV_sensitivity_analysis.m. 

To reproduce all the figures, you can use Fig***.R and ***plot.R. Note that the output figures here are meta data without embellishment.
