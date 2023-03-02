# Zeroing_out_emerging_contagions
The Dataset.zip provides epidemiological information, records of public health interventions, and a subset of auxiliary data. The auxiliary data comprises control variables obtained from GEE, namely temperature and humidity, and is available in the form of zip files: Temp***.zip and Humidity***.zip.

The ajusted reported cases and estimated Rt are stored in Estimated_Rt.csv.

The code for estimating Rt is located in the files merge_data.R, 2_NPI_group, and 3_Rt_estimation. The 4 auxdata_addition, 5 NPI_smoothed and 6 Bayesian_estimated_variant were used to run the analysis based on Bayesian inference model. Additionally, please download stanmodel.zip. 
