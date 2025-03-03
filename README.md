# SGR_Flows_AMP

This analysis was used to inform the San Gabriel River adaptive management plan
for flow management in the urbanized San Gabriel River, California, USA where 
modification of discharge of treated wastewater is desired to promote water reuse, 
without adversely affecting habitat and sensitive species reliant on the flows. 

The objectives of the AMP are to ensure continuation of the pre-Project conditions (overall quality and quantity) 
of the habitat influenced by treatment plant discharges.  

To accomplish this, it is important to understand:
•	Relationships between changes in WRP discharge and habitat quality
•	Sensitivity of the selected stressor s to changes in SGR river flow
•	Relative influence of changes in WRP discharge vs. other stressor variables on habitat
•	Intensity and duration of monitoring necessary to detect habitat changes associated with flow modifications

Note that there were many iterations during the model development, which is contained
in the following scripts. Final models for both stressor indicators as in 
Irving et al xxx is located in scripts 5, 5a, 6 & 6a

## Script 00_format_data.R

data needed 
input_data/RawData_SGR.csv
input_data/Daily_Data_Dec2022.csv

Formats monitoring data - position, substrate season and variable names
Assigns groups to gages and effluent outfalls
Plots data to visualise
Formats climate data from PRISM - 34.0283_-118.0331 (https://prism.oregonstate.edu/)
calculates monthly values from hydrology and climate data
Adds in distance from outfall

Output is 00_bio_Q_matched_groups_distance.RData

## Scripts 01-04

All scripts are experimental aligning with questions and comments from Habitat Management Committee
Can be ignored

## Scripts 05-06 are final models for both indicators, absolute & relative models

## Script 05_SWP_Raw_final_model.R

Data needed
00_bio_Q_matched_groups_distance.RData

Stem Water Potential on observed values - 

Random Forest (RF) - variable importance
Mixed Model (MM) - significance of variables
Predictive Model (PM) - probability of stress 
Power analysis (PA) - duration of monitoring

RF need import of external functions - Code/Functions/Functions.R
RF cross validation at end of script (takes a long time to run)
RF & MM run with and without damaged trees
PM run on dry years and overall - names delta models in code
PM Final figure is 05_probability_of_stress_overall_current_delta_range_ALL.jpg

## Script 05a_SWP_Relative_final_model.R

Data needed
00_bio_Q_matched_groups_distance.RData

Stem Water Potential on relative values - 

Random Forest (RF) - variable importance
Mixed Model (MM) - significance of variables

RF need import of external functions - Code/Functions/Functions.R
RF cross validation at end of script (takes a long time to run)
RF & MM run with and without damaged trees

## Script 06_CV_Raw_final_model.R

Data needed
00_bio_Q_matched_groups_distance.RData

Stem Water Potential on observed values - 

Random Forest (RF) - variable importance
Mixed Model (MM) - significance of variables
Predictive Model (PM) - probability of stress 
Power analysis (PA) - duration of monitoring

RF need import of external functions - Code/Functions/Functions.R
RF cross validation at end of script (takes a long time to run)
RF & MM run with and without damaged trees
PM run on dry years and overall - names delta models in code
PM Final figure is 06_CV_probability_of_stress_overall_current_delta_range.jpg

## Script 06a_CV_Relative_final_model.R

Data needed
00_bio_Q_matched_groups_distance.RData

Stem Water Potential on relative values - 

Random Forest (RF) - variable importance
Mixed Model (MM) - significance of variables

RF need import of external functions - Code/Functions/Functions.R
RF cross validation at end of script (takes a long time to run)
RF & MM run with and without damaged trees

## Script 07_SWP_proportional_flow_sensitivity.R

Sensitivity of proportional flow to group 5 for SWP
relative im portance from RF tested
main plot - 07_relative_importance_SWP_comparison_plot.jpg

## Script 08_CV_proportional_flow_sensitivity.R

Sensitivity of proportional flow to group 5 for CV
relative importance from RF tested
main plot - 08_relative_importance_CV_comparison_plot.jpg

## Script 09_figures_for_slides.R

Figures for SCCWRP symposium slide deck
For animation purposes

## Script citations.R

Citations for packages used





