# xenopus_tadpole_food_availability
Code and data associated with the journal article "Food availability early in life impacts among and within individual variation in behaviour".

This readme file was generated on [22-02-2023] by [x]

GENERAL INFORMATION

Title of Dataset: Xenopus_2020_FA


Information about funding sources that supported the collection of the data: 
This project was funded by a NERC doctoral training partnership grant (X) awarded to X.  The NERC Professional Internship Programme funded collaboration with XX at X University. Funding was also provided by The School of Biology, The University of X through funding received by XX.

SHARING/ACCESS INFORMATION

Licenses/restrictions placed on the data: None

Links to publications that cite or use the data: None

Links to other publicly accessible locations of the data: 


Links/relationships to ancillary data sets: None

Was data derived from another source? No
If yes, list source(s):

Recommended citation for this dataset: 

[Author names]. 2024. Data from: Food availability early in life impacts among and within individual variation in behaviour

DATA & FILE OVERVIEW

File List:

1) Xenopus_2020.R (R script for data analysis)
2) Xenopus_ACT_EXP.txt (txt file with data to be read into R)
3) tadpole_noegg_model_yrep.rds (stan model of tadpole behavior)

Relationship between files, if important: 

Additional related data collected that was not included in the current data package: 

Are there multiple versions of the dataset?
If yes, name of file(s) that was updated: No
Why was the file updated? 
When was the file updated? 

METHODOLOGICAL INFORMATION

Description of methods used for collection/generation of data: 
We reared Xenopus laevis frog tadpoles under high or low food availability. We looked at the effect of feeding treatment on behavior by measuring the distance tadpoles swam in a familiar and unfamiliar context eight times during development.  We used a multivariate double hierarchical mixed effect model to investigate the effect of treatment on population mean behaviour and predictability as well as among individual variation in personality, plasticity and predictability.  We also looked how tadpole behaviour was correlated across contexts by looking at behavioural, plasticity and predictability syndromes. 

Methods for processing the data: Data was processed using raw data in RStudio 2021.09.1+372 

Instrument- or software-specific information needed to interpret the data: RStudio 2021.09.1+372

Standards and calibration information, if appropriate: 

Environmental/experimental conditions: 

Describe any quality-assurance procedures performed on the data: 

People involved with sample collection, processing, analysis and/or submission: 
XX (data collection, processing, analysis and submission)
XX (analysis, submission)
XX (submission)
XX (submission)

DATA-SPECIFIC INFORMATION FOR: Xenopus_ACT_EXP.txt

Number of variables: 8

Number of cases/rows: 813

Variable List:

Treatment: Name of feeding treatment - High or Low
Egg_MassID: ID of egg mass
Set: Batch tadpole testing in (1 or 2)
TadpoleID: Tadpole identification number
Rep: Trial number (1-8)
SVL: Tadpole snout vent length = body size (mm)
ACT_Total_Dist: Distance travelled in pixels in the activity assay (familar context)
EXP_Total_Dist: Distance travelled in pixels in the exploration assay (unfamilar context)

Missing data codes: None

Specialized formats or other abbreviations used: 
ACT = activity assay (familar context)
EXP = exploration assay (unfamilar context)
