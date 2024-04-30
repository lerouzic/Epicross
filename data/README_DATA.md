# Detecting directional epistasis and dominance from cross-line analyses in Alpine populations of *Arabidopsis thaliana*

by Arnaud Le Rouzic, Marie Roumet, Alex Widmer, and Josselin Clo

## Data collection

The timing and procedure for data collection is extensively described in the Methods section of the article, which should serve as a reference. 

## Workflow

There are 3 original ("raw") data files, **Siliq_count_final.Rform.csv**, **Phenotype_data_final_Rform.csv**, and **Croisements.csv**, that were created during data collection, and untouched since then. **Siliq_count** and **Phenotype_data** files were simply transformed from Excel format (.xls) to comma-separated (.csv) with the library ``readxl`` in R (Hadley Wickham and Jennifer Bryan (2023). readxl: Read Excel Files. R package version 1.4.3. https://CRAN.R-project.org/package=readxl). 

These 3 files are processed by a R script, **import_data.R**, which (i) imports the three files, (ii) cleans the data (removal of missing values and duplicated lines, fixing typos in the population names), and (iii) stores a summary data file which merges information from the 3 raw data files into a new tab-separated text file, **data_clean.txt**. The content of this summary file is described in the mail README of the project. 

## Original data files

### Phenotype_data_final_Rform.csv

````
Gen  Mother  Cross  Rep     Position  Flower_bud_1  Flower_1  Imature_fruit_1  Mature_fruit_1  Comment_collection
F2   SA4-1   W-Pop  rep_1   G1_P1_1   41718         41724     41729            41745           NA
F2   SA2-2   W-Alt  rep_1   G1_P1_2   41742         41745     41750            41770           NA
F2   SA2-6   W-Pop  rep_1   G1_P1_3   41729         41734     41739            41759           NA
F2   SA17-6  B-Alt  rep_1   G1_P1_4   41753         41761     41767            41797           NA
F2   SA4-2   W-Alt  rep_1   G1_P1_5   41740         41743     41750            41770           NA
...
````

* ``Gen``: first- or second-generation cross. 
* ``Mother``: Lineage of the mother
* ``Cross``: either ``SELF`` (selfing), ``W-POP`` (within population), ``W-Alt`` (between populations of the same altitude), ``B-ALT`` (between populations of different altitudes)
* ``Rep``: replicate (1 to 7)
* ``Position``: position in the greenhouse
* ``Flower_bud_1``: day of the first flower bud (in Excel format, days from Jan 1st 1900)
* ``Flower_1``: day of the first flower (in Excel format)
* ``Immature_fruit_``: day of the first immature fruit (Excel format)
* ``Mature_fruit_1``: day of the first mature fruit (Excel format)
* ``Comment_collection``: optional comment (mostly "DEAD" or "*Collection immature*")

###  Siliq_count_final.Rform.csv

````
    Position   Total_length  Weight_3branch  Weight_tot  nb_siliq_branch1  nb_siliq_branch2  nb_siliq_branch3
1   G1_P14_9   65            0.0462          0.34766     40                45                38
2   G1_P14_4   68            0.03955         0.68871     50                35                33
3   G1_P14_15  65            0.121           0.41712     43                52                67
4   G1_P2_19   57            0.10055         0.53557     60                63                49
5   G1_P14_11  48            0.03346         0.32378     33                23                43
...

Length_branch1  Length_branch2  Length_branch3  date_data_collection  Comment
22              24              25              2014-07-23            NA
26              16              19              2014-07-23            NA
29              38              44              2014-07-23            NA
28              29              22              2014-07-23            NA
18              14              22              2014-07-23            NA
...

nbsiliq_total_reel   nb_branch_total_estim  nb_siliq_per_g    nb_siliq_per_cm   nb_siliq_total_estim  nb_siliq_total_estim/reel
NA                   7.52510822510823       2662.33766233766  1.73239436619718  925.588311688312      925.588311688312
NA                   17.4136536030341       2983.56510745891  1.9344262295082   2054.81112515803      2054.81112515803
NA                   3.44727272727273       1338.84297520661  1.45945945945946  558.458181818182      558.458181818182
NA                   5.32640477374441       1710.5917454003   2.17721518987342  916.141621084038      916.141621084038
NA                   9.67662881052003       2958.7567244471   1.83333333333333  957.986252241483      957.986252241483
...
````

* ``Position``: position in the greenhouse (makes the link with file **Phenotype_data_final_Rform.csv**)
* ``Total_length``: plant height (mm) (not used in the analysis)
* ``Weight_3branch``: weight (g) of the 3 surveyed branches
* ``Weight tot``: total weight (g)
* ``nb_siliq_branch_i``: number of siliques counted on branch i (i=1, 2, or 3)
* ``Length_branch_i``: length of branch i (i=1, 2, or 3)
* ``date_data_collection``
* ``Comment``: optional comment
* ``nbsiliq_total_reel``: NA most of the time, real number of siliques when the plant was small / less than 4 branches
* ``nb_branch_total_estim``: = ``Weight tot`` / ``Weight_3branch``
* ``nb_siliq_per_g``: = Sum(``nb_siliq_branch_i``) / ``Weight_3branch``
* ``nb_siliq_per_cm ``: = Sum(``nb_siliq_branch_i``) / Sum(``Length_branch_i``)
* ``nb_siliq_total_estim``: = ``nb_siliq_per_g`` * ``Weight tot``
* ``nb_siliq_total_estim/reel``: ``nbsiliq_total_reel`` or ``nb_siliq_total_estim``, depending on the number of branches

Note: 
* Formulae were present in the original Excel file, only the end result of the calculation is provided in the .csv

### Croisements.csv

List of the 99 types of crosses in the experiment. 

````
Mere    type_cross    Pere
SA11-11 WG            SA16-8
SA11-11 W-Pop         SA11-16
SA11-11 W-Alt         SA17-4
SA11-11 B-Alt         SA4-2
...
````

* ``Mere``: Lineage of the mother
* ``type_cross`` Same as the ``Cross`` column in **Phenotype_data_final_Rform.csv**
* ``Pere``: Lineage of the father
