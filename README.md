# Detecting directional epistasis and dominance from cross-line analyses in Alpine populations of *Arabidopsis thaliana*

by Arnaud Le Rouzic, Marie Roumet, Alex Widmer, and Josselin Clo


## Organization

This repository contains:
* **src**, the R scripts
* **results**, where figures and tables are created
* **data**, containing the raw data for analysis

## Data

* **data/data_clean.txt** is the data file from which all analyses are based. **data_clean.txt** was generated automatically by calling the script **data/import_data.R**, which merges and fixes information from the raw data files **data/Siliq_count_final.Rform.csv** and **data/Phenotype_data_final_Rform.csv**. More information about these raw data files are in **data/README_DATA.md**. 

**data/data_clean.txt** is a tab-separated data file, the first row contains the column names, the first (unnamed) column stands for (unused) row lines. 

````{verbatim}
            Mother_line   Mother_pop   Father_line   Father_pop    Gen     Weight  Fitness
G1_P14_9    SA2-21        SA2          SA11-18       SA11          F2      0.348   926
G1_P14_4    SA2-1         SA2          SA3-6         SA3           F2      0.689   2055
G1_P14_15   SA11-18       SA11         SA11-15       SA11          F2      0.417   558
G1_P2_19    SA4-11        SA4          SA16-12       SA16          F2      0.536   916
...
````

Each line corresponds to a cross (1615 crosses documented). ``Mother_line`` and ``Mother_pop`` stand for the lineage and the population of the mother (note the redundant notation for ``Mother_line``), the same for ``Father_line`` and ``Father_pop``. The column ``Gen`` can be either ``F1`` or ``F2``, and indicates whether the line corresponds to a first or second generation cross. ``Weight`` is the dry biomass measurement, and the ``Fitness`` column gives the estimated number of siliques. 


## Scripts

* **src/fig-theor.R**     -> Produces Fig1AB.pdf and Fig3.pdf
* **src/data_analysis.R** -> Produces Fig1CD.pdf, FigS2.pdf, FigS3a.pdf, FigS3b.pdf, FigS4
* **src/run-model.R**     -> Produces Fig2.pdf, Table1.txt, TableS5.txt

The analysis was run with R version 4.1.2 (2021-11-01). 

It requires the external library ``multcomp`` and its dependencies (Hothorn T, Bretz F, Westfall P (2008). “Simultaneous Inference in General Parametric Models.” _Biometrical Journal_, *50*(3), 346-363.)

