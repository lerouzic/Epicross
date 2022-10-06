# Epicross

Estimating epistasis and dominance in line cross data. 

### Genetic model

Usual line cross models (e.g. Lynch & Walsh 1998) are derived to compare two populations. Lynch & Walsh 1998 p 209 propose to use the F2 population as a reference in a line cross between populations $P_1$ and $P_2$:

| Population | Phenotype |
| ---------- | ----------- |
| $P_1$         | $z_{F2} + A - D + AA - AD + DD $ |
| $P_2$         | $z_{F2} - A - D + AA + AD + DD $ |
| $F1$          | $z_{F2} + D + DD $               | 
| $F2$          | $z_{F2}$                         |

Note that in this setting (no backcross populations), the additive-by-dominance epistatic effect AD cannot be distinguished from the additive effect A. Both will thus be merged in the rest of the analysis. 

When dealing with more than 2 populations, the reference cannot be the F2 between two arbitrary populations, and the setting gets closer to a diallel model. The following considers $\bar z$, the grand mean of the phenotype in the dataset, as the new reference (intercept of the model). Additive effects are population-specific, there is one $A_i$ for each population $i$, and it is possible that $A_i < 0$. Dominance and epistatic coefficients become specific of a pair of populations. The model becomes:

| Population | Phenotype |
| ---------- | ----------- |
| $P_i$         | $\bar z + A_i$ |
| $P_j$         | $\bar z + A_j$ |
| $F1_{ij}$     | $\bar z + \frac{1}{2} A_i + \frac{1}{2} A_j + 2D_{ij} - AA_{ij}$          |
| $F2_{ij}$     | $\bar z + \frac{1}{2} A_i + \frac{1}{2} A_j +  D_{ij} - AA_{ij} - DD_{ij}$ |

This full model requires 5 parameters for each pair of populations (P1, P2, F1, F2), and thus cannot be fit as such. Two strategies have been used: 
* Neglect the DxD epistatic interaction and fit a linear model with only 4 parameters.
* Assume that the directionality of epistasis $\varepsilon$ is the same for AxA and DxD epistasis, and pose $AA_{ij} = \varepsilon A_i A_j$, and $DD_{ij} = \varepsilon D_{ij}^2$. This leads to a 4-parameter non-linear model. 

