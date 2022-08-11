# Epicross

Estimating epistasis and dominance in line cross data. 

## Data

Line cross data consists in a series of phenotypic measurements performed in two (or more) parental populations (P1 and P2) and in their first (F1) and second-generation (F2) crosses.

The number or measurements for each cross is $n_{P1}$, $n_{P2}$, $n_{F1}$, and $n_{F2}$, respectively. The phenotypic trait is noted $z$, so that $z_{P1_i}$ is the measurement in the $i$ th individual from the P1 sample. Average phenotypes ( $\bar z_{P1}$,  $\bar z_{P2}$,  $\bar z_{F1}$, and  $\bar z_{F2}$) are noted $z_{P1}$,  $z_{P2}$,  $z_{F1}$, and  $z_{F2}$ for simplicity. 

## Genetic Model

### 2 populations

It will be assumed that parental populations P1 and P2 diverge by $K$ loci, at which allele A is fixed in population P1 and allele B is fixed at population P2. We will consider the $F_\infty$ model setting, in which the reference genotype (phenotype $z_R$) is by definition the mid-phenotype between P1 and P2. In this context, considering one locus, the genotypic values are $z_R - a$ for the genotype AA, $z_R + d$ for the genotype AB, and $z_R + a$ for genotype BB, where $a$ and $d$ stand for the additive and dominance effects at this locus. With 2 loci, extra epistasis terms (noted $aa$, $ad$, $da$, and $dd$ for $A\times A$, $A\times D$, $D\times A$ and $D\times D$ epistasis)  need to be considered:


| | $A_1 A_1$ | $A_1 B_1$ | $B_1 B_1$ | 
| --- | --- | --- | --- |
| $A_2 A_2$ | $z_R - a_1 - a_2$ | $z_R + d_1 - a_2 - da- aa$ | $z_R + a_1 - a_2 - 2aa$ |
| $A_2 B_2$ | $z_R - a_1 + d_2 - ad - aa$ | $z_R + d_1 + d_2 + dd - aa$ | $z_R + a_1 + d_2 + ad - aa$ |
| $B_2 B_2$ | $z_R - a_1 + a_2 - 2aa$ | $z_R + d_1 + a_2 + da - aa$ | $z_R + a_1 + a_2$ |

This can be extended for $K$ loci, accounting for the fact that there exists an epistatic coefficient for each pair $(k,l)$ of loci. The rest will neglect epistatic coefficients of order > 2 (i.e., only pairwise interactions are considered). 

| Population | Composition | Phenotype |
| ---------- | ----------- | --------- |
| P1         | 100% $A_k A_k$ | $z_R - \sum_k a_k$ |
| P2         | 100% $B_k B_k$ | $z_R + \sum_k a_k$ |
| F1         | 100% $A_k B_k$ | $z_R + \sum_k d_k + \sum_{k,l > k} dd_{kl} - \sum_{k,l > k} aa_{kl}$ |
| F2         | 1/16 $A_k A_k, A_l A_l$ etc. | $z_R + \frac{1}{2} \sum_k d_k + \frac{1}{4} \sum_{k,l > k} dd_{kl} - \sum_{k,l > k} aa_{kl}$ |

This model illustrates the identifyability problem of epistatic terms: there are 4 data points (the mean of each cross), and 5 independent genetic quantities to estimate ( $z_R$, $A = \sum a_k$, $D=\sum d_k$, $AA = \sum aa_{kl}$, and $DD = \sum dd_{kl}$). This cannot be easily fixed by measuring additional crosses (such as backcrosses), because these will introduce new additive $\times$ dominance interaction parameters to estimate). 

The proposed solution consists in estimating the directionality of epistasis $\varepsilon$ instead of individual epistatic components. For a pair of loci $(k,l)$, $aa_{kl} = \varepsilon_{kl} a_k a_l$. $\varepsilon_{kl}$ is a constant characteristic of the pair of loci. 

In a line cross analysis, there is no information about the number of loci, the distribution of additive and epistatic effects among loci, etc. For simplicity, we will assume that $\varepsilon$ is constant, and that all loci have identical effects ( $a_1 = a_2 = \dots a_K$). Then, $\sum_{k,l>k} aa_{kl} = \sum_{k,l>k} \varepsilon_{kl} a_k a_l = \frac{1}{2} \varepsilon K(K-1) a_k^2 \simeq \frac{1}{2} \varepsilon A^2$, where $A = \sum a_k$, and assuming $K$ is large. Phenotypic values thus become:

| Population | Phenotype |
| ---------- | ----------- |
| P1         | $z_R - A$ |
| P2         | $z_R + A$ |
| F1         | $z_R + D + \frac{1}{2} \varepsilon D^2 - \frac{1}{2} \varepsilon A^2$ |
| F2         | $z_R + \frac{1}{2} D  + \frac{1}{8} \varepsilon D^2 - \frac{1}{2} \varepsilon A^2$ |

There are now 4 parameters to estimate, from 4 data points. There is no easy analytical solution to this non-linear system, but numerical solving is trivial. 

