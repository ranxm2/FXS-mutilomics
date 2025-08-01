# Fragile X syndrome

Fragile X syndrome (FXS), also known as Martin-Bell syndrome, is a genetic disorder that manifests mainly through mental retardation and developmental problems martin-Bell syndrome fragile X syndrome (FXS), is a genetic disorder that manifests itself mainly through mental retardation and developmental problems.

Source: https://www.invitra.com/en/fragile-x-syndrome/

# What is the Fragile X Syndrome?
Fragile X syndrome is an X-linked genetic disease caused by a mutation in the FMR-1gene. It is a change in the DNA sequence that causes the absence of a protein called FMRP1 (Fragile X Mental Retardation Protein). This protein is essential for normal brain development, which explains the problems in the development of people affected by this syndrome.

The name of the syndrome is due to the fact that this mutation causes a small break in one of the arms of the X chromosome known as the fragile site. This breakage is only apparent when the chromosome is treated under specific conditions in the laboratory.

![Comparison between normal chromosome and fragile X chromosome](./images/fragile-x-chromosome.png)

Source: https://www.invitra.com/en/fragile-x-syndrome/fragile-x-chromosome/

## Causes of Fragile X Syndrome

This genetic disease is caused by a mutation in the FMR-1 gene sequence. This change causes the protein for which this gene codes(FMRP1) is not produced and does not perform its function.

Within the FMR-1 gene there is a specific segment characterized by the repetition of a specific genetic sequence: the CGG trinucleotide. This sequence is present between 5 and 40 times in people not affected by Fragile X syndrome.

![Number of CGG triplet repeats and protein synthesis in Fragile X](./images/number-repetitions-syndrome-x-fragil.png)

Source: https://www.invitra.com/en/fragile-x-syndrome/number-repetitions-syndrome-x-fragil/

# Analysis of the Fragile X Syndrome

## Naive Test

We compare the GSVA score across different disease situaion and the drug treatment. The GSVA score is a measure of the gene expression level of the one specific pathway.

![GSVA score](./images/GSVA.png)

And we use t-test to compare the GSVA score between different disease situation 




## Mixed-Effects Model

We model **GSVA scores** using a **mixed-effects model**, incorporating **fixed effects** for disease status (FXS), drug treatments, and their **interactions**, along with a **random intercept** to account for individual-level variability:

$$
Y_{ij} = \beta_0 + \beta_1 X_{FXS} + \beta_2 X_{Drug\text{-}BAY} + \beta_3 X_{Drug\text{-}BPN} + \beta_4 X_{Drug\text{-}BP} + \beta_5 (X_{FXS} \times X_{Drug\text{-}BAY}) + \beta_6 (X_{FXS} \times X_{Drug\text{-}BPN}) + \beta_7 (X_{FXS} \times X_{Drug\text{-}BP}) + b_i + \varepsilon_{ij}
$$

### Model Components

- $Y_{ij}$: GSVA score for individual $i$ under condition $j$.

- **Fixed Effects:**
  - $X_{FXS}$: Disease status indicator (1 for **FXS**, 0 for **CTRL**).
  - $X_{Drug\text{-}BAY}$: BAY drug treatment indicator.
  - $X_{Drug\text{-}BPN}$: BPN drug treatment indicator.
  - $X_{Drug\text{-}BP}$: BP drug treatment indicator.

- **Interaction Terms:**
  - $X_{FXS} \times X_{Drug\text{-}BAY}$: Interaction between FXS and BAY.
  - $X_{FXS} \times X_{Drug\text{-}BPN}$: Interaction between FXS and BPN.
  - $X_{FXS} \times X_{Drug\text{-}BP}$: Interaction between FXS and BP.

- **Random Intercept:**
  - $b_i \sim \mathcal{N}(0, \tau^2)$: Individual-specific random intercept.

- **Error Term:**
  - $\varepsilon_{ij} \sim \mathcal{N}(0, \sigma^2)$: Residual error.

This model helps estimate both the **main effects** and **interactions** between disease status and treatment, enabling evaluation of **disease-specific drug responses**.


### Hypothesis Testing: Drug Effects Compared to CTRL in the FXS Background

### **Combined Effect Hypothesis Testing**

To assess the impact of each drug specifically within the **FXS background**, we tested whether the **combined effect** of disease status and drug treatment differs significantly from the CTRL baseline. The null hypothesis is formulated as:

$$
H_0: \beta_{\text{FXS}} + \beta_{\text{Drug}} + \beta_{\text{FXS:Drug}} = 0
$$

This hypothesis evaluates whether the **overall effect of a drug in FXS individuals** is statistically indistinguishable from the CTRL group receiving vehicle treatment — effectively testing for a **lack of rescue effect** in the FXS condition.


# BioDomain results

## LLM labeling

We used a **large language model (LLM)** to label the BioDomain for each GO term. The LLM was trained on a diverse dataset of biological literature, enabling it to understand the context and relationships between different biological entities. The model was fine-tuned to recognize and categorize GO terms based on their functional annotations.

## Mixed-Effects Model result on BioDomain
We aggregated the mixed-effects model results across all BioDomains, focusing on the **combined effect** of disease status and drug treatment. We first selet the GO term with signficant on the disease effect (FXS). Then we used the **Benjamini-Hochberg procedure** to control the false discovery rate (FDR) at a threshold of 0.05. Then we consider the difference with the interaction term (FXS:Drug) to see if the drug treatment can rescue the disease effect. The results are also visualized in the following figure.

<img src="./images/Bubble_FXS-AI.png" alt="Naive Test" width="80%">

<img src="./images/Bubble_FXS_+_BP_Rescue-AI.png" alt="Naive Test" width="80%">

<img src="./images/Bubble_FXS_+_BPN_Rescue-AI.png" alt="Naive Test" width="80%">

<img src="./images/Bubble_FXS_+_BAY_Rescue-AI.png" alt="Naive Test" width="80%">
