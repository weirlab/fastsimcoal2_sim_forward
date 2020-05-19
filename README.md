# What this pipeline does

This pipeline takes output from a fastSIMCOAL2 model (run separate from this pipeline) and simulates SNPs under that model to a variety of time points both past and future. The data is analyzed using Hudsons Fst and ADMIXTURE at each time point allowing the user to see how genome-wide Fst and admixture values are expected to change through time.

---

# Limitations

The pipeline is currently coded only for a fastSIMCOAL2 model with 5 populations. Recording is necessary in order to allow other numbers of populations.


---

#Setup


![OrthoFinder workflow](assets/Workflow.png)
*Figure 1: Automatic OrthoFinder analysis*

## What does OrthoFinder do?
OrthoFinder is a fast, accurate and comprehensive platform for comparative genomics. It finds **orthogroups** and **orthologs**, infers **rooted gene trees** for all orthogroups and identifies all of the **gene duplication events** in those gene trees. It also infers a **rooted species tree** for the species being analysed and maps the gene duplication events from the gene trees to
