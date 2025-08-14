# Master Thesis 2025 - From predicting drug response in cancer cell lines to personalized oncology

Simon PENELLE

📄 [Download the full thesis (PDF)](penelle_master-thesis_2025.pdf)

<embed src="penelle_master-thesis_2025.pdf" type="application/pdf" width="100%" height="600px" />

ABSTRACT

The main objective of this thesis is to explore whether integrating structural biology
features into predictive pipelines can improve our ability to predict cancer cell
sensitivity to treatment. Using data from existing genomic datasets like GDSC, this
work aims to go one step further by adding protein-level structural information using
modern deep learning tools. In this thesis, we focus on the BRAF gene, a
well-studied oncogene with clinically relevant mutations such as V600E, to assess
whether structural differences between wild-type and mutated forms can explain
variations in drug response. To achieve this goal, I developed a computational
pipeline combining open-source tools such as Alpha Fold 3 (for structure prediction),
PyMOL (for visualization), and Boltz 2 (for estimating ligand–protein binding affinity).
A second analysis layer based on gene expression levels was added to allow the
contribution of multiple biological dimensions (genomic, transcriptomic, and
structural) to the observed variability in drug response. The results of this analysis
show a modest correlation of 0.255 between our predictor and experimental data,
mainly due to limitations of the deep learning models used and the complex,
multidimensional and evolving nature of cancer. Ultimately, this thesis aims to
evaluate the performance of deep learning–based structural predictors in a
multi-omics drug discovery context, and to determine to what extent these
approaches can contribute to identifying, in silico, cancer cell lines that are sensitive
or resistant to therapy.


**Keywords:** Bioinformatics, Structural biology, Personalized oncology, Multi-omics,
BRAF mutations, Protein-ligand interaction, GDSC dataset, Deep learning
