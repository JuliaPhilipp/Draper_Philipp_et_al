# Isoform-specific translational control is evolutionarily conserved in primates
Jolene Draper, Julia Philipp, Zach Neeb, Richard Thomas, Solomon Katzman, Sofie Salama, David Haussler,  Jeremy R. Sanford


This is the repository containing the code for the paper mnentioned above. The paper can be found on bioarchive under https://www.biorxiv.org/content/10.1101/2023.04.21.537863v2

### Abstract
Alternative splicing (AS) alters messenger RNA (mRNA) coding capacity, localization, stability, and translation. Here we use comparative transcriptomics to identify cis-acting elements coupling AS to translational control (AS-TC). We sequenced total cytosolic and polyribosome-associated mRNA from human, chimpanzee, and orangutan induced pluripotent stem cells (iPSCs), revealing thousands of transcripts with splicing differences between subcellular fractions. We found both conserved and species-specific polyribosome association patterns for orthologous splicing events. Intriguingly, alternative exons with similar polyribosome profiles between species have stronger sequence conservation than exons with lineage-specific ribosome association. These data suggest that sequence variation underlies differences in the polyribosome association. Accordingly, single nucleotide substitutions in luciferase reporters designed to model exons with divergent polyribosome profiles are sufficient to regulate translational efficiency. We used position specific weight matrixes to interpret exons with species-specific polyribosome association profiles, finding that polymorphic sites frequently alter recognition motifs for trans-acting RNA binding proteins. Together, our results show that AS can regulate translation by remodeling the cis-regulatory landscape of mRNA isoforms.

### JunctionCounts pipeline
JunctionCounts - a pipeline developed in the Sanfordlab - was run for event detection and quantification. The full code can be found here: https://github.com/ajw2329/junctionCounts. Information on how the pipeline was run can be found in `run_all.sh`

### Code for Figures & Supplementary Material

|Figure|Script|Comment|    
|---|---|---|
|Fig 1|filter_for_ASTC.Rmd | done |
|Fig 2|remap_sym_heatmaps.Rmd | check |
|Fig 3|remap_background_conservation.Rmd, remap_conservation_plots.Rmd| check |
|Fig 5|interactive_rbp_plots.Rmd| needs work |
|Supp Fig 1a|manhattan_distance_dist.Rmd| check |
|Supp Fig 1b|map_stats.R| check |
|Supp Fig 1c|genebody_coverage_plot.R | check |
|Supp Fig 1d| | missing |
|Supp Table 1| supp_tables.Rmd| check |
|Supp Table 8-11 | make_RBP_tables.Rmd | check |

### Data
Data for this publication can be found on GEO under accession number GSE230441 or [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE230441). The data is currently private and is scheduled to be released on Apr 30, 2024 or upon publication.


