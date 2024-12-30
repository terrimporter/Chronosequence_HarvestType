# README

This repository contains the supplementary digital resources, scripts, and files used in the preprint Porter et al., 2024 "Soil communities following clearcut and salvage harvest have different early successional dynamics compared with post-wildfire patterns."

## How to cite

Please acknowledge the preprint:  
Porter, T. M., Morris, D. M., Smenderovac, E., Emilson, E. J. S., & Venier, L. (2024). Soil communities following clearcut and salvage harvest have different early successional dynamics compared with post-wildfire patterns. BioRxiv, https://doi.org/10.1101/2024.11.10.622867

## Digital Resources

The code for the MetaWorks bioinformatic pipeline that we developed to process multi-marker metabarcoding data is available from https://github.com/terrimporter/MetaWorks, installation and usage instructions are available here https://terrimporter.github.io/MetaWorksSite/, and our publication is https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0274260 .  

The following curated and trained naive Bayesian classifiers are available from:  
  ITS - https://github.com/terrimporter/UNITE_ITSClassifier, the original UNITE reference set is available from https://unite.ut.ee/repository.php  
  COI - https://github.com/terrimporter/CO1Classifier (our publication is here https://www.nature.com/articles/s41598-018-22505-4)  
  16S - comes with the Ribosomal Database Project Naive Bayesian Classifier, the original code for this is available from https://sourceforge.net/projects/rdp-classifier/  

## Infiles

## R Scripts

## References

Abarenkov, Kessy; Zirk, Allan; Piirmann, Timo; Pöhönen, Raivo; Ivanov, Filipp; Nilsson, R. Henrik; Kõljalg, Urmas (2021): Full UNITE+INSD dataset for eukaryotes. Version 10.05.2021. UNITE Community. https://doi.org/10.15156/BIO/1281567

Nilsson RH, Larsson K-H, Taylor AFS, Bengtsson-Palme J, Jeppesen TS, Schigel D, Kennedy P, Picard K, Glöckner FO, Tedersoo L, Saar I, Kõljalg U, Abarenkov K. 2018. The UNITE database for molecular identification of fungi: handling dark taxa and parallel taxonomic classifications. Nucleic Acids Research, DOI: 10.1093/nar/gky1022

If you use this dataflow or any of the provided scripts, please cite the MetaWorks paper:
Porter, T. M., & Hajibabaei, M. (2022). MetaWorks: A flexible, scalable bioinformatic pipeline for high-throughput multi-marker biodiversity assessments. PLOS ONE, 17(9), e0274260. doi: 10.1371/journal.pone.0274260

You can also site this repository: Teresita M. Porter. (2020, June 25). MetaWorks: A Multi-Marker Metabarcode Pipeline (Version v1.10.0). Zenodo. http://doi.org/10.5281/zenodo.4741407

If you use this dataflow for making COI taxonomic assignments, please cite the COI classifier publication:
Porter, T. M., & Hajibabaei, M. (2018). Automated high throughput animal CO1 metabarcode classification. Scientific Reports, 8, 4226.

If you use the pseudogene filtering methods, please cite the pseudogene publication: 
Porter, T.M., & Hajibabaei, M. (2021). Profile hidden Markov model sequence analysis can help remove putative pseudogenes from DNA barcoding and metabarcoding datasets. BMC Bioinformatics, 22: 256.

If you use the RDP classifier, please cite the publication:
Wang, Q., Garrity, G. M., Tiedje, J. M., & Cole, J. R. (2007). Naive Bayesian Classifier for Rapid Assignment of rRNA Sequences into the New Bacterial Taxonomy. Applied and Environmental Microbiology, 73(16), 5261–5267. doi:10.1128/AEM.00062-07

## Acknowledgements

I would like to acknowledge funding from the Canadian government from the Genomics Research and Development Initiative (GRDI), Metagenomics-Based Ecosystem Biomonitoring (Ecobiomics) project.

Last updated: October 27, 2024
