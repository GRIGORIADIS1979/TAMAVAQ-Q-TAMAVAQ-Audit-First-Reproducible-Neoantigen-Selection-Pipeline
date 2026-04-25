# Citation alignment

## NeoFox

Lang et al., *NeoFox: annotating neoantigen candidates with neoantigen features*, Bioinformatics 37(22):4246–4247, 2021, DOI: 10.1093/bioinformatics/btab344.

Use in TAMAVAQ package:

- annotate candidates with presentation and recognition features;
- preserve biological feature columns as first-class evidence;
- avoid reducing all biology to one opaque rank.

## NeoAgDT

Mösch et al., *NeoAgDT: optimization of personal neoantigen vaccine composition by digital twin*, Bioinformatics 40(5):btae205, 2024, DOI: 10.1093/bioinformatics/btae205.

Use in TAMAVAQ package:

- treat panel selection as constrained optimization;
- support tumor heterogeneity / tumor-cell coverage features;
- prefer exact ILP for production when tumor-cell simulation data are available.

## neoDesign

Yu et al., *NeoDesign: a computational tool for optimal selection of polyvalent neoantigen combinations*, Bioinformatics 40(10):btae585, 2024, DOI: 10.1093/bioinformatics/btae585.

Use in TAMAVAQ package:

- after peptide selection, optimize peptide order into a polyvalent construct;
- minimize junctional neoantigens, functional-domain-like segments, structural compactness, and linker usage;
- pass final protein sequence to lambda evaluation or LinearDesign-like downstream mRNA design.
