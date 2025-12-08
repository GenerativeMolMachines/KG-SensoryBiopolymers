# KG-SensoryBiopolymers

This project aims to develop a multimodal model based on a knowledge graph of sensing biopolymers. It will be used to structure data to establish relationships between biological entities used in sensing, such as proteins and nucleic acids, and their relationships and activities. In contrast to the classical approach, where embeddings of monomodal models are applied, our solution will allow us to generate a holistic representation of systems covering different aspects of their functioning.

## Projects overview (submodules)

### [BiopolymersKG](https://github.com/GenerativeMolMachines/BiopolymersKG)
Development of a system using knowledge graph-based embeddings and machine learning methods.

### [BindingVsEnergyPredictionModels](https://github.com/GenerativeMolMachines/BindingVsEnergyPredictionModels)
API applications of models for predicting binding, mfe and affinity between different structures and complexes to populate the graph.

### [siRNA_RPAII](https://github.com/GenerativeMolMachines/siRNA_RPAII)
Machine learning on autoencoder- and LLM-derived embeddings for the design of highly effective chemically modified siRNAs for gene knockdown

### [AntibodyAptamerGeneration](https://github.com/GenerativeMolMachines/AntibodyAptamerGeneration)
Development of a model to generate optimal aptamer sequences based on antibody and target data, offering a faster and more reliable alternative to antibodies for medical and research applications.

### [CyclicPeptideDesign](https://github.com/GenerativeMolMachines/CyclicPeptideDesign)
A multimodal deep learning approach for cyclic peptide-based antibiotic discovery.

## Latest updates
 - A scalable end-to-end system for link prediction in a multi-relational biopolymer knowledge graph was developed, including data extraction (Neo4j), preprocessing, training, inference, and evaluation. The framework integrates PyTorch-BigGraph (primary) and PyKEEN (experimental), supports multiple operators and negative sampling strategies, and provides reproducible configuration via DVC and Hydra. Comprehensive ranking and binary metrics (MRR, MAP, Hits@K, NDCG@K and Precision/Recall/F1) are implemented to assess model quality across test subsets.
 - BindingVsEnergyPredictionModels - API applications for third-party filling of the graph with new links using ML and DL.
 - BiopolymersKG - codes for database operation and statistics
 - [src/learning/KGE.ipynb](https://github.com/GenerativeMolMachines/KG-SensoryBiopolymers/blob/main/src/learning/KGE.ipynb) - preview version of the model for link prediction

\+ updates on personal projects

## Acknowledgements

The research was supported by ITMO University Research Projects in AI Initiative (RPAII) (project #640100).
