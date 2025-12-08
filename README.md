# KG-SensoryBiopolymers

This project aims to develop a multimodal model based on a knowledge graph of sensing biopolymers. It will be used to structure data to establish relationships between biological entities used in sensing, such as proteins and nucleic acids, and their relationships and activities. In contrast to the classical approach, where embeddings of monomodal models are applied, our solution will allow us to generate a holistic representation of systems covering different aspects of their functioning.

## Projects overview (submodules)

### [BiopolymersKG](https://github.com/GenerativeMolMachines/BiopolymersKG)
Development of a system using knowledge graph-based embeddings and machine learning methods.

### [MultimodalMols](https://github.com/GenerativeMolMachines/MultimodalMols)
An updated architectures for processing point clouds and graphs.

### [BindingVsEnergyPredictionModels](https://github.com/GenerativeMolMachines/BindingVsEnergyPredictionModels)
API applications of models for predicting binding, mfe and affinity between different structures and complexes to populate the graph.

### [siRNA_RPAII](https://github.com/GenerativeMolMachines/siRNA_RPAII)
Machine learning on autoencoder- and LLM-derived embeddings for the design of highly effective chemically modified siRNAs for gene knockdown

### [AntibodyAptamerGeneration](https://github.com/GenerativeMolMachines/AntibodyAptamerGeneration)
Development of a model to generate optimal aptamer sequences based on antibody and target data, offering a faster and more reliable alternative to antibodies for medical and research applications.

### [AptaBench](https://github.com/GenerativeMolMachines/AptaBench)
An unified benchmark suite for aptamer–small molecule interaction prediction

### [CyclicPeptideDesign](https://github.com/GenerativeMolMachines/CyclicPeptideDesign)
A multimodal deep learning approach for cyclic peptide-based antibiotic discovery.

## Latest updates
 - Update BiopolymersKG repository: 
 
 A scalable end-to-end system for link prediction in a multi-relational biopolymer knowledge graph was developed, including data extraction (Neo4j), preprocessing, training, inference, and evaluation. The framework integrates PyTorch-BigGraph (primary) and PyKEEN (experimental), supports multiple operators and negative sampling strategies, and provides reproducible configuration via DVC and Hydra. Comprehensive ranking and binary metrics (MRR, MAP, Hits@K, NDCG@K and Precision/Recall/F1) are implemented to assess model quality across test subsets.
 - Added AptaBench repository:

The project now includes AptaBench, a unified benchmark suite for aptamer–small molecule interaction prediction.  It provides a standardized dataset, a complete preprocessing pipeline, leakage-free evaluation splits, and updated benchmarking results for both classical and deep learning models. Reproducible training and evaluation scripts have been integrated, ensuring a consistent and comparable benchmarking workflow across experiments.

 - Added MultimodalMols repository:

We have trained architectures for processing molecular data in the format of three-dimensional point clouds and in the format of molecular graphs, so we will be able to obtain appropriate embeddings for these modalities.
 - BindingVsEnergyPredictionModels - API applications for third-party filling of the graph with new links using ML and DL.

\+ updates on personal projects

## Acknowledgements

The research was supported by ITMO University Research Projects in AI Initiative (RPAII) (project #640100).
