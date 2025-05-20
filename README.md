# KG-SensoryBiopolymers

This project aims to develop a multimodal model based on a knowledge graph of sensing biopolymers. It will be used to structure data to establish relationships between biological entities used in sensing, such as proteins and nucleic acids, and their relationships and activities. In contrast to the classical approach, where embeddings of monomodal models are applied, our solution will allow us to generate a holistic representation of systems covering different aspects of their functioning. In addition to knowledge graph embeddings, the following data structures and modalities are envisioned:

- molecular graphs and linear context-dependent sequences - to describe monomer linkages and structural patterns;
- diffractograms - to describe the spatial arrangement of monomers;
- tabular data - to characterize theoretically calculated and empirically measured properties and activities depending on input conditions;
- spectral data - to characterize optical properties of compounds and features of electron density distribution;
- fluorescence and electron microscopy images - to characterize more specific therapeutic and morphological properties. 

## Projects overview (submodules)

### [BiopolymersKG](https://github.com/GenerativeMolMachines/BiopolymersKG)
Development of a system using knowledge graph-based embeddings and machine learning methods for peptide analysis.

### [BindingVsEnergyPredictionModels](https://github.com/GenerativeMolMachines/BindingVsEnergyPredictionModels)
Api applications of models for predicting binding, mfe and affinity between different structures and complexes to populate the graph.

### [siRNA_RPAII](https://github.com/GenerativeMolMachines/siRNA_RPAII)
Machine learning on autoencoder- and LLM-derived embeddings for the design of highly effective chemically modified siRNAs for gene knockdown

### [AntibodyAptamerGeneration](https://github.com/GenerativeMolMachines/AntibodyAptamerGeneration)
Development of a generative model of aptamers.

### [CyclicPeptideDesign](https://github.com/GenerativeMolMachines/CyclicPeptideDesign)
A multimodal deep learning approach for cyclic peptide-based antibiotic discovery.


## Overview
[KG_embeddings.ipynb](https://github.com/GenerativeMolMachines/KG-SensoryBiopolymers/blob/main/src/learning/KG_embeddings.ipynb), executed in Google Colab using GPU Tesla T4, represents a complete pipeline for processing a Knowledge Graph (KG) based on triplet data (subject, predicate, object) loaded from a .CSV file. It includes data preparation, model training, quality evaluation, missing link prediction, embedding clustering, and graph structure analysis. 

### Environment Setup and Data Preparation 
In the first stage, the environment is set up. The dataset is loaded from Google Drive and converted into the TriplesFactory structure from pykeen, which maps unique identifiers for entities and relations. The dataset is then split into training, validation, and test sets. 

### Training the Knowledge Graph Embedding Model  
After data preparation, a knowledge graph embedding model is trained. The RotatE model is used, trained on triplets with an embedding dimension of 300, using the Adam optimizer with lr=0.001, batch_size=128, and num_epochs=100. As a result, the model learns KG embeddings, which can be used to predict missing links in the graph. 

### Model Evaluation using Rank-Based Metrics
The model’s performance is evaluated using the RankBasedEvaluator, which computes key metrics: 
- Mean Rank (MR) - the average rank of predicted triplets. 
- Mean Reciprocal Rank (MRR) - the mean reciprocal rank of predictions. 
- Hits@1, Hits@10 - the proportion of correct predictions in the top-1 and top-10 results. 

The model is tested on two example missing link prediction tasks: 
1. "Aptamers", subclass_of  ?
2. "Lysozyme", interacts_with, ?

In both cases, pykeen.predict is used to generate the top-10 most probable objects for the missing entity. 

### Clustering Entity Embeddings
After model evaluation, clustering of entity embeddings is performed. The learned embeddings are extracted and dimensionality reduction techniques are applied: 
- Truncated SVD 
- Principal Component Analysis

Next, the K-Means algorithm is used to group entities into clusters, followed by visualization of clustering results. The script also outputs lists of entities belonging to each cluster. 

### Saving the Trained Model
The trained pykeen model is saved to Google Drive, along with its metadata, for future use. 

### Knowledge Graph Structure Analysis 
In the final stage, the script performs a structural analysis of the knowledge graph, including: 
- Counting the number of entities, triplets, and unique relation types. 
- Identifying the number of classes and subclasses. 
- Calculating the degree of nodes, including the most connected entity. 
- Analyzing the distribution of connections between entities. 
- Computing the average number of connections per entity.

[create_dataset.ipynb](https://github.com/GenerativeMolMachines/KG-SensoryBiopolymers/blob/main/src/learning/create_dataset.ipynb) executed in Google Colab, where a dataset is formed based on CSV files containing data on antibodies, aptamers, molecules, and their interactions. The final dataset is structured in the form of triplets (subject, predicate, object) and saved.

Next, a list of triplets is generated:

- Relationships between antibodies and aptamers with their classes and subclasses are recorded.
- Interactions of antibodies with ligands and aptamers with target molecules are added as triplets with the predicate interacts_with.
- For molecules derived from antibodies and aptamers, relationships with their respective classes are created.
- Additionally, three fundamental relationships are manually added: Aptamers subclass_of Nucleic Acids, DNA subclass_of Nucleic Acids, RNA subclass_of Nucleic Acids.


An analysis of the resulting dataset is performed:

- General information about the dataframe is displayed, including the number of rows and data types.
- The unique number of predicates, objects, and subjects is determined.
- The frequency of each predicate in the dataset is counted.
- Unique values of objects for the class_of and subclass_of predicates are identified and displayed.

## Acknowledgements

The research was supported by ITMO University Research Projects in AI Initiative (RPAII) (project #640100).
