Collaborators: None
Libaries used: https://docs.rs/rustworkx-core/latest/rustworkx_core/index.html

Name: Sean Cottrell

This project was designed to examine transcription factor and gene interaction data in order to determine what factors were ubiquitous, and which are gene/pathway specific. 

Data was sourced from http://bioinfo.life.hust.edu.cn/hTFtarget#!/download

Initially, I wanted to use eigenvector centrality in my analysis of graph sparsity, however, I quickly found out that the sparseness of the graph was way too great to give any meaningful results using eigenvector centrality. Therefore, I used bread first search.
Additionally, because the graph was incredibly sparse, I decided on using Katz centrality, as it measures "indirect influence": this is meaningful for the dataset as it relates to overall specific pathways(I.E. a TF with a large katz centrality score would most likely belong to a large network that it may not directly connect to all of the components in)
The rustworkx_core petgraph extension was used for its katz centrality function.
I calculated degree centrality to measure direct connections for TFs. More direct connections = more interactions with different genes overall, therefore indicating a ubiquitous factor (present in many different systems throughout the body).

Unsurprisingly, factors with large katz centrality are not those with large degree centrality, as katz centrality factors are most likely more pathway specific than ubiquitous factors(those with large degree centrality).

