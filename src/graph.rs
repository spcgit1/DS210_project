use std::collections::HashMap;

pub type Vertex = String;
pub type OutEdges = Vec<(Vertex,Vertex)>;

// since my data is strings, and I want to work with integers within the adjacency matrix, I will create sequential integers for each string that I can look up to find the data.
pub fn create_mapping(edges: &OutEdges) -> (Vec<Vertex>, HashMap<Vertex, usize>) {
    // create hashmap (will help with both mapping and finding the number of unique vertices)
    let mut indices = HashMap::new();
    // create vec for reverse mapping lookups
    let mut index_to_node = Vec::new();
    
    // assign index values to TF and GENE strings. Keys will be strings, values will be indices (integers)
    for (tf, gene) in edges {
        // check whether TF or gene is in the indices hashmap. The "!" denotes that if the edge isn't in the hashmap, the code inside the block will run.
        if !indices.contains_key(tf) {
            // if it isn't in the hashmap, insert it (must make a clone because hashmaps will take ownership). The key's value will be the length of the HashMap.
            indices.insert(tf.clone(), indices.len());
            // then, take note of the TF's index in the vector. The index in the vector is the key's value in the hashmap 
            index_to_node.push(tf.clone());
        }

        // now do it for the genes
        if !indices.contains_key(gene) {
            indices.insert(gene.clone(), indices.len());
            index_to_node.push(gene.clone());
        }
    }
    (index_to_node, indices)
}

pub fn create_adjacency_matrix(edges: &OutEdges, indices: &HashMap<Vertex, usize>) -> Vec<Vec<u32>> {
    //the size of the adjacency matrix should be the length of the hashmap(which only has unique values) x length of hashmap
    let mut adjacency_matrix = vec![vec![0; indices.len()]; indices.len()];
    //create directed graph from TF to gene, not gene to TF)
    for (tf, gene) in edges {
        // map the connections from the same index as the key's(the TF or gene) value in the hashmap. This assures that I can pull up whatever gene or TF it is later.
        let tf_index = indices[tf];
        let gene_index = indices[gene];
        adjacency_matrix[tf_index][gene_index] = 1;
    }
    adjacency_matrix
}
