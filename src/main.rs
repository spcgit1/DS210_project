mod readfile;
mod graph;

// I will the rustworkx_core, which is a "graph algorithm crate built on top of petgraph", to calculate katz centrality
use rustworkx_core::petgraph::graph::{DiGraph, NodeIndex};
use rustworkx_core::centrality::{katz_centrality};
use std::collections::{VecDeque, HashMap};

// need a function that returns if something is a gene or not to differentiate between genes and tfs. Return true if a gene or false if its a tf.
fn is_gene(node_name: &str, edges: &[(String, String)]) -> bool {
    // create an iterator that will return true if the gene in the edges slice matches the node name.
    // basically, look through all edges and check if the gene name appears as the second element in any edge tuple. If so, it returns true.
    // closure indicates to skip over the first element, and that gene is the name of the second element of the tuple.
    // the .any method will use the closure for each edge in the edges slice(because we call .iter()) 
    edges.iter().any(|(_, gene)| gene == node_name)
}

// graph may be sparse (which affects what centrality I want to measure), so I will implement BFS
// created the compute_sparsity function, which should use BFS to traverse the graph, count the edges and nodes, and then calculate the sparsity ratio
fn compute_sparsity(start_index: usize, graph: &DiGraph<String, ()>, adjacency_matrix: &[Vec<u32>]) -> f64 {
    let mut visited = vec![false; adjacency_matrix.len()];
    // initialize queue, which I will use to keep track of nodes that need to be visited next.
    // Using a vector deque(double ended queue) because it allows adding or removing elements from both ends of the queue.
    // within the double ended queue are node indices (type from rustworkx_core)
    let mut queue: VecDeque<NodeIndex> = VecDeque::new();

    // setup for BFS(mark starting node as true)
    visited[start_index] = true;
    // add the element to the end of the queue
    queue.push_back(NodeIndex::new(start_index));

    let mut edge_count = 0;
    let mut node_count = 0;

    // while loop that continues as long as elements are in the queue
    // remove and return the first element in the queue until it returns none
    while let Some(node) = queue.pop_front() {
        // keep track of node count
        node_count += 1;
        // get the row of the adjacency matrix of the current node. Each element in the row will represent whether there is an edge from te node to another node.
        for (adj, &is_edge) in adjacency_matrix[node.index()].iter().enumerate() {
            // check if there is an edge from the current node(node) to the adjacent node (adj)
            if is_edge == 1 && !visited[adj] {
                // when the adjacent node hasn't been visited and there is an edge from the current node to adjacent node, mark the adjacent node as visited, add it to the queue, and add 1 to edge_count 
                visited[adj] = true;
                queue.push_back(NodeIndex::new(adj));
                edge_count += 1;
            }
        }
    }

    // sparsity will be calculated as the ratio of existing edges to all possible edges
    let possible_edges = node_count * (node_count - 1);
    // make sure division by zero doesn't occur
    let sparsity = if possible_edges > 0 {
        edge_count as f64 / possible_edges as f64
    } else {
        0.0
    };
    sparsity
}
// create degree centrality function that will take the adjacency matrix
fn degree_centrality(adj_matrix: Vec<Vec<u32>>) -> Vec<f64> {
    let num_nodes = adj_matrix.len();
    let mut centrality_scores = vec![0;num_nodes];

    // get centrality scores
    for i in 0..num_nodes {
        centrality_scores[i] = adj_matrix[i].iter().sum::<u32>(); //sum calculated as a u32 integer 
    }
    // normalize centrality scores
    let max = (num_nodes - 1) as f64;
    // convert centrality scores into an iterator. into_iter is used to take ownership of the vector, and then pass it ot the map functio
    centrality_scores.into_iter().map(|score| score as f64 / max).collect()
}

fn main() {
    // read the data file into an adjacency list
    let edges = readfile::read_file_edges("DS210Data.txt");

    // create the reverse mapping and adjacency matrix
    let (index_to_node, indices) = graph::create_mapping(&edges);
    let adjacency_matrix = graph::create_adjacency_matrix(&edges, &indices);

    // convert the adjacency matrix to a petgraph graph
    let mut graph = DiGraph::<String, ()>::new();
    // create a hashmap to track relationship between node index and its corresponding NodeIndex in the DiGraph
    let mut node_indices = HashMap::<usize, NodeIndex>::new();
    // add_node needs to own the string, so clone is used to add the nodes to the graph
    for i in 0..index_to_node.len() {
        let node_index = graph.add_node(index_to_node[i].clone());
        // will map original index to the nodeindex in the graph
        node_indices.insert(i, node_index);
    }
    // this code is just to add the edges to the petgraph graph structure
    for (i, row) in adjacency_matrix.iter().enumerate() {
        for (j, &value) in row.iter().enumerate() {
            if value == 1 {
                let start_index = node_indices[&i];
                let dest_index = node_indices[&j];
                // the () at the end means that the edge isnt associated with a weight
                graph.add_edge(start_index, dest_index, ());
            }
        }
    }

    // ue rustworkx's katz_centrality function call. 
    let katz_result = katz_centrality(
        &graph,
        |_| Ok::<f64, ()>(1.0), // assumes all edges have the same weight. closure must return a Result<f64, ()>, f64 is the weight of the edge
        Some(0.1), // alpha
        None, // beta_map
        Some(1.0), // beta_scalar
        Some(100), // max iterations
        Some(1e-6), // tolerance
    );

    // Initialize top centrality values
    let mut top_gene_katz = ("", 0.0);
    let mut top_tf_katz = ("", 0.0);

  

    // katz centrality result is the variable holding the result from katz centrality calculation
    // it holds a Return<Option<Vec<f64>, Error>. Will either return 'Ok' value, Error, or None
    // check if katz_results is holding 'Some' value. if it is, it will make centrality_values bind to a ref of the result values(vector of centrality scores) inside the Some
    if let Ok(Some(ref centrality_values)) = katz_result {
        // combine two iterators into singular iterator of pairs -- combine  the node indices with the corresponding centrality values.
        // node_indices returns an iterator of node indices
        for (node_index, &centrality) in graph.node_indices().zip(centrality_values.iter()) {
            let node_name = &graph[node_index];
            // is the node a gene or TF
            if is_gene(node_name, &edges) {
                // if the centrality is greater than the current centrality, replace it
                if centrality > top_gene_katz.1 {
                    top_gene_katz = (node_name, centrality);
                }
            } else {
                if centrality > top_tf_katz.1 {
                    top_tf_katz = (node_name, centrality);
                }
            }
        }
    }


    println!("Top Gene (Katz Centrality): {}, Centrality: {}", top_gene_katz.0, top_gene_katz.1);
    println!("Top TF (Katz Centrality): {}, Centrality: {}", top_tf_katz.0, top_tf_katz.1);

    // now calculate degree centrality
    // clone is used because the function will take ownership
    let centrality_scores = degree_centrality(adjacency_matrix.clone());
    // will find the 5 TFs and genes with highest centrality
    let mut tf_scores = Vec::new();
    let mut gene_scores = Vec::new();
    // iterate over each index and centrality score of each node
    for (i, &score) in centrality_scores.iter().enumerate() {
        // reverse mapping
        let node_name = &index_to_node[i];
        // separate genes and TFs
        if is_gene(node_name, &edges) {
            gene_scores.push((node_name, score));
        } else {
            tf_scores.push((node_name, score));
        }
    }
    
    // sort centrality values
    // partial_cmp used because the scores are f64
    // compare a.1 (the centrality score of a, and b.1, the centrality score of b)
    // sorts in descending order
    tf_scores.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
    gene_scores.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
    
    // iterate through and take top 10
    let top_tf_scores = tf_scores.iter().take(10);
    let top_gene_scores = gene_scores.iter().take(10);
    
    println!("Top 10 TFs (Degree Centrality):");
    for (tf, score) in top_tf_scores {
        println!("{}: {}", tf, score);
    }
    
    println!("\nTop 10 Genes (Degree Centrality):");
    for (gene, score) in top_gene_scores {
        println!("{}: {}", gene, score);
    } 

    let mut count_zero_sparsity = 0;
    let mut count_nonzero_sparsity = 0;

    // compute sparsity of graph
    // 
    for i in 0..index_to_node.len() {
        // call the compute_sparsity function
        let sparsity = compute_sparsity(i, &graph, &adjacency_matrix);
        // count how many nodes have a sparsity of 0 and how many don't
        if sparsity == 0.0 {
            count_zero_sparsity += 1;
        } else {
            count_nonzero_sparsity += 1;
        }
    }

    println!("Number of nodes with zero sparsity: {}", count_zero_sparsity);
    println!("Number of nodes with non-zero sparsity: {}", count_nonzero_sparsity);
}

// Top Gene (Katz Centrality): TMEM107, Centrality: 0.015957796537280987
// Top TF (Katz Centrality): AEBP2, Centrality: 0.00978345689245355
// suggests the gene TMEM107 and the TF AEBP2 are the nodes with the highest "relatve" influence/importance according to my chosen attenuation factor.
// katz identifies "friends of friends" and incorporates indirect influence, unlike degree centrality. 

//interestingly, the gene and tf with the highest katz centrality score isn't in the top ten of gene and tf nodes for degree centrality
//most likely, the katz centrality genes may be apart of the largest network of genes/tfs, while the genes and tfs w/ the highest degree centrality scores directly interact with more nodes.
// therefore, because the dataset is from a large variety of different tissues, the tfs with the most direct interactions may indicate a "ubiquitous status", while the ones with the highest katz centrality score may indicate ones that are part of a large, specific network (tissue/gene specific network)

// the BSF returned a significant amount of nodes with zero sparsity, suggesting a lot lof nodes have extremely low connectivity relative to the possible number of edges.
// this makes sense -> there may be more recurring TFs then genes in the data set. Since the graph is directional and bipartite - genes don't connect to genes - there are a lot less edges than vertices, meaning the graph is very sparse

// compile and run the code inside the module only when running tests
#[cfg(test)]
mod tests {
    use crate::degree_centrality;
    use crate::readfile::read_file_edges; 

    // make a test to see the readfile fn works as intended
    #[test]
    fn test_read_file_edges() {
        let test_file_content = "TF GENE TISSUE\nTF1 GENE1 TISSUE\nTF2 GENE2 TISSUE";
        let test_path = "test_data.txt";
        // write my test data, unwrap called so that if the result is 'err' the test will fail.
        std::fs::write(test_path, test_file_content).unwrap();

        // call the function on the test file
        let edges = read_file_edges(test_path);
        // assert if the two expressions are equal
        assert_eq!(edges, vec![("TF1".to_string(), "GENE1".to_string()), ("TF2".to_string(), "GENE2".to_string())]);
        // remove the file. if can't be removed, will cause the test to fail
        std::fs::remove_file(test_path).unwrap();
    }
    // test the degree centrality function
    #[test]
    fn test_degree_centrality() {
        // test adjacency matrix
        let adj_matrix = vec![
            vec![0, 1, 1],  // denoted as A's connections
            vec![0, 0, 1],  // denoted as B's connections
            vec![0, 0, 0],  // denoted as C's connections
        ];


        // expected scores: A's score: 2/2, B's score: 1/2, C's score: 0/2
        let expected_scores = vec![1.0, 0.5, 0.0];

        // calculate degree centrality
        let centrality_scores = degree_centrality(adj_matrix);

        // assert the results, because floating point I am use epsilon as tolerance level
        let epsilon = 1e-6;
        // create and combine the double iterator (did this in previous code) to make pairs of values
        for (calculated, expected) in centrality_scores.iter().zip(expected_scores.iter()) {
            assert!((calculated - expected).abs() < epsilon);
        }
    }
}

