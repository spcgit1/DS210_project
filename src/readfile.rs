use std::fs::File;
use std::io::{BufRead, BufReader};

pub type Vertex = String; // TF and GENE names are strings

// read file into an adjacency list
pub fn read_file_edges(path: &str) -> Vec<(Vertex, Vertex)> {
    let file = File::open(path).expect("Couldn't open the file");
    // initiate bufreader
    let buf_reader = BufReader::new(file);
    // create empty adjacency list
    let mut edges = Vec::new();
    // get rid of the header
    for (index, line) in buf_reader.lines().enumerate() {
        if index == 0 {
            continue; 
        }
        // check for problems while reading
        let line = line.expect("Error reading line");
        let parts: Vec<&str> = line.trim().split_whitespace().collect();

        // another check to make sure each line has at least two entries
        if parts.len() < 2 {
            continue;
        }
        // first column should be TF
        // second column should be gene
        // push the tf gene pair into edges as a tuple 
        let tf = parts[0].to_string(); 
        let gene = parts[1].to_string(); 
        edges.push((tf, gene));
    }
    edges
}
