//The goal of this program is to find the genes in a graph of unitigs, starting from

#![allow(non_snake_case)]
//a set of unitigs that are known to be in the same component.
//Output: The list of the genes associated with the component and the subgraph
//of the neighborhood of the component.


fn main() {
    //First let's read the arguments
    let args: Vec<String> = std::env::args().collect();
    std::env::set_var("RUST_BACKTRACE", "1");
    use std::io::BufRead;
    use std::io::Write;
    if args.len() != 5 {
        eprintln!("Usage: {} comp_id.txt graph.annotated_nodes graph.edges output_dir", args[0]);
        std::process::exit(1);
    }

    //Read the graph.annotated_nodes file; it contains 8 columns separated by a tab:
    //1. Unitig ID (integer)
    //2. sequence (string)
    //3. Distance from the center of the component (integer)
    //4. Weight of the node (integer)
    //5. List of the transposable elements from Dfam database that intersect the unitigs
    //6. List of the transposable elements from Repbase database that intersect the unitigs
    //7. List of the genes that intersect the unitigs
    //8. List of the chromosome where the unitigs is aligned (and not in any genes)
    //9. Abundance of the unitig (integer)
    //10. Best alignment Score of the gene(s)
    //11. Alignment Score of the transposable element(s) from Dfam


    //Read and store the comp_id file in a set. Each line of the file is a unitig ID
    let mut comp_id = std::collections::HashSet::new();
    //Read the lines of the file, split it according to /t and store the first element, i.e.
    // the unitig ID in the set.
    let file = std::fs::File::open(&args[1]).unwrap();
    let reader = std::io::BufReader::new(file);
    for line in reader.lines() {
        let id = line.unwrap().split('\t').next().unwrap().parse::<u32>().unwrap();
        comp_id.insert(id);
    }

    //Read and store the graph.annotated_nodes file in dictionary with the unitig ID as key and the 7 others columns as value
    let mut nodes = std::collections::HashMap::new();
    let mut reader = csv::ReaderBuilder::new().delimiter(b'\t').has_headers(false).from_path(&args[2]).unwrap();
    for result in reader.deserialize() {
        let (id, seq, dist, weight, genes, dfam, repbase, intergenes, abundance, ASgen, ASdfam): (u32, String, i32, i32, String, String, String, String, i32, i32, i32) = result.unwrap();
        nodes.insert(id, (seq, dist, weight, dfam, repbase, genes, intergenes, abundance, ASgen, ASdfam));
    }

    //Read and store the graph.edges file in a dictionary with the unitig ID as key and the list of the tuples (unitigs,edge type) it is connected to as value
    //edge_type is two letters (FF, FR, RF, RR) that indicates the orientation of the edge between the two unitigs
    let mut edges = std::collections::HashMap::new();
    let mut reader = csv::ReaderBuilder::new().delimiter(b'\t').has_headers(false).from_path(&args[3]).unwrap();
    for result in reader.deserialize() {
        let (id1, id2, edge_type): (u32, u32, String) = result.unwrap();
        if !edges.contains_key(&id1) {
            edges.insert(id1, Vec::new());
        }
        edges.get_mut(&id1).unwrap().push((id2, edge_type));
    }

    //Set of the initial genes
    let mut genes = std::collections::HashSet::new();

    //Set of the expressed genes
    let mut exp_genes = std::collections::HashSet::new();

    //Set of the TE of the comp
    let mut te = std::collections::HashSet::new();

    //Set of leaves nodes
    let mut leaves = std::collections::HashSet::new();

    //Set of loops nodes
    let mut loops = std::collections::HashSet::new();

    //Set of the unitigs that are already visited
    let mut visited = std::collections::HashSet::new();

    //Stack of the unitigs to visit
    let mut stack = Vec::new();

    //Add only one unitigs and a way of the component to the stack
    let first = comp_id.iter().next().unwrap();
    stack.push((*first,true));

    //While the stack is not empty
    while !stack.is_empty() {
        //Pop the last element of the stack and its way
        let (id, way) = stack.pop().unwrap();
        //If the unitig is not visited
        if !visited.contains(&id) {
            //Add the unitig to the visited set
            visited.insert(id);
            //If the unitig is in the component
            if comp_id.contains(&id) {
                //Add the genes of the unitig to the genes set and the dfam TE to the te set
                let (_, _, _, dfam_str, _, genes_str, _, _, _, _) = nodes.get(&id).unwrap();
                for gene in genes_str.split(" ;") {
                    genes.insert(gene.to_string());
                }
                for el in dfam_str.split(" ;") {
                    te.insert(el.to_string());
                }
            }

            //If the unitig is a leaf
            //A valid neighbor is an adjacent node such that its edge type is FR/FF if way == true
            //or RF/FF if way == false
            let mut neighbors = Vec::new();
            if edges.contains_key(&id) {
                for (nb, edge_type) in edges.get(&id).unwrap() {
                    if (way && (edge_type == "FR" || edge_type == "FF")) || (!way && (edge_type == "RF" || edge_type == "FF")) {
                        let new_way :bool;
                        if edge_type == "FF" || edge_type == "RF" {
                            new_way = true;
                        }
                        else {
                            new_way = false;
                        }
                        neighbors.push((*nb, new_way));
                    }
                }
            }
            if neighbors.len() == 0 {
                leaves.insert(id);
            }

            //Before looking for neighbors, check if the node is an expressed gene
            //Such a node is a node such that there is only one gene in the genes field
            //such that its Alignment Score is maximal and strictly greater than the other genes
            let (_, _, _, _, _, genes_str, _, _, _, _) = nodes.get(&id).unwrap();
            let mut max_AS = 0;
            let mut max_gene = "";
            let mut max_gene_AS = 0;
            for gene in genes_str.split(" ;") {
                // A gene is a string of the form "gene_name@position€AS%coverage_rate"
                let (gene_name, rest) = gene.split_once('@').unwrap();
                let (_position, rest2) = rest.split_once('€').unwrap();
                let (AS, _coverage_rate) = rest2.split_once('%').unwrap();
                let AS = AS.parse::<i32>().unwrap();
                if AS > max_AS {
                    max_AS = AS;
                    max_gene = gene_name;
                    max_gene_AS = 1;
                }
                if AS == max_AS {
                    if gene_name != max_gene {
                        max_gene_AS+=1;
                    }
                }
            }
            //Case there is only one gene with the maximal AS
            if max_gene_AS == 1 {
                exp_genes.insert(max_gene);
                continue;
            }
            //Otherwise, the node is not an expressed gene and we continue the exploration

            //Add the neighbors of the unitig to the stack
            for (nb, new_way) in neighbors {
                stack.push((nb, new_way));
            }
        }
        else {
            //Case of a loop / LTR
            loops.insert(id);
        }
    }


    //FINAL : Save the five sets into 5 files of the output directory
    let output_dir = &args[4];
    let mut file = std::fs::File::create(format!("{}/initial_genes.txt", output_dir)).unwrap();
    for gene in genes {
        writeln!(file, "{}", gene).unwrap();
    }
    //close the file
    drop(file);

    file = std::fs::File::create(format!("{}/expressed_genes.txt", output_dir)).unwrap();
    for gene in exp_genes {
        writeln!(file, "{}", gene).unwrap();
    }
    //close the file
    drop(file);

    file = std::fs::File::create(format!("{}/TE.txt", output_dir)).unwrap();
    for te in te {
        writeln!(file, "{}", te).unwrap();
    }
    //close the file
    drop(file);

    file = std::fs::File::create(format!("{}/leaves.txt", output_dir)).unwrap();
    for leaf in leaves {
        writeln!(file, "{}", leaf).unwrap();
    }
    //close the file
    drop(file);

    let mut file = std::fs::File::create(format!("{}/loops.txt", output_dir)).unwrap();
    for l in loops {
        writeln!(file, "{}", l).unwrap();
    }
    //close the file
    drop(file);


}
