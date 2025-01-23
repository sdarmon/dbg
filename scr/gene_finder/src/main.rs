//The goal of this program is to find the genes in a graph of unitigs, starting from
//a set of unitigs that are known to be in the same component.
//Output: The list of the genes associated with the component and the subgraph
//of the neighborhood of the component.


#![allow(non_snake_case)]
fn main() {
    //First let's read the arguments
    let args: Vec<String> = std::env::args().collect();
    std::env::set_var("RUST_BACKTRACE", "1");
    use std::io::BufRead;
    use std::io::Write;
    if args.len() != 6 {
        eprintln!("Usage: {} comp_id.txt graph.annotated_nodes graph.edges output_dir k", args[0]);
        std::process::exit(1);
    }

    //Read the graph.annotated_nodes file; it contains 10 columns separated by a tab:
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

    let K = args[5].parse::<usize>().unwrap();
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

    //First, by looping over the nodes, let compute the % of unitigs covering TE
    let mut total = 0;
    let mut covered = 0;
    for idc in &comp_id {
        total += 1;
        let (_, _, _, dfam_str, _, _, _, _, _, _) = nodes.get(&idc).unwrap();
        if dfam_str != "*" {
            covered += 1;
        }
    }
    let percent = (covered as f64 / total as f64) * 100.0;

    //Set of the initial genes
    let mut genes = std::collections::HashSet::new();

    //Set of the expressed genes
    let mut exp_genes = std::collections::HashSet::new();

    //Set of the TE of the comp
    let mut te = std::collections::HashSet::new();

    //Set of leaves nodes
    let mut leaves = std::collections::HashSet::new();

    //Set of loops nodes
    //let mut loops = std::collections::HashSet::new();

    //Set of the unitigs that are already visited
    let mut visited = std::collections::HashSet::new();

    //Stack of the unitigs to visit
    let mut stack = Vec::new();

    //Let's loop over comp_id and add the first element to the stack
    for idc in comp_id.clone() {
        if !visited.contains(&idc) {
            stack.push((idc, true,0));
            //While the stack is not empty we do some bfs
            while !stack.is_empty() {
                //Pop the last element of the stack and its way
                let  (id, way,mut depth) = stack.pop().unwrap();
                //If the unitig is not visited
                if !visited.contains(&id) {
                    //Add the unitig to the visited set
                    visited.insert(id);
                    //If the unitig is in the component
                    if comp_id.contains(&id) {
                        depth = 0;
                        //Add the genes of the unitig to the genes set and the dfam TE to the te set
                        let (_, _, _, dfam_str, _, genes_str, _, _, _, _) = nodes.get(&id).unwrap();
                        //test if genes_str is equal to "*"
                        if genes_str != "*" {
                            for gene in genes_str.split("; ") {
                                genes.insert(gene.to_string());
                            }
                        }
                        if dfam_str != "*" {
                            for el in dfam_str.split("; ") {
                                te.insert(el.to_string());
                            }
                    }
                }

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
                    //If the unitig is a leaf
                    if neighbors.len() == 0 {
                        leaves.insert(id);
                    }

                    //Before looking for neighbors, check if the node is an expressed gene
                    //Such a node is a node such that there is only one gene in the genes field
                    //such that its Alignment Score is maximal and strictly greater than the other genes
                    let (_, _, _, _, _, genes_str, _, _, AS, _) = nodes.get(&id).unwrap();
                    //Test if genes_str is equal to "*"
                    if genes_str != "*" {

                        let mut max_AS = 0;
                        let mut max_gene = "";
                        let mut max_gene_AS = 0;
                        for gene in genes_str.split("; ") {
                            // A gene is a string of the form "gene_name@position€AS%coverage_rate"
                            let (gene_name, _rest) = gene.split_once('@').unwrap();
                            //Uncommenting the options allowing different AS scores.
                            //let (_position, rest2) = rest.split_once('€').unwrap();
                            //let (AS, _coverage_rate) = rest2.split_once('%').unwrap();
                            //let AS = AS.parse::<i32>().unwrap();
                            // if AS > max_AS {
                            //     max_AS = AS;
                            //     max_gene = gene_name;
                            //     max_gene_AS = 1;
                            // }
                            // if AS == max_AS {
                            //     if gene_name != max_gene {
                            //         max_gene_AS+=1;
                            //     }
                            // }
                            if max_gene == "" {
                                max_gene = gene_name;
                                max_gene_AS = 1;
                            }
                            else if gene_name != max_gene {
                                max_gene_AS+=1;
                            }
                        }
                        //Case there is only one gene with the maximal AS
                        if max_gene_AS == 1 {
                            exp_genes.insert((genes_str,depth));
                            continue;
                        }
                    }
                    //Otherwise, the node is not an expressed gene and we continue the exploration

                    //Add the neighbors of the unitig to the stack
                    for (nb, new_way) in neighbors {
                        //length of nb
                        let mut d = nodes.get(&nb).unwrap().0.len() - K + 1;
                        stack.push((nb, new_way,depth+d));
                    }
                }
                else {
                    //Case of a loop / LTR
                    //loops.insert(id);
                }
            }

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


    let mut number_genes = 0;
    let mut distances_genes_total = 0;
    let mut number_intron_genes = 0;
    let mut number_exon_genes = 0;

    file = std::fs::File::create(format!("{}/expressed_genes.txt", output_dir)).unwrap();
    for gene in exp_genes {
        let (gene_str,depth) = gene;
        //If it is a list (separator "; ") of genes, we split it and add each gene to the file
        //without duplication
        if gene_str.contains("; ") {
            for gene_el in gene_str.split("; ") {
                writeln!(file, "{}", gene_el).unwrap();
                number_genes += 1;
                distances_genes_total += depth;
                if gene_el.contains("intron") {
                    number_intron_genes += 1;
                }
                else {
                    number_exon_genes += 1;
                }
            }
        }

        writeln!(file, "{}", gene_str).unwrap();
        number_genes += 1;
        distances_genes_total += depth;
        if gene_str.contains("intron") {
            number_intron_genes += 1;
        }
        else {
            number_exon_genes += 1;
        }
    }
    //close the file
    drop(file);

    let average_distance_genes = distances_genes_total as f64 / number_genes as f64;

    file = std::fs::File::create(format!("{}/TE.txt", output_dir)).unwrap();
    for te_el in &te {
        writeln!(file, "{}", te_el).unwrap();
    }
    //close the file
    drop(file);

    file = std::fs::File::create(format!("{}/leaves.txt", output_dir)).unwrap();
    for leaf in leaves {
        writeln!(file, "{}", leaf).unwrap();
    }
    //close the file
    drop(file);
    /*
    let mut file = std::fs::File::create(format!("{}/loops.txt", output_dir)).unwrap();
    for l in loops {
        writeln!(file, "{}", l).unwrap();
    }
    //close the file
    drop(file);

     */

    //Let compute the number of distinct TE as the size of te hash set
    let mut number_te = 0;
    for _ in &te {
        number_te += 1;
    }

    //Let compute the proportion of exon genes
    let proportion_exon_genes = (number_exon_genes as f64 / number_genes as f64) * 100.0;

    //Print the percentage of unitigs covering TE
    println!("Percentage of unitigs covering TE: {:.2}%", percent);
    //Print the number of distinct TE
    println!("Number of distinct TE: {}", number_te);
    //print the mean distance of the genes
    println!("Mean distance of the genes: {:.2}", average_distance_genes);
    //print the percentage of exon genes
    println!("Percentage of exon genes: {:.2}%", proportion_exon_genes);

    //Saving those values into a file "gene_summary.txt"
    let mut file = std::fs::File::create(format!("{}/gene_summary.txt", output_dir)).unwrap();
    writeln!(file, "{:.2} \t % (Percentage of unitigs covering TE)", percent).unwrap();
    writeln!(file, "{} \t (Number of distinct TE)", number_te).unwrap();
    writeln!(file, "{:.2} \t (Mean distance of the genes)", average_distance_genes).unwrap();
    writeln!(file, "{:.2} \t % (Percentage of exon genes)", proportion_exon_genes).unwrap();
    //close the file
    drop(file);

}
