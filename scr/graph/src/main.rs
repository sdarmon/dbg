//This code is for loading a De Bruijn graph corresponding to a composante of repeats and genes.
//Its goal is to find the subgraphs of each genes and then to see what and where gene intersects with what repeat/genes.

fn main() {
    //First let's read the arguments
    let args: Vec<String> = std::env::args().collect();
    std::env::set_var("RUST_BACKTRACE", "1");
    if args.len() != 3 {
        eprintln!("Usage: {} graph.annotated_nodes graph.edges", args[0]);
        std::process::exit(1);
    }

    //Read the graph.annotated_nodes file; it contains 8 columns separated by a tab:
    //1. Unitig ID (integer)
    //2. sequence (string)
    //3. Distance from the center of the component (integer)
    //4. Weight of the node (integer)
    //7. List of the genes that intersect the unitigs
    //5. List of the transposable elements from Dfam database that intersect the unitigs
    //6. List of the transposable elements from Repbase database that intersect the unitigs
    //8. List of the chromosome where the unitigs is aligned (and not in any genes)
    //9. Abundance of the unitig (integer)
    //10. Alignment Score of the gene(s)
    //11. Alignment Score of the transposable element(s) from Dfam

    //Read and store the graph.annotated_nodes file in dictionary with the unitig ID as key and the 7 others columns as value
    let mut nodes = std::collections::HashMap::new();
    let mut reader = csv::ReaderBuilder::new().delimiter(b'\t').has_headers(false).from_path(&args[1]).unwrap();
    for result in reader.deserialize() {
    let (id, seq, dist, weight, genes, dfam, repbase, _intergenes, abundance, _ASgen, _ASdfam): (u32, String, i32, i32, String, String, String,String, i32, i32, i32) = result.unwrap();
    nodes.insert(id, (seq, dist, weight, dfam, repbase, genes, abundance));
}

    //Read and store the graph.edges file in a dictionary with the unitig ID as key and the list of the tuples (unitigs,edge type) it is connected to as value
    //edge_type is two letters (FF, FR, RF, RR) that indicates the orientation of the edge between the two unitigs
    let mut edges = std::collections::HashMap::new();
    let mut reader = csv::ReaderBuilder::new().delimiter(b'\t').has_headers(false).from_path(&args[2]).unwrap();
    for result in reader.deserialize() {
        let (id1, id2, edge_type): (u32, u32, String) = result.unwrap();
        if !edges.contains_key(&id1) && nodes.contains_key(&id1) {
            edges.insert(id1, Vec::new());
        }
        if nodes.contains_key(&id1) && nodes.contains_key(&id2) {
            edges.get_mut(&id1).unwrap().push((id2, edge_type));
        }
    }


    //Now do the same for dfam repeats
    let mut repeats_dfam = std::collections::HashSet::new();
    for (_, (_, _, _, dfam, _, _, _)) in &nodes {
        for repeat in dfam.split("; ") {
            if repeat != "*" {
                repeats_dfam.insert(repeat.to_string());
            }
        }
    }

    for repeat in &repeats_dfam {
        get_repeats_dfam_unitigs(repeat, &nodes, &edges);
    }

    //Print a horizontal line
    println!("----------------------------------------");

    //Now do the same for repbase repeats
    let mut repeats_rb = std::collections::HashSet::new();
    for (_, (_, _, _, _, repbase, _, _)) in &nodes {
        for repeat in repbase.split("; ") {
            if repeat != "*" {
                repeats_rb.insert(repeat.to_string());
            }
        }
    }

    for repeat in &repeats_rb {
        get_repeats_repbase_unitigs(repeat, &nodes, &edges);
    }

    //Print a horizontal line
    println!("----------------------------------------");

    //Scan the nodes to collect the genes names
    let mut genes = std::collections::HashSet::new();
    for (_, (_, _, _, _, _, genes_list, _)) in &nodes {
        for gene in genes_list.split("; ") {
            let gene_split = gene.split("@").collect::<Vec<&str>>();
            let gene_name = gene_split[0];
            if gene_name != "*" {
                genes.insert(gene_split[0].to_string());
            }
        }
    }

    //For each gene, output the statistics of the subgraphs of that gene
    for gene in &genes {
        //get_gene_unitigs(gene, &nodes, &edges);
    }


}

//Define the function that given a gene name, output the statistics of the subgraphs of that gene
fn get_gene_unitigs(
    gene: &str,
    nodes: &std::collections::HashMap<u32, (String, i32, i32, String, String, String, i32)>,
    edges: &std::collections::HashMap<u32, Vec<(u32, String)>>
)
{
    let mut gene_unitigs = std::collections::HashSet::new();
    for (id, (_, _, _, _, _, genes_list, _)) in nodes {
        if genes_list.contains(gene) {
            gene_unitigs.insert(*id);
        }
    }

    //Now find the connected component of the gene unitigs by fixing a direction of an unitig
    //and then following the edges in that direction (and in reverse) until finding the sinks and
    //the sources of the component

    //Define a vector that will contain the unitigs of each connected component
    let mut components = Vec::new();
    let mut sinks_of_comp = Vec::new();
    let mut sources_of_comp = Vec::new();
    let mut visited = std::collections::HashSet::new();
    for gene_unitig in &gene_unitigs {
        if visited.contains(gene_unitig) {
            continue;
        }
        visited.insert(*gene_unitig);
        let mut component = std::collections::HashSet::new();
        let mut sinks = std::collections::HashSet::new();
        let mut sources = std::collections::HashSet::new();

        //First going to the sinks
        let mut stack = vec![(*gene_unitig, true)];
        //true means forward, false means reverse
        while !stack.is_empty() {
            let (id, forward) = stack.pop().unwrap();
            //Test if id is in the gene_unitigs set
            if !gene_unitigs.contains(&id) {
                continue;
            } else
            {
                component.insert(id);
            }
            let mut found = false;
            for (id2, edge_type) in edges.get(&id).unwrap() {
                if forward && edge_type == "FF" || !forward && edge_type == "RF" {
                    found = true;
                    if !visited.contains(id2) {
                        visited.insert(*id2);
                        stack.push((*id2, true));
                    }
                } else if forward && edge_type == "FR" || !forward && edge_type == "RR" {
                    found = true;
                    if !visited.contains(id2) {
                        visited.insert(*id2);
                        stack.push((*id2, false));
                    }
                }
            }
            if !found { //means we are at a sink
                sinks.insert((id, forward));
            }
        }

        //Then going to the sources
        let mut stack = vec![(*gene_unitig, false)];
        while !stack.is_empty() {
            let (id, forward) = stack.pop().unwrap();
            //Test if id is in the gene_unitigs set
            if !gene_unitigs.contains(&id) {
                continue;
            } else {
                component.insert(id);
            }
            let mut found = false;
            for (id2, edge_type) in edges.get(&id).unwrap() {
                if forward && edge_type == "FF" || !forward && edge_type == "RF" {
                    found = true;
                    if !visited.contains(id2) {
                        visited.insert(*id2);
                        stack.push((*id2, true));
                    }
                } else if forward && edge_type == "FR" || !forward && edge_type == "RR" {
                    found = true;
                    if !visited.contains(id2) {
                        visited.insert(*id2);
                        stack.push((*id2, false));
                    }
                }
            }
            if !found { //means we are at a source
                sources.insert((id, forward));
            }
        }

        components.push(component);
        sinks_of_comp.push(sinks);
        sources_of_comp.push(sources);
    }

    //Now try to see if there exists oriented paths from sinks/sources from different components
    //to the sources/sinks of other components. If it is the case, then merge the components into one
    //and update the sinks and sources of the new component

    //First create a map that given an unitig ID returns the index of the component it belongs to
    let mut unitig_to_comp = std::collections::HashMap::new();
    for (i, comp) in components.iter().enumerate() {
        for unitig in comp {
            unitig_to_comp.insert(*unitig, i);
        }
    }

    //Define a paths of interest vector that will contain the pairs of vertices that have a path from one the other
    let mut paths_of_interest: Vec<(usize, u32)> = Vec::new();

    //Create a map that will contain the sorthest parent of the unitigs from the sinks. Its type is i32, bool

    let mut parents: std::collections::HashMap<u32, (i32, bool)> = std::collections::HashMap::new();

    //Then using a stack, do a BFS in the whole graph to find the oriented shortest paths between the components
    //start from the sinks of the components, ignore the vertices of the components and check if we
    //reach another component. It should be either a source or a sink of another component; if it's
    //a source, then we can merge the two components; if it's a sink, then we can merge the two components
    //Only if we invert the source and sink of the second component. Do not forget to update the unitig_to_comp
    //map and the components, sinks and sources vectors.
    for i in 0..components.len() {
        //Create a stack with the sinks of the component i
        let mut stack = Vec::new();
        //Create a set that will contain the visited unitigs
        let mut visited = std::collections::HashSet::new();
        for (id, forward) in &sinks_of_comp[i] {
            stack.push((*id, *forward));
            parents.insert(*id, (-1, true));
        }
        while !stack.is_empty() {
            let (id, forward) = stack.pop().unwrap();
            visited.insert(id);
            //Check if we reach another component sinks or sources using unitig_to_comp
            let j;
            if let Some(j_ref) = unitig_to_comp.get(&id) {
                j = *j_ref;
            } else {
                //We have a gap or a cycle; add it to the paths of interest
                paths_of_interest.push((i, id));
                continue;
            }

            if j != i {
                //Add to the path of interest vector
                paths_of_interest.push((i, id));
                //Check if it's a source or a sink
                if sources_of_comp[j].contains(&(id, forward)) {
                    //Merge the two components
                    let unitigs: Vec<_> = components[j].iter().cloned().collect();
                    for unitig in unitigs {
                        unitig_to_comp.insert(unitig, i);
                        components[i].insert(unitig);
                    }
                    for (id, forward) in sinks_of_comp[j].clone() {
                        sinks_of_comp[i].insert((id, forward));
                    }
                    for (id, forward) in sources_of_comp[j].clone() {
                        sources_of_comp[i].insert((id, forward));
                    }
                    components[j].clear();
                    sinks_of_comp[j].clear();
                    sources_of_comp[j].clear();
                } else {
                    //Merge the two components
                    for unitig in components[j].clone() {
                        unitig_to_comp.insert(unitig, i);
                        components[i].insert(unitig);
                    }
                    for (id, forward) in sinks_of_comp[j].clone() {
                        sources_of_comp[i].insert((id, !(forward)));
                    }
                    for (id, forward) in sources_of_comp[j].clone() {
                        sinks_of_comp[i].insert((id, !(forward)));
                    }
                    components[j].clear();
                    sinks_of_comp[j].clear();
                    sources_of_comp[j].clear();
                }
                //Add the new sinks to the stack
                for (id, forward) in sinks_of_comp[i].clone() {
                    stack.push((id, forward));
                    parents.insert(id, (-1, true));
                }
                continue;
            }
            //Then compute the next vertices
            for (id2, edge_type) in edges.get(&id).unwrap() {
                if (forward && edge_type == "FF" || !forward && edge_type == "RF") && !visited.contains(id2) {
                    visited.insert(*id2);
                    stack.push((*id2, true));
                    parents.insert(*id2, (id.try_into().unwrap(), true));
                } else if (forward && edge_type == "FR" || !forward && edge_type == "RR") && !visited.contains(id2) {
                    visited.insert(*id2);
                    stack.push((*id2, false));
                    parents.insert(*id2, (id.try_into().unwrap(), false));
                }
            }
        }
        //Now let's do the same for sources
        let mut stack = Vec::new();
        let mut visited = std::collections::HashSet::new();
        for (id, forward) in &sources_of_comp[i] {
            stack.push((*id, *forward));
            parents.insert(*id, (-1, false));
        }
        while !stack.is_empty() {
            let (id, forward) = stack.pop().unwrap();
            visited.insert(id);
            let j ;
            if let Some(j_ref) = unitig_to_comp.get(&id) {
                j = *j_ref;
            } else {
                paths_of_interest.push((i, id));
                continue;
            }
            if j != i {
                paths_of_interest.push((i, id));
                if sinks_of_comp[j].contains(&(id, forward)) {
                    for unitig in components[j].clone() {
                        unitig_to_comp.insert(unitig, i);
                        components[i].insert(unitig);
                    }
                    for (id, forward) in sinks_of_comp[j].clone() {
                        sinks_of_comp[i].insert((id, forward));
                    }
                    for (id, forward) in sources_of_comp[j].clone() {
                        sources_of_comp[i].insert((id, forward));
                    }
                    components[j].clear();
                    sinks_of_comp[j].clear();
                    sources_of_comp[j].clear();
                } else {
                    for unitig in components[j].clone() {
                        unitig_to_comp.insert(unitig, i);
                        components[i].insert(unitig);
                    }
                    for (id, forward) in sources_of_comp[j].clone() {
                        sinks_of_comp[i].insert((id, !(forward)));
                    }
                    for (id, forward) in sinks_of_comp[j].clone() {
                        sources_of_comp[i].insert((id, !(forward)));
                    }
                    components[j].clear();
                    sinks_of_comp[j].clear();
                    sources_of_comp[j].clear();
                }
                for (id, forward) in sources_of_comp[i].clone() {
                    stack.push((id, forward));
                    parents.insert(id, (-1, false));
                }
            continue;
            }
            for (id2, edge_type) in edges.get(&id).unwrap() {
                if (forward && edge_type == "FF" || !forward && edge_type == "RF") && !visited.contains(id2) {
                    visited.insert(*id2);
                    stack.push((*id2, true));
                    parents.insert(*id2, (id.try_into().unwrap(), true));
                } else if (forward && edge_type == "FR" || !forward && edge_type == "RR") && !visited.contains(id2) {
                    visited.insert(*id2);
                    stack.push((*id2, false));
                    parents.insert(*id2, (id.try_into().unwrap(), false));
                }
            }
        }
    }
    //Now we have the components, the sinks and the sources of the components and the paths of interest
    //We can compute some statistics on each non empty component
    //Print the name of the gene and the number of non empty components in one line.
    //Then for each non empty component, print the number of different genes intersecting with the unitigs
    //and the number of different repeats intersecting with the unitigs.

    //let mut non_empty_components = 0;
    let mut comp_index = 0;

    let mut genes_intron = std::collections::HashSet::new();
    let mut genes_utr = std::collections::HashSet::new();
    let mut genes_cds = std::collections::HashSet::new();
    let mut genes_NCE = std::collections::HashSet::new();
    let mut repeats_dfam_intron = std::collections::HashSet::new();
    let mut repeats_dfam_utr = std::collections::HashSet::new();
    let mut repeats_dfam_cds = std::collections::HashSet::new();
    let mut repeats_dfam_NCE = std::collections::HashSet::new();
    let mut repeats_rb_intron = std::collections::HashSet::new();
    let mut repeats_rb_utr = std::collections::HashSet::new();
    let mut repeats_rb_cds = std::collections::HashSet::new();
    let mut repeats_rb_NCE = std::collections::HashSet::new();
    let mut genes_paths_intron = std::collections::HashSet::new();
    let mut genes_paths_utr = std::collections::HashSet::new();
    let mut genes_paths_cds = std::collections::HashSet::new();
    let mut genes_paths_NCE = std::collections::HashSet::new();
    let mut repeats_dfam_paths_intron = std::collections::HashSet::new();
    let mut repeats_dfam_paths_utr = std::collections::HashSet::new();
    let mut repeats_dfam_paths_cds = std::collections::HashSet::new();
    let mut repeats_dfam_paths_NCE = std::collections::HashSet::new();
    let mut repeats_rb_paths_intron = std::collections::HashSet::new();
    let mut repeats_rb_paths_utr = std::collections::HashSet::new();
    let mut repeats_rb_paths_cds = std::collections::HashSet::new();
    let mut repeats_rb_paths_NCE = std::collections::HashSet::new();

    for comp in &components {
        if comp.is_empty() {
            continue;
        }
        //non_empty_components += 1;
        for unitig in comp {
            if nodes.contains_key(unitig) {
                let (_, _, _, dfam, repbase, genes_list, _) = nodes.get(unitig).unwrap();
                //First find the gene to split it up in intron, UTR and CDS
                //gene should like "gene@intron", "gene@UTR" or "gene@CDS"
                let mut position : String = "".to_string();
                for genee in genes_list.split("; ") {
                    let genee_split = genee.split("@").collect::<Vec<&str>>();
                    if genee_split[0] == gene {
                        position = genee_split[1].to_string().split("%").collect::<Vec<&str>>()[0].to_string();
                        break;
                    }
                }
                //Test if position is uninitialized, there is an error print a message and continue
                if position == "" {
                    eprintln!("Error: gene {} not found in unitig {}",gene,unitig);
                    //Print the comp
                    for unitig in comp {
                        eprint!("{} ",unitig);
                    }
                    continue;
                }
                //Then add the gene to the corresponding set depending on the position
                for genee in genes_list.split("; ") {
                    let genee_split = genee.split("@").collect::<Vec<&str>>();
                    if genee_split[0] != gene && genee_split[0] != "*" {
                        if position == "intron" {
                            genes_intron.insert(genee_split[0].to_string());
                        } else if position == "UTR" {
                            genes_utr.insert(genee_split[0].to_string());
                        } else if position == "CDS" {
                            genes_cds.insert(genee_split[0].to_string());
                        } else if position == "NCE" {
                            genes_NCE.insert(genee_split[0].to_string());
                        }
                    }
                }
                //Then add the repeats to the corresponding set depending on the position
                for repeat in dfam.split("; ") {
                    if repeat == "*" {
                        continue;
                    }
                    if position == "intron" {
                        repeats_dfam_intron.insert(repeat.to_string());
                    } else if position == "UTR" {
                        repeats_dfam_utr.insert(repeat.to_string());
                    } else if position == "CDS" {
                        repeats_dfam_cds.insert(repeat.to_string());
                    } else if position == "NCE" {
                        repeats_dfam_NCE.insert(repeat.to_string());
                    }
                }
                for repeat in repbase.split("; ") {
                    if repeat == "*" {
                        continue;
                    }
                    if position == "intron" {
                        repeats_rb_intron.insert(repeat.to_string());
                    } else if position == "UTR" {
                        repeats_rb_utr.insert(repeat.to_string());
                    } else if position == "CDS" {
                        repeats_rb_cds.insert(repeat.to_string());
                    } else if position == "NCE" {
                        repeats_rb_NCE.insert(repeat.to_string());
                    }
                }
            } else {
                eprintln!("Error: unitig {} not found in nodes",unitig);
                continue;
            }
        }

        //Check if there is a path of interest between unitigs of this component and if so
        //recover the unitigs of the path.
        let mut paths: Vec<u32> = Vec::new();
        for (i, id) in &paths_of_interest {
            if *i == comp_index {
                paths.push(*id);
                while let Some((id2, _ )) = parents.get(&paths[paths.len() - 1]) {
                    if *id2 == -1 {
                        break;
                    }
                    paths.push((*id2).try_into().unwrap());
                }
            }
            //Get the gene position from the last unitig of the path
            let (_, _, _, _, _, genes_list, _) = nodes.get(&paths[paths.len() - 1]).unwrap();
            let mut position : String = "".to_string();
            for genee in genes_list.split("; ") {
                let genee_split = genee.split("@").collect::<Vec<&str>>();
                if genee_split[0] == gene {
                    position = genee_split[1].to_string().split("%").collect::<Vec<&str>>()[0].to_string();
                    break;
                }
            }
            //If position is empty, empty the paths vector and continue
            if position == "" {
                paths.clear();
                continue;
            }

            //Loop over the unitigs of the path to add the genes and repeats to the corresponding sets
            for unitig in &paths {
                let (_, _, _, dfam, repbase, genes_list, _) = nodes.get(unitig).unwrap();
                for genee in genes_list.split("; ") {
                    let genee_split = genee.split("@").collect::<Vec<&str>>();
                    if genee_split[0] != gene && genee_split[0] != "*" {
                        if position == "intron" {
                            genes_paths_intron.insert(genee_split[0].to_string());
                        } else if position == "UTR" {
                            genes_paths_utr.insert(genee_split[0].to_string());
                        } else if position == "CDS" {
                            genes_paths_cds.insert(genee_split[0].to_string());
                        } else if position == "NCE" {
                            genes_paths_NCE.insert(genee_split[0].to_string());
                        }
                    }
                }
                for repeat in dfam.split("; ") {
                    if repeat == "*" {
                        continue;
                    }
                    if position == "intron" {
                        repeats_dfam_paths_intron.insert(repeat.to_string());
                    } else if position == "UTR" {
                        repeats_dfam_paths_utr.insert(repeat.to_string());
                    } else if position == "CDS" {
                        repeats_dfam_paths_cds.insert(repeat.to_string());
                    } else if position == "NCE" {
                        repeats_dfam_paths_NCE.insert(repeat.to_string());
                    }
                }
                for repeat in repbase.split("; ") {
                    if repeat == "*" {
                        continue;
                    }
                    if position == "intron" {
                        repeats_rb_paths_intron.insert(repeat.to_string());
                    } else if position == "UTR" {
                        repeats_rb_paths_utr.insert(repeat.to_string());
                    } else if position == "CDS" {
                        repeats_rb_paths_cds.insert(repeat.to_string());
                    } else if position == "NCE" {
                        repeats_rb_paths_NCE.insert(repeat.to_string());
                    }
                }
            }
            //Clear the paths vector
            paths.clear();
        }

        //Check if everything is equal to 0 and if it's the case continue
        if genes_intron.len() == 0 && genes_utr.len() == 0 && genes_cds.len() == 0 && repeats_dfam_intron.len() == 0 && repeats_dfam_utr.len() == 0 && repeats_dfam_cds.len() == 0 && repeats_rb_intron.len() == 0 && repeats_rb_utr.len() == 0 && repeats_rb_cds.len() == 0 && genes_paths_intron.len() == 0 && genes_paths_utr.len() == 0 && genes_paths_cds.len() == 0 && repeats_dfam_paths_intron.len() == 0 && repeats_dfam_paths_utr.len() == 0 && repeats_dfam_paths_cds.len() == 0 && repeats_rb_paths_intron.len() == 0 && repeats_rb_paths_utr.len() == 0 && repeats_rb_paths_cds.len() == 0 {
            continue;
        }



        println!(">{}\t Int \t UTR \t CDS",gene);
        println!("Genes                 \t {} \t {} \t {} \t {}",genes_intron.len(),genes_utr.len(),genes_cds.len(), genes_NCE.len());
        println!("Repeats Dfam          \t {} \t {} \t {} \t {}",repeats_dfam_intron.len(),repeats_dfam_utr.len(),repeats_dfam_cds.len(), repeats_dfam_NCE.len());
        println!("Repeats Repbase       \t {} \t {} \t {} \t {}",repeats_rb_intron.len(),repeats_rb_utr.len(),repeats_rb_cds.len(), repeats_rb_NCE.len());
        println!("Genes in paths        \t {} \t {} \t {} \t {}",genes_paths_intron.len(),genes_paths_utr.len(),genes_paths_cds.len(), genes_paths_NCE.len());
        println!("Repeats Dfam in paths \t {} \t {} \t {} \t {}",repeats_dfam_paths_intron.len(),repeats_dfam_paths_utr.len(),repeats_dfam_paths_cds.len(), repeats_dfam_paths_NCE.len());
        println!("Repeats Repbase paths \t {} \t {} \t {} \t {}",repeats_rb_paths_intron.len(),repeats_rb_paths_utr.len(),repeats_rb_paths_cds.len(), repeats_rb_paths_NCE.len());


        comp_index += 1;
    }
}


//Now let define the function get_repeats_dfam_unitigs that given a repeat name from the dfam data, output the statistics of the subgraphs of that repeat
fn get_repeats_dfam_unitigs(
    repeat: &str,
    nodes: &std::collections::HashMap<u32, (String, i32, i32, String, String, String, i32)>,
    edges: &std::collections::HashMap<u32, Vec<(u32, String)>>
)
{
    let mut repeat_unitigs = std::collections::HashSet::new();
    for (id, (_, _, _, dfam, _, _, _)) in nodes {
        if dfam.contains(repeat) {
            repeat_unitigs.insert(*id);
        }
    }

    // Similar logic for finding connected components as in get_gene_unitigs
    let mut components = Vec::new();
    let mut sinks_of_comp = Vec::new();
    let mut sources_of_comp = Vec::new();
    let mut visited = std::collections::HashSet::new();
    for repeat_unitig in &repeat_unitigs {
        if visited.contains(repeat_unitig) {
            continue;
        }
        visited.insert(*repeat_unitig);
        let mut component = std::collections::HashSet::new();
        let mut sinks = std::collections::HashSet::new();
        let mut sources = std::collections::HashSet::new();

        // First go to the sinks
        let mut stack = vec![(*repeat_unitig, true)];
        while !stack.is_empty() {
            let (id, forward) = stack.pop().unwrap();
            if !repeat_unitigs.contains(&id) {
                continue;
            } else {
                component.insert(id);
            }
            let mut found = false;
            for (id2, edge_type) in edges.get(&id).unwrap() {
                if forward && edge_type == "FF" || !forward && edge_type == "RF" {
                    found = true;
                    if !visited.contains(id2) {
                        visited.insert(*id2);
                        stack.push((*id2, true));
                    }
                } else if forward && edge_type == "FR" || !forward && edge_type == "RR" {
                    found = true;
                    if !visited.contains(id2) {
                        visited.insert(*id2);
                        stack.push((*id2, false));
                    }
                }
            }
            if !found {
                sinks.insert((id, forward));
            }
        }

        // Then go to the sources
        let mut stack = vec![(*repeat_unitig, false)];
        while !stack.is_empty() {
            let (id, forward) = stack.pop().unwrap();
            if !repeat_unitigs.contains(&id) {
                continue;
            } else {
                component.insert(id);
            }
            let mut found = false;
            for (id2, edge_type) in edges.get(&id).unwrap() {
                if forward && edge_type == "FF" || !forward && edge_type == "RF" {
                    found = true;
                    if !visited.contains(id2) {
                        visited.insert(*id2);
                        stack.push((*id2, true));
                    }
                } else if forward && edge_type == "FR" || !forward && edge_type == "RR" {
                    found = true;
                    if !visited.contains(id2) {
                        visited.insert(*id2);
                        stack.push((*id2, false));
                    }
                }
            }
            if !found {
                sources.insert((id, forward));
            }
        }

        components.push(component);
        sinks_of_comp.push(sinks);
        sources_of_comp.push(sources);
    }

    // Now merge components if necessary as in the previous function
    let mut unitig_to_comp = std::collections::HashMap::new();
    for (i, comp) in components.iter().enumerate() {
        for unitig in comp {
            unitig_to_comp.insert(*unitig, i);
        }
    }

    let mut paths_of_interest: Vec<(usize, u32)> = Vec::new();
    let mut parents: std::collections::HashMap<u32, (i32, bool)> = std::collections::HashMap::new();

    // BFS to find paths between components
    for i in 0..components.len() {
        let mut stack = Vec::new();
        let mut visited = std::collections::HashSet::new();
        for (id, forward) in &sinks_of_comp[i] {
            stack.push((*id, *forward));
            parents.insert(*id, (-1, true));
        }
        while !stack.is_empty() {
            let (id, forward) = stack.pop().unwrap();
            visited.insert(id);
            let j;
            if let Some(j_ref) = unitig_to_comp.get(&id) {
                j = *j_ref;
            } else {
                paths_of_interest.push((i, id));
                continue;
            }
            if j != i {
                paths_of_interest.push((i, id));
                if sources_of_comp[j].contains(&(id, forward)) {
                    let unitigs: Vec<_> = components[j].iter().cloned().collect();
                    for unitig in unitigs {
                        unitig_to_comp.insert(unitig, i);
                        components[i].insert(unitig);
                    }
                    for (id, forward) in sinks_of_comp[j].clone() {
                        sinks_of_comp[i].insert((id, forward));
                    }
                    for (id, forward) in sources_of_comp[j].clone() {
                        sources_of_comp[i].insert((id, forward));
                    }
                    components[j].clear();
                    sinks_of_comp[j].clear();
                    sources_of_comp[j].clear();
                } else {
                    for unitig in components[j].clone() {
                        unitig_to_comp.insert(unitig, i);
                        components[i].insert(unitig);
                    }
                    for (id, forward) in sinks_of_comp[j].clone() {
                        sources_of_comp[i].insert((id, !(forward)));
                    }
                    for (id, forward) in sources_of_comp[j].clone() {
                        sinks_of_comp[i].insert((id, !(forward)));
                    }
                    components[j].clear();
                    sinks_of_comp[j].clear();
                    sources_of_comp[j].clear();
                }
                for (id, forward) in sinks_of_comp[i].clone() {
                    stack.push((id, forward));
                    parents.insert(id, (-1, true));
                }
                continue;
            }
            for (id2, edge_type) in edges.get(&id).unwrap() {
                if (forward && edge_type == "FF" || !forward && edge_type == "RF") && !visited.contains(id2) {
                    visited.insert(*id2);
                    stack.push((*id2, true));
                    parents.insert(*id2, (id.try_into().unwrap(), true));
                } else if (forward && edge_type == "FR" || !forward && edge_type == "RR") && !visited.contains(id2) {
                    visited.insert(*id2);
                    stack.push((*id2, false));
                    parents.insert(*id2, (id.try_into().unwrap(), false));
                }
            }
        }

        let mut stack = Vec::new();
        let mut visited = std::collections::HashSet::new();
        for (id, forward) in &sources_of_comp[i] {
            stack.push((*id, *forward));
            parents.insert(*id, (-1, false));
        }
        while !stack.is_empty() {
            let (id, forward) = stack.pop().unwrap();
            visited.insert(id);
            let j;
            if let Some(j_ref) = unitig_to_comp.get(&id) {
                j = *j_ref;
            } else {
                paths_of_interest.push((i, id));
                continue;
            }
            if j != i {
                paths_of_interest.push((i, id));
                if sinks_of_comp[j].contains(&(id, forward)) {
                    for unitig in components[j].clone() {
                        unitig_to_comp.insert(unitig, i);
                        components[i].insert(unitig);
                    }
                    for (id, forward) in sinks_of_comp[j].clone() {
                        sinks_of_comp[i].insert((id, forward));
                    }
                    for (id, forward) in sources_of_comp[j].clone() {
                        sources_of_comp[i].insert((id, forward));
                    }
                    components[j].clear();
                    sinks_of_comp[j].clear();
                    sources_of_comp[j].clear();
                } else {
                    for unitig in components[j].clone() {
                        unitig_to_comp.insert(unitig, i);
                        components[i].insert(unitig);
                    }
                    for (id, forward) in sources_of_comp[j].clone() {
                        sinks_of_comp[i].insert((id, !(forward)));
                    }
                    for (id, forward) in sinks_of_comp[j].clone() {
                        sources_of_comp[i].insert((id, !(forward)));
                    }
                    components[j].clear();
                    sinks_of_comp[j].clear();
                    sources_of_comp[j].clear();
                }
                for (id, forward) in sources_of_comp[i].clone() {
                    stack.push((id, forward));
                    parents.insert(id, (-1, false));
                }
            }
            for (id2, edge_type) in edges.get(&id).unwrap() {
                if (forward && edge_type == "FF" || !forward && edge_type == "RF") && !visited.contains(id2) {
                    visited.insert(*id2);
                    stack.push((*id2, true));
                    parents.insert(*id2, (id.try_into().unwrap(), true));
                } else if (forward && edge_type == "FR" || !forward && edge_type == "RR") && !visited.contains(id2) {
                    visited.insert(*id2);
                    stack.push((*id2, false));
                    parents.insert(*id2, (id.try_into().unwrap(), false));
                }
            }
        }
    }
    let mut repeats_gene_intron = std::collections::HashSet::new();
    let mut repeats_gene_utr = std::collections::HashSet::new();
    let mut repeats_gene_cds = std::collections::HashSet::new();
    let mut repeats_gene_NCE = std::collections::HashSet::new();

    let mut repeats_gene_paths_intron = std::collections::HashSet::new();
    let mut repeats_gene_paths_utr = std::collections::HashSet::new();
    let mut repeats_gene_paths_cds = std::collections::HashSet::new();
    let mut repeats_gene_paths_NCE = std::collections::HashSet::new();
    
    let mut repeats_dfam = std::collections::HashSet::new();
    let mut repeats_rb = std::collections::HashSet::new();
    
    let mut repeats_dfam_paths = std::collections::HashSet::new();
    let mut repeats_rb_paths = std::collections::HashSet::new();
    
    let mut comp_index = 0;

    for comp in &components {
        if comp.is_empty() {
            continue;
        }
        for repeat_c in comp {
            let (_, _, _, dfam, rb, genes_list, _) = nodes.get(repeat_c).unwrap();
            let mut position: String ;
            for genee in genes_list.split("; ") {
                let genee_split = genee.split("@").collect::<Vec<&str>>();
                if genee_split[0] != "*" {
                    position = genee_split[1].to_string().split("%").collect::<Vec<&str>>()[0].to_string();
                    if position == "intron" {
                        repeats_gene_intron.insert(genee_split[0].to_string());
                    } else if position == "UTR" {
                        repeats_gene_utr.insert(genee_split[0].to_string());
                    } else if position == "CDS" {
                        repeats_gene_cds.insert(genee_split[0].to_string());
                    } else if position == "NCE" {
                        repeats_gene_NCE.insert(genee_split[0].to_string());
                    }
                }
            }
            for repeat_it in dfam.split("; ") {
                if repeat_it == "*" || repeat_it == repeat {
                    continue;
                }
                repeats_dfam.insert(repeat_it.to_string());
            }

            for repeat_it in rb.split("; ") {
                if repeat_it == "*" || repeat_it == repeat {
                    continue;
                }
                repeats_rb.insert(repeat_it.to_string());
            }
        }
    }


    let mut paths: Vec<u32> = Vec::new();

    for (i, id) in &paths_of_interest {
        if *i == comp_index {
            paths.push(*id);
            while let Some((id2, _ )) = parents.get(&paths[paths.len() - 1]) {
                if *id2 == -1 {
                    break;
                }
                paths.push((*id2).try_into().unwrap());
            }
        }
        let (_, _, _, dfam, rb, genes_list, _) = nodes.get(&paths[paths.len() - 1]).unwrap();
        let mut position : String ;
        for genee in genes_list.split("; ") {
            let genee_split = genee.split("@").collect::<Vec<&str>>();
            if genee_split[0] != "*" {
                position = genee_split[1].to_string().split("%").collect::<Vec<&str>>()[0].to_string();
                if position == "intron" {
                    repeats_gene_paths_intron.insert(genee_split[0].to_string());
                } else if position == "UTR" {
                    repeats_gene_paths_utr.insert(genee_split[0].to_string());
                } else if position == "CDS" {
                    repeats_gene_paths_cds.insert(genee_split[0].to_string());
                } else if position == "NCE" {
                    repeats_gene_paths_NCE.insert(genee_split[0].to_string());
                }
            }
        }
        for repeat_it in dfam.split("; ") {
            if repeat_it == "*" || repeat_it == repeat {
                continue;
            }
            repeats_dfam_paths.insert(repeat_it.to_string());
        }

        for repeat_it in rb.split("; ") {
            if repeat_it == "*" || repeat_it == repeat {
                continue;
            }
            repeats_rb_paths.insert(repeat_it.to_string());
        }
        comp_index += 1;
    }

    //Check if everything is equal to 0 and if it's the case continue
    if repeats_gene_intron.len() == 0 && repeats_gene_utr.len() == 0 && repeats_gene_cds.len() == 0 && repeats_dfam.len() == 0 && repeats_rb.len() == 0 && repeats_gene_paths_intron.len() == 0 && repeats_gene_paths_utr.len() == 0 && repeats_gene_paths_cds.len() == 0 && repeats_dfam_paths.len() == 0 && repeats_rb_paths.len() == 0 {
        return;
    }


    println!(">{}\n                      \t Int \t UTR \t CDS \t NCE \t TE",repeat);
    println!("Genes                 \t {} \t {} \t {} \t {}",repeats_gene_intron.len(),repeats_gene_utr.len(),repeats_gene_cds.len(), repeats_gene_NCE.len());
    println!("Gene paths            \t {} \t {} \t {} \t {}",repeats_gene_paths_intron.len(),repeats_gene_paths_utr.len(),repeats_gene_paths_cds.len(), repeats_gene_paths_NCE.len());
    println!("Repeats Dfam          \t   \t   \t   \t   \t {}",repeats_dfam.len());
    println!("Repeats Repbase       \t   \t   \t   \t   \t {}",repeats_rb.len());
    println!("Repeats Dfam paths    \t   \t   \t   \t   \t {}",repeats_dfam_paths.len());
    println!("Repeats Repbase paths \t   \t   \t   \t   \t {}\n",repeats_rb_paths.len());
}

//Idem for the repbase repeats
fn get_repeats_repbase_unitigs(
    repeat: &str,
    nodes: &std::collections::HashMap<u32, (String, i32, i32, String, String, String, i32)>,
    edges: &std::collections::HashMap<u32, Vec<(u32, String)>>
)
{
    let mut repeat_unitigs = std::collections::HashSet::new();
    for (id, (_, _, _, _, repbase, _, _)) in nodes {
        if repbase.contains(repeat) {
            repeat_unitigs.insert(*id);
        }
    }

    // Similar logic for finding connected components as in get_gene_unitigs
    let mut components = Vec::new();
    let mut sinks_of_comp = Vec::new();
    let mut sources_of_comp = Vec::new();
    let mut visited = std::collections::HashSet::new();
    for repeat_unitig in &repeat_unitigs {
        if visited.contains(repeat_unitig) {
            continue;
        }
        visited.insert(*repeat_unitig);
        let mut component = std::collections::HashSet::new();
        let mut sinks = std::collections::HashSet::new();
        let mut sources = std::collections::HashSet::new();

        // First go to the sinks
        let mut stack = vec![(*repeat_unitig, true)];
        while !stack.is_empty() {
            let (id, forward) = stack.pop().unwrap();
            if !repeat_unitigs.contains(&id) {
                continue;
            } else {
                component.insert(id);
            }
            let mut found = false;
            for (id2, edge_type) in edges.get(&id).unwrap() {
                if forward && edge_type == "FF" || !forward && edge_type == "RF" {
                    found = true;
                    if !visited.contains(id2) {
                        visited.insert(*id2);
                        stack.push((*id2, true));
                    }
                } else if forward && edge_type == "FR" || !forward && edge_type == "RR" {
                    found = true;
                    if !visited.contains(id2) {
                        visited.insert(*id2);
                        stack.push((*id2, false));
                    }
                }
            }
            if !found {
                sinks.insert((id, forward));
            }
        }

        // Then go to the sources
        let mut stack = vec![(*repeat_unitig, false)];
        while !stack.is_empty() {
            let (id, forward) = stack.pop().unwrap();
            if !repeat_unitigs.contains(&id) {
                continue;
            } else {
                component.insert(id);
            }
            let mut found = false;
            for (id2, edge_type) in edges.get(&id).unwrap() {
                if forward && edge_type == "FF" || !forward && edge_type == "RF" {
                    found = true;
                    if !visited.contains(id2) {
                        visited.insert(*id2);
                        stack.push((*id2, true));
                    }
                } else if forward && edge_type == "FR" || !forward && edge_type == "RR" {
                    found = true;
                    if !visited.contains(id2) {
                        visited.insert(*id2);
                        stack.push((*id2, false));
                    }
                }
            }
            if !found {
                sources.insert((id, forward));
            }
        }

        components.push(component);
        sinks_of_comp.push(sinks);
        sources_of_comp.push(sources);
    }

    // Now merge components if necessary as in the previous function
    let mut unitig_to_comp = std::collections::HashMap::new();
    for (i, comp) in components.iter().enumerate() {
        for unitig in comp {
            unitig_to_comp.insert(*unitig, i);
        }
    }

    let mut paths_of_interest: Vec<(usize, u32)> = Vec::new();
    let mut parents: std::collections::HashMap<u32, (i32, bool)> = std::collections::HashMap::new();

    // BFS to find paths between components
    for i in 0..components.len() {
        let mut stack = Vec::new();
        let mut visited = std::collections::HashSet::new();
        for (id, forward) in &sinks_of_comp[i] {
            stack.push((*id, *forward));
            parents.insert(*id, (-1, true));
        }
        while !stack.is_empty() {
            let (id, forward) = stack.pop().unwrap();
            visited.insert(id);
            let j;
            if let Some(j_ref) = unitig_to_comp.get(&id) {
                j = *j_ref;
            } else {
                paths_of_interest.push((i, id));
                continue;
            }
            if j != i {
                paths_of_interest.push((i, id));
                if sources_of_comp[j].contains(&(id, forward)) {
                    let unitigs: Vec<_> = components[j].iter().cloned().collect();
                    for unitig in unitigs {
                        unitig_to_comp.insert(unitig, i);
                        components[i].insert(unitig);
                    }
                    for (id, forward) in sinks_of_comp[j].clone() {
                        sinks_of_comp[i].insert((id, forward));
                    }
                    for (id, forward) in sources_of_comp[j].clone() {
                        sources_of_comp[i].insert((id, forward));
                    }
                    components[j].clear();
                    sinks_of_comp[j].clear();
                    sources_of_comp[j].clear();
                } else {
                    for unitig in components[j].clone() {
                        unitig_to_comp.insert(unitig, i);
                        components[i].insert(unitig);
                    }
                    for (id, forward) in sinks_of_comp[j].clone() {
                        sources_of_comp[i].insert((id, !(forward)));
                    }
                    for (id, forward) in sources_of_comp[j].clone() {
                        sinks_of_comp[i].insert((id, !(forward)));
                    }
                    components[j].clear();
                    sinks_of_comp[j].clear();
                    sources_of_comp[j].clear();
                }
                for (id, forward) in sinks_of_comp[i].clone() {
                    stack.push((id, forward));
                    parents.insert(id, (-1, true));
                }
                continue;
            }
            for (id2, edge_type) in edges.get(&id).unwrap() {
                if (forward && edge_type == "FF" || !forward && edge_type == "RF") && !visited.contains(id2) {
                    visited.insert(*id2);
                    stack.push((*id2, true));
                    parents.insert(*id2, (id.try_into().unwrap(), true));
                } else if (forward && edge_type == "FR" || !forward && edge_type == "RR") && !visited.contains(id2) {
                    visited.insert(*id2);
                    stack.push((*id2, false));
                    parents.insert(*id2, (id.try_into().unwrap(), false));
                }
            }
        }

        let mut stack = Vec::new();
        let mut visited = std::collections::HashSet::new();
        for (id, forward) in &sources_of_comp[i] {
            stack.push((*id, *forward));
            parents.insert(*id, (-1, false));
        }
        while !stack.is_empty() {
            let (id, forward) = stack.pop().unwrap();
            visited.insert(id);
            let j;
            if let Some(j_ref) = unitig_to_comp.get(&id) {
                j = *j_ref;
            } else {
                paths_of_interest.push((i, id));
                continue;
            }
            if j != i {
                paths_of_interest.push((i, id));
                if sinks_of_comp[j].contains(&(id, forward)) {
                    for unitig in components[j].clone() {
                        unitig_to_comp.insert(unitig, i);
                        components[i].insert(unitig);
                    }
                    for (id, forward) in sinks_of_comp[j].clone() {
                        sinks_of_comp[i].insert((id, forward));
                    }
                    for (id, forward) in sources_of_comp[j].clone() {
                        sources_of_comp[i].insert((id, forward));
                    }
                    components[j].clear();
                    sinks_of_comp[j].clear();
                    sources_of_comp[j].clear();
                } else {
                    for unitig in components[j].clone() {
                        unitig_to_comp.insert(unitig, i);
                        components[i].insert(unitig);
                    }
                    for (id, forward) in sources_of_comp[j].clone() {
                        sinks_of_comp[i].insert((id, !(forward)));
                    }
                    for (id, forward) in sinks_of_comp[j].clone() {
                        sources_of_comp[i].insert((id, !(forward)));
                    }
                    components[j].clear();
                    sinks_of_comp[j].clear();
                    sources_of_comp[j].clear();
                }
                for (id, forward) in sources_of_comp[i].clone() {
                    stack.push((id, forward));
                    parents.insert(id, (-1, false));
                }
            }
            for (id2, edge_type) in edges.get(&id).unwrap() {
                if (forward && edge_type == "FF" || !forward && edge_type == "RF") && !visited.contains(id2) {
                    visited.insert(*id2);
                    stack.push((*id2, true));
                    parents.insert(*id2, (id.try_into().unwrap(), true));
                } else if (forward && edge_type == "FR" || !forward && edge_type == "RR") && !visited.contains(id2) {
                    visited.insert(*id2);
                    stack.push((*id2, false));
                    parents.insert(*id2, (id.try_into().unwrap(), false));
                }
            }
        }
    }
    let mut repeats_gene_intron = std::collections::HashSet::new();
    let mut repeats_gene_utr = std::collections::HashSet::new();
    let mut repeats_gene_cds = std::collections::HashSet::new();
    let mut repeats_gene_NCE = std::collections::HashSet::new();

    let mut repeats_gene_paths_intron = std::collections::HashSet::new();
    let mut repeats_gene_paths_utr = std::collections::HashSet::new();
    let mut repeats_gene_paths_cds = std::collections::HashSet::new();
    let mut repeats_gene_paths_NCE = std::collections::HashSet::new();

    let mut repeats_dfam = std::collections::HashSet::new();
    let mut repeats_rb = std::collections::HashSet::new();

    let mut repeats_dfam_paths = std::collections::HashSet::new();
    let mut repeats_rb_paths = std::collections::HashSet::new();

    let mut comp_index = 0;

    for comp in &components {
        if comp.is_empty() {
            continue;
        }
        for repeat_c in comp {
            let (_, _, _, dfam, rb, genes_list, _) = nodes.get(repeat_c).unwrap();
            let mut position: String ;
            for genee in genes_list.split("; ") {
                let genee_split = genee.split("@").collect::<Vec<&str>>();
                if genee_split[0] != "*" {
                    position = genee_split[1].to_string().split("%").collect::<Vec<&str>>()[0].to_string();
                    if position == "intron" {
                        repeats_gene_intron.insert(genee_split[0].to_string());
                    } else if position == "UTR" {
                        repeats_gene_utr.insert(genee_split[0].to_string());
                    } else if position == "CDS" {
                        repeats_gene_cds.insert(genee_split[0].to_string());
                    } else if position == "NCE" {
                        repeats_gene_NCE.insert(genee_split[0].to_string());
                    }
                }
            }
            for repeat_it in dfam.split("; ") {
                if repeat_it == "*" || repeat_it == repeat {
                    continue;
                }
                repeats_dfam.insert(repeat_it.to_string());
            }

            for repeat_it in rb.split("; ") {
                if repeat_it == "*" || repeat_it == repeat {
                    continue;
                }
                repeats_rb.insert(repeat_it.to_string());
            }
        }
    }

    let mut paths: Vec<u32> = Vec::new();

    for (i, id) in &paths_of_interest {
        if *i == comp_index {
            paths.push(*id);
            while let Some((id2, _ )) = parents.get(&paths[paths.len() - 1]) {
                if *id2 == -1 {
                    break;
                }
                paths.push((*id2).try_into().unwrap());
            }
        }
        let (_, _, _, dfam, rb, genes_list, _) = nodes.get(&paths[paths.len() - 1]).unwrap();
        let mut position : String ;
        for genee in genes_list.split("; ") {
            let genee_split = genee.split("@").collect::<Vec<&str>>();
            if genee_split[0] != "*" {
                position = genee_split[1].to_string().split("%").collect::<Vec<&str>>()[0].to_string();
                if position == "intron" {
                    repeats_gene_paths_intron.insert(genee_split[0].to_string());
                } else if position == "UTR" {
                    repeats_gene_paths_utr.insert(genee_split[0].to_string());
                } else if position == "CDS" {
                    repeats_gene_paths_cds.insert(genee_split[0].to_string());
                } else if position == "NCE" {
                    repeats_gene_paths_NCE.insert(genee_split[0].to_string());
                }
            }
        }
        for repeat_it in dfam.split("; ") {
            if repeat_it == "*" || repeat_it == repeat {
                continue;
            }
            repeats_dfam_paths.insert(repeat_it.to_string());
        }

        for repeat_it in rb.split("; ") {
            if repeat_it == "*" || repeat_it == repeat {
                continue;
            }
            repeats_rb_paths.insert(repeat_it.to_string());
        }
        comp_index += 1;
    }

    //Check if everything is equal to 0 and if it's the case continue
    if repeats_gene_intron.len() == 0 && repeats_gene_utr.len() == 0 && repeats_gene_cds.len() == 0 && repeats_dfam.len() == 0 && repeats_rb.len() == 0 && repeats_gene_paths_intron.len() == 0 && repeats_gene_paths_utr.len() == 0 && repeats_gene_paths_cds.len() == 0 && repeats_dfam_paths.len() == 0 && repeats_rb_paths.len() == 0 {
        return;
    }

    println!(">{}\n                      \t Int \t UTR \t CDS \t NCE \t TE",repeat);
    println!("Genes                 \t {} \t {} \t {} \t {}",repeats_gene_intron.len(),repeats_gene_utr.len(),repeats_gene_cds.len(), repeats_gene_NCE.len());
    println!("Gene paths            \t {} \t {} \t {} \t {}",repeats_gene_paths_intron.len(),repeats_gene_paths_utr.len(),repeats_gene_paths_cds.len(), repeats_gene_paths_NCE.len());
    println!("Repeats Dfam          \t   \t   \t   \t   \t {}",repeats_dfam.len());
    println!("Repeats Repbase       \t   \t   \t   \t   \t {}",repeats_rb.len());
    println!("Repeats Dfam paths    \t   \t   \t   \t   \t {}",repeats_dfam_paths.len());
    println!("Repeats Repbase paths \t   \t   \t   \t   \t {}\n",repeats_rb_paths.len());
}