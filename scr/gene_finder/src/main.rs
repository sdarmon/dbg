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
if args.len() != 6 && args.len() != 7 {
        eprintln!("Usage: {} comp_id.txt/comps_prefix graph.annotated_nodes graph.edges output_dir k [nb_comps]", args[0]);
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

    //If there are only 6 arguments, we assume that the component is the first one
    if args.len() == 6 {
        //Read the lines of the file, split it according to /t and store the first element, i.e.
        // the unitig ID in the set.
        let file = std::fs::File::open(&args[1]).unwrap();
        let reader = std::io::BufReader::new(file);
        let mut comp_id = std::collections::HashSet::new();
        for line in reader.lines() {
            let id = line.unwrap().split('\t').next().unwrap().parse::<u32>().unwrap();
            comp_id.insert(id);
        }
        //First, by looping over the nodes, let compute the % of unitigs covering TE
        let mut total = 0;
        let mut covered = 0;
        let mut max_weight = 0;
        let mut max_id = 0;
        for idc in &comp_id {
            total += 1;
            let (_, _, w, dfam_str, _, _, _, _, _, _) = nodes.get(&idc).unwrap();
            if dfam_str != "*" {
                covered += 1;
            }
            if *w > max_weight {
                max_weight = *w;
                max_id = *idc;
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


        #[derive(Eq, PartialEq)]
        struct StackElem(i32, u32, bool, i32, bool);

        impl Ord for StackElem {
            fn cmp(&self, other: &Self) -> std::cmp::Ordering { //True when self.0 is less than other.0
                // Decreasing order by the first element (depth)
                other.0.cmp(&self.0) //True when other.0 is greater than self.0
                    .then_with(|| self.1.cmp(&other.1)) // For tie-breaking, compare by id
            }
        }

        impl PartialOrd for StackElem {
            fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
                Some(self.cmp(other))
            }
        }

        let mut stack = std::collections::BinaryHeap::new();
        stack.push(StackElem(0, max_id, true, 0, true));
        stack.push(StackElem(0, max_id, false, 0, false));
        //While the stack is not empty we do some bfs based on the depth

        while !stack.is_empty() {
            //Pop the first element of the stack and its way
            let StackElem(d, id, way, o_a, start_way) = stack.pop().unwrap();
            let mut depth = d;
            let mut old_ab = o_a;
            //If the unitig is not visited
            if !visited.contains(&(id,way)) {
                //Add the unitig to the visited set
                visited.insert((id,way));
                //If the unitig is in the component
                if comp_id.contains(&id) {
                    depth = 0;
                    //Add the genes of the unitig to the genes set and the dfam TE to the te set
                    let (_, _, _, dfam_str, _, genes_str, _, ab, _, _) = nodes.get(&id).unwrap();
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
                    old_ab = *ab;
                }

                //A valid neighbor is an adjacent node such that its edge type is FR/FF if way == true
                //or RF/RR if way == false
                let mut neighbors = Vec::new();

                if edges.contains_key(&id) {
                    for (nb, edge_type) in edges.get(&id).unwrap() {
                        if (way && (edge_type == "FR" || edge_type == "FF")) || (!way && (edge_type == "RF" || edge_type == "RR")) {
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
                    leaves.insert((id, depth, old_ab));
                }

                //Before looking for neighbors, check if the node is an expressed gene
                //Such a node is a node such that there is only one gene in the genes field
                //such that its Alignment Score is maximal and strictly greater than the other genes
                let (seq, _, _, dfam_str, _, genes_str, inter_str, ab, _, _) = nodes.get(&id).unwrap();
                //Test if genes_str is equal to "*"
                if dfam_str != "*" {
                    for el in dfam_str.split("; ") {
                        //Add the TE of the unitig to the te set with
                        //the distance of the unitig from the center of the component
                        te.insert(format!("{}\t{}", el, depth));
                    }
                }
                if genes_str != "*" {

                    //let mut max_AS = 0;
                    let mut max_gene = "";
                    let mut max_gene_AS = 0;
                    for gene in genes_str.split("; ") {
                        // A gene is a string of the form "gene_name@position€AS%coverage_rate"
                        let (gene_name, _rest) = gene.split_once('@').unwrap();
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
                        //Depending on the way (true or false) and the last character of the gene name
                        // (+ or -), we add either the first gene (true and + or false and -) or
                        //the last gene (false and + or true and -) to the set of expressed genes
                        let gene;
                        if way {
                            if max_gene.ends_with('+') {
                                //Insert the first element of genes_str.split("; ")
                                gene = genes_str.split("; ").next().unwrap();
                            } else {
                                //Insert the last element of genes_str.split("; ")
                                let l = genes_str.split("; ").collect::<Vec<&str>>();
                                gene = l[l.len()-1];
                            }
                        } else {
                            if max_gene.ends_with('-') {
                                //Insert the first element of genes_str.split("; ")
                                gene = genes_str.split("; ").next().unwrap();
                            } else {
                                //Insert the last element of genes_str.split("; ")
                                let  l = genes_str.split("; ").collect::<Vec<&str>>();
                                gene = l[l.len()-1];
                            }
                        }
                        exp_genes.insert((id,gene,depth,ab,start_way));
                        if old_ab>0 {
                            leaves.insert((id, depth, *ab));
                        }
                        continue; //NE PAS COMMENTER CETTE LIGNE!!!!!!
                    }
                }

                if inter_str != "*" {
                    let l = inter_str.split("; ").collect::<Vec<&str>>();
                    if l.len() == 1 {
                        exp_genes.insert((id,inter_str,depth,ab,start_way));
                        if old_ab>0 {
                            leaves.insert((id, depth, *ab));
                        }
                        continue;
                    }
                }
                //Otherwise, the node is not an expressed gene and we continue the exploration

                //Add the neighbors of the unitig to the stack

                let d:i32 = (seq.len() - K + 1).try_into().unwrap();
                let actual_ab : i32 = nodes.get(&id).unwrap().7;
                for (nb, new_way) in neighbors {
                    let ab_n : i32 = nodes.get(&nb).unwrap().7;
                    //length of nb
                    if (old_ab>0) && (ab_n * 2 >= actual_ab) {
                        //If the abundance of the neighbor is greater than twice the abundance of the current node
                        //we push it to the stack with the new way, depth + d and max_ab
                        stack.push(StackElem(depth+d, nb, new_way,actual_ab, start_way));
                    } else {
                        //Otherwise we push it to the stack with the new way, depth + d and 0
                        stack.push(StackElem(depth+d, nb, new_way,0, start_way));
                    }
                }
            }
            else {
                //Remove the abundance of the node from the flux
                //flux -= old_ab;
            }
        }
        //} //End of the loop over comp_id


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
        //let mut number_intron_genes = 0;
        let mut number_exon_genes = 0;

        file = std::fs::File::create(format!("{}/expressed_genes.txt", output_dir)).unwrap();
        for gene in exp_genes {
            let (seq_id,gene_str,depth,ab,start_way) = gene;
            //Write seq_id and gene_str in the file, separated by a tab

            let start_str = if start_way { "forward" } else { "reverse" };
            writeln!(file, "{}\t{}\t{}\t{}\t{}", seq_id, gene_str,depth,ab,start_str).unwrap();
            number_genes += 1;
            distances_genes_total += depth;
            if gene_str.contains("intron") {
                //number_intron_genes += 1; //never used
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

        let mut nb_leaves = 0;
        file = std::fs::File::create(format!("{}/leaves.txt", output_dir)).unwrap();
        let mut fileS = std::fs::File::create(format!("{}/selected_leaves.txt", output_dir)).unwrap();
        for (leaf,_,old_ab) in &leaves {
            //write the nodes[leaf] in the file
            let (seq, dist, weight, dfam_str, repbase_str, genes_str, inter_str, abundance, ASgen, ASdfam) = nodes.get(&leaf).unwrap();
            if *old_ab > 0 {
                nb_leaves += 1;
                //If the abundance of the leaf is greater than twice the abundance of the node
                //we write the leaf with the depth and old_ab
                writeln!(file, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", leaf, seq, dist, weight, dfam_str, repbase_str, genes_str, inter_str, abundance, ASgen, ASdfam).unwrap();
                writeln!(fileS, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", leaf, seq, dist, weight, dfam_str, repbase_str, genes_str, inter_str, abundance, ASgen, ASdfam).unwrap();

            } else {
                //If the abundance of the leaf is not greater than twice the abundance of the node
                //we write the leaf with the depth and 0
                writeln!(file, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", leaf, seq, dist, weight, dfam_str, repbase_str, genes_str, inter_str, abundance, ASgen, ASdfam).unwrap();
            }
        }
        //close the file
        drop(file);
        drop(fileS);


        //Let compute the number of distinct TE as the size of te hash set
        let mut number_te = 0;
        for el in &te {
            //Check if there is a /t in the string
            //and if so ignore it
            if el.contains("\t") {
                continue;
            }
            else {
                number_te += 1;
            }
        }

        //Let compute the proportion of exon genes
        let proportion_exon_genes = (number_exon_genes as f64 / number_genes as f64) * 100.0;


        //Let loop over the comp_id and store the abundance of the unitigs in a sorted vector
        let mut abundances = Vec::new();
        for idc in &comp_id {
            let (_, _, _, _, _, _, _, ab, _, _) = nodes.get(&idc).unwrap();
            abundances.push(*ab);
        }
        //Sort the vector
        abundances.sort();

        //Now, let loop over the abundances and find the ratios of consecutive abundances (a_{n+1}/a_n)
        //that are greater than 2 and store then in vector ratios

        let mut ratios = Vec::new();
        let  mut max_ratio : f64 = 1.0;

        for i in 0..abundances.len()-1 {
            let ratio = abundances[i+1] as f64 / abundances[i] as f64;
            if ratio >= 2.0 {
                ratios.push(ratio);
                if ratio > max_ratio {
                    max_ratio = ratio;
                }
            }
        }



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
        //print the number of leaves
        writeln!(file, "{} \t {} \t (Leaves conserved and Total number of leaves)", nb_leaves , leaves.len()).unwrap();
        //print the max ratio and the ratios vector
        writeln!(file, "{:.2} \t (Max ratio of abundances greater than 2) \t {:?} \t (Ratios of abundances greater than 2)", max_ratio, ratios).unwrap();
        //close the file
        drop(file);
    } else if args.len() == 7 {
        //If there are 7 arguments, we assume that the component is the one specified by the user
        let nb_comps = args[6].parse::<u32>().unwrap();
        if nb_comps == 0 {
            eprintln!("Error: The number of components must be greater than 0");
            std::process::exit(1);
        }
        //Define a vector of comps_id of size nb_comps
        let mut comp_id_vec = Vec::with_capacity(nb_comps as usize);

        //Define a set of all_comp_id
        let mut all_comp_id = std::collections::HashSet::new();


        //The comps are in file prefix_comps + i + ".txt" for i in [0, nb_comps-1]
        for i in 0..nb_comps {
            let file_name = format!("{}{}.txt", &args[1], i);
            let file = std::fs::File::open(&file_name).unwrap();
            let reader = std::io::BufReader::new(file);
            let mut comp_id = std::collections::HashSet::new();
            for line in reader.lines() {
                let id = line.unwrap().split('\t').next().unwrap().parse::<u32>().unwrap();
                comp_id.insert(id);
                all_comp_id.insert(id);
            }
            comp_id_vec.push(comp_id);
        }

        //Define max weight and max_id empty vector
        let mut max_weight = Vec::new();
        let mut max_id = Vec::new();
        //Loop over the comp_id_vec and compute the max weight and max_id for each component
        for comp_id in &comp_id_vec {
            let mut max_w = 0;
            let mut max_i = 0;
            for idc in comp_id {
                let (_, _, w, _, _, _, _, _, _, _) = nodes.get(&idc).unwrap();
                if *w > max_w {
                    max_w = *w;
                    max_i = *idc;
                }
            }
            max_weight.push(max_w);
            max_id.push(max_i);
        }
        //Set of the expressed transcript
        let mut exp_transcript = std::collections::HashSet::new();

        //Set of unitigs connecting components
        let mut connecting_unitigs = std::collections::HashSet::new();

        #[derive(Eq, PartialEq)]
        struct StackElem(i32, u32, bool, i32, bool, usize); // depth, id, way, old_ab, start_way, initial_comp

        impl Ord for StackElem {
            fn cmp(&self, other: &Self) -> std::cmp::Ordering { //True when self.0 is less than other.0
                // Decreasing order by the first element (depth)
                other.0.cmp(&self.0) //True when other.0 is greater than self.0
                    .then_with(|| self.1.cmp(&other.1)) // For tie-breaking, compare by id
            }
        }
        impl PartialOrd for StackElem {
            fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
                Some(self.cmp(other))
            }
        }

        let mut stack = std::collections::BinaryHeap::new();

        //Add the heaviest unitig of each component to the stack
        for (i, comp_id) in comp_id_vec.iter().enumerate() {
            if !comp_id.is_empty() {
                stack.push(StackElem(0, max_id[i], true, 0, true, i.try_into().unwrap()));
                stack.push(StackElem(0, max_id[i], false, 0, false, i.try_into().unwrap()));
            }
        }

        //Set of the unitigs that are already visited
        let mut visited = std::collections::HashSet::new();

        //Vector of visited nodes for each component
        let mut visited_in_comp = Vec::new();
        for _ in 0..nb_comps {
            visited_in_comp.push(std::collections::HashSet::new());
        }

        //Dictionary of the comp where the unitig have been visited : visited_loc[id] = Vec<i32>
        let mut visited_loc = std::collections::HashMap::new();


        //Define the flux for each component
        let mut flux = Vec::new();
        let mut leaves = Vec::new();
        for _ in 0..nb_comps {
            flux.push(0);
            leaves.push(0);
        }


        //Looping over the stack
        while !stack.is_empty() {
            //Pop the first element of the stack and its way
            let StackElem(d, id, way, o_a, start_way, i) = stack.pop().unwrap();
            let mut depth = d;
            let mut old_ab = o_a;
            let actual_ab: i32 = nodes.get(&id).unwrap().7;
            //If the unitig is not visited
            //Check if the unitig is in the component
            let seen_id_way = visited.contains(&(id, way));
            let seen_id_in_comp = visited_in_comp[i].contains(&id);
            let seen_id_other_way = visited.contains(&(id, !way));

            if !seen_id_way && (seen_id_in_comp || !seen_id_other_way){
                visited.insert((id, way)); //Add the unitig to the visited set
                visited_in_comp[i].insert(id); //Add the unitig to the visited_in_comp set
                visited_loc.entry(id).or_insert_with(Vec::new).push(i); //Add the unitig to the visited_loc set

                //If the unitig is in the component
                if comp_id_vec[i].contains(&id) {
                    depth = 0;
                    old_ab = nodes.get(&id).unwrap().7; //Get the abundance of the unitig
                }

                //A valid neighbor is an adjacent node such that its edge type is FR/FF if way == true
                //or RF/RR if way == false
                let mut neighbors = Vec::new();

                if edges.contains_key(&id) {
                    for (nb, edge_type) in edges.get(&id).unwrap() {
                        if (way && (edge_type == "FR" || edge_type == "FF")) || (!way && (edge_type == "RF" || edge_type == "RR")) {
                            let new_way: bool;
                            if edge_type == "FF" || edge_type == "RF" {
                                new_way = true;
                            } else {
                                new_way = false;
                            }
                            neighbors.push((*nb, new_way));
                        }
                    }
                }
                //If the unitig is a leaf
                if neighbors.len() == 0 {
                    exp_transcript.insert((id, depth, old_ab,i));
                    flux[i] += actual_ab;
                    if old_ab > 0 {
                        leaves[i]+= 1;
                    }
                }
                //Add the neighbors of the unitig to the stack

                let d : i32 = (nodes.get(&id).unwrap().0.len() - K + 1).try_into().unwrap();
                let mut max_ab = 0;
                let mut best_el = 0;
                for (nb, _) in &neighbors {
                    let ab_n: i32 = nodes.get(&nb).unwrap().7;
                    if ab_n > max_ab {
                        max_ab = ab_n;
                        best_el = *nb;
                    }
                }
                for (nb, new_way) in neighbors {
                    //stack.push(StackElem(depth + d, nb, new_way,  actual_ab, start_way, i.try_into().unwrap()));
                    if (old_ab > 0) && (nb == best_el) {
                        //If the abundance of the neighbor is greater than twice the abundance of the current node
                        //we push it to the stack with the new way, depth + d and max_ab
                        stack.push(StackElem(depth + d, nb, new_way,  actual_ab, start_way, i.try_into().unwrap()));
                    } else {
                        //Otherwise we push it to the stack with the new way, depth + d and 0
                        stack.push(StackElem(depth + d, nb, new_way, 0, start_way, i.try_into().unwrap()));
                    }
                }
            }
            if seen_id_way {
                //Case where the node has already been visited

                //First case : it has been visited from another component
                if !seen_id_in_comp {
                    //Add the node to the visited_in_comp set
                    visited_in_comp[i].insert(id);
                    visited_loc.entry(id).or_insert_with(Vec::new).push(i);
                    //Add the node to the connecting_unitigs set
                    connecting_unitigs.insert(id);
                    //flux[i] += old_ab;
                } else {
                    //Second case : it has been visited from the same component
                    //We have a bubble
                }
            }

            if !seen_id_in_comp && seen_id_other_way { //Id seen in the other way and from another component
                //If the node has not been visited in the current component
                visited_loc.entry(id).or_insert_with(Vec::new).push(i);
                //Add the node to the connecting_unitigs set
                connecting_unitigs.insert(id);
                //flux[i] += old_ab;
            }
        }

        //Compute for each component the number of other components it is connected to
        //by connecting_unitigs
        let mut comp_connections = vec![std::collections::HashSet::new(); nb_comps as usize];
        for id in &connecting_unitigs {
            let comp_ids = visited_loc.get(id).unwrap();
            //Add the pairwise combinations of comps in comp_ids to the comp_connections set
            for i in 0..comp_ids.len() {
                for j in (i + 1)..comp_ids.len() {
                    let comp_i = comp_ids[i];
                    let comp_j = comp_ids[j];
                    if comp_i != comp_j {
                        comp_connections[comp_i as usize].insert(comp_j);
                        comp_connections[comp_j as usize].insert(comp_i);
                    }
                }
            }
        }


        //output the results
        //Output the expressed transcript to a file

        let output_dir = &args[4];
        let mut file = std::fs::File::create(format!("{}/expressed_transcripts.txt", output_dir)).unwrap();
        for (id, depth, abundance, i) in exp_transcript {
            //Write the id, depth and abundance in the file, separated by a tab
            let seq = nodes.get(&id).unwrap().0.clone();
            writeln!(file, "{}\t{}\t{}\t{}\t{}", id, i, seq, depth, abundance).unwrap();
        }
        //close the file
        drop(file);

        //Output the connecting unitigs to a file
        file = std::fs::File::create(format!("{}/connecting_unitigs.txt", output_dir)).unwrap();
        for id in &connecting_unitigs {
            //Write the id in the file
            let seq = nodes.get(&id).unwrap().0.clone();
            let comp_ids = visited_loc.get(id).unwrap();
            let comp_ids_str = comp_ids.iter().map(|x| x.to_string()).collect::<Vec<String>>().join(",");
            writeln!(file, "{}\t{}\t{}", id, seq, comp_ids_str).unwrap();
        }
        //close the file
        drop(file);

        //transcript_summary_comps.txt
        file = std::fs::File::create(format!("{}/transcript_summary_comps.txt", output_dir)).unwrap();
        //Write the number of expressed transcripts for each component
        for (i, comp_id) in comp_id_vec.iter().enumerate() {
            let mut max_ab = 0;
            for id in comp_id {
                let ab = nodes.get(&id).unwrap().7;
                if ab > max_ab {
                    max_ab = ab;
                }
            }
            writeln!(file, "Component {}: {} flux \t {} max ab \t {} leaves \t {} connected comp", i, flux[i], max_ab, leaves[i], comp_connections[i].len()).unwrap();
        }
        //close the file
        drop(file);
    }
}
