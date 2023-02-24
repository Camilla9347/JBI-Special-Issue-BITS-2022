#R version 4.1.1 (2021-08-10) -- "Kick Things"
#Copyright (C) 2021 The R Foundation for Statistical Computing
#Platform: aarch64-apple-darwin20 (64-bit)
require(readr)
source("procedure_functions.R")
options(max.print = .Machine$integer.max)
#options(warn=-1)
options("digits"=15)
path <- getwd()

args = commandArgs(trailingOnly=TRUE)


if (length(args)==0 | length(args)==1 | length(args)==2) {
  stop("At least 3 argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==3) {
  
  
  expansion_list_dir <- c(paste(path, paste("/", args[1], sep = "", collapse = NULL) , sep = "", collapse = NULL))
  pc_output_dir <- paste(path, "/pc_output.txt", sep = "", collapse = NULL)
  
  # Hs EXPANSION LIST PROCESSING
  
  if(args[2] == "Hs"){
    
    # Hs EXPANSION LIST LOADING
    expansion_list <- read.csv(expansion_list_dir, header = TRUE, sep = ",", dec = ".", skip = 1)
    
    # FANTOM MATRIX LOADING
    print("Loading fantom matrix...")
    matrix_dir <- paste(path, "/fantom_mat.csv", sep = "", collapse = NULL)
    input_matrix <- load_matrix(matrix_dir)
    print("Loaded fantom matrix.")
    
    # PC INPUT SETTING
    header <- read.table(args[1], header = F, sep = ",", nrows = 1)
    target_gene <- header[1, "V4"]
    target_isoform <- header[1, "V5"]
    target_tid <- strsplit(target_gene, "-")[[1]][1]
    pc_input_matrix <- set_pc_input(input_matrix, target_gene, expansion_list, args[3], args[2])
    
    # PC APPLICATION
    
    apply_pc_parallel(pc_input_matrix, pc_output_dir)
    
    # PC OUTPUT PROCESSING FOR VISUALIZATION
    cpdag <- txt_to_df(pc_output_dir)
    file.remove(pc_output_dir)
    print("Loading Hs annotation information")
    fantom_anno_table <- read.csv("anno-hsf5.csv", header=FALSE, sep = ",", dec = ".")
    print("Loaded Hs annotation information")
    gene_of_interest <- data.frame(
      rank = "0",
      ID = fantom_anno_table$V1[fantom_anno_table$V1 == target_tid],
      Fabs = NA, 
      Frel = 1,
      association_with_transcript = fantom_anno_table$V2[fantom_anno_table$V1 == target_tid],
      entrezgene_id = fantom_anno_table$V3[fantom_anno_table$V1 == target_tid], 
      hgnc_id = fantom_anno_table$V4[fantom_anno_table$V1 == target_tid],
      uniprot_id = fantom_anno_table$V5[fantom_anno_table$V1 == target_tid],
      gene_name = fantom_anno_table$V6[fantom_anno_table$V1 == target_tid],
      description = fantom_anno_table$V7[fantom_anno_table$V1 == target_tid], 
      type = fantom_anno_table$V8[fantom_anno_table$V1 == target_tid]
    )
    expansion_list <- rbind(gene_of_interest, expansion_list)
    nodes_list <- get_nodes(cpdag, expansion_list)
    check_deleted_nodes(nodes_list, colnames(pc_input_matrix), expansion_list)
    nodes_list <- format_Hs(nodes_list)
    nodes_list <- nodes_list[, c(2,5,9,6,7,8,10,1,4,11)]
    nodes_list[is.na(nodes_list)] <- ""
    cpdag <- get_pearson_corr (cpdag, pc_input_matrix) # INITIAL PEARSON CORRELATION (EDGES)
    nodes_dir <- paste(path, paste(target_isoform, "_nodes.csv", sep = "", collapse = NULL) , sep = "/Hs/", collapse = NULL)
    write.csv(nodes_list, nodes_dir, row.names=FALSE, quote=FALSE)
    cpdag_dir <- paste(path, paste(target_isoform, "_edges.csv", sep = "", collapse = NULL) , sep = "/Hs/", collapse = NULL)
    write.csv(cpdag, cpdag_dir, row.names=FALSE, quote=FALSE)
    
    #Vv EXPANSION LIST PROCESSING
    
  } else if (args[2] == "Vv") {
    
    #Vv EXPANSION LIST LOADING
    #expansion_list <- read.csv(expansion_list_dir, header = TRUE, sep = ",", dec = ".", skip = 1)
    expansion_list <- read_csv(expansion_list_dir, skip = 1, progress =FALSE, show_col_types = FALSE)
    
    colnames(expansion_list)[2] <- "ID"
    colnames(expansion_list)[11] <- "gene_name_synonyms"
    expansion_list[ , 16:17] <- list(NULL)
    
    # VESPUCCI MATRIX LOADING
    print("Loading vespucci matrix...")
    matrix_dir <- paste(path, "/vespucci_mat.csv", sep = "", collapse = NULL)
    input_matrix <- load_matrix(matrix_dir)
    print("Loaded vespucci matrix.")
    
    # PC INPUT SETTING
    target_gene <- read.table(args[1], header = F, sep = ",", nrows = 1)
    target_vit <- target_gene[1, "V4"]
    pc_input_matrix <- set_pc_input(input_matrix, target_vit, expansion_list, args[3], args[2])
    
    # PC APPLICATION
    
    apply_pc_parallel(pc_input_matrix, pc_output_dir)
    
    # PC OUTPUT PROCESSING (VISUALISATION)
    cpdag <- txt_to_df(pc_output_dir)
    file.remove(pc_output_dir)
    print("Loading Vv annotation information")
    vespucci_anno_table_dir <- paste(path, "/anno-vvv2.csv", sep = "", collapse = NULL)
    vespucci_anno_table <- read.csv(vespucci_anno_table_dir, header=FALSE, sep = ",", dec = ".", skip=1)
    print("Loaded Vv annotation information")
    
    gene_of_interest <- data.frame(
      rank = "0",
      ID = vespucci_anno_table$V2[vespucci_anno_table$V2 == target_vit],
      Fabs = NA, 
      Frel = 1,
      gene_symbol = vespucci_anno_table$V3[vespucci_anno_table$V2 == target_vit],
      v3 = vespucci_anno_table$V4[vespucci_anno_table$V2 == target_vit], 
      group = vespucci_anno_table$V5[vespucci_anno_table$V2 == target_vit], 
      subgroup = vespucci_anno_table$V6[vespucci_anno_table$V2 == target_vit], 
      pathway_or_family = vespucci_anno_table$V7[vespucci_anno_table$V2 == target_vit],
      full_name = vespucci_anno_table$V8[vespucci_anno_table$V2 == target_vit],
      gene_name_synonyms = vespucci_anno_table$V9[vespucci_anno_table$V2 == target_vit],
      type = vespucci_anno_table$V10[vespucci_anno_table$V2 == target_vit],
      reaction_EC = vespucci_anno_table$V11[vespucci_anno_table$V2 == target_vit],
      Vitisnet_Network1 = vespucci_anno_table$V12[vespucci_anno_table$V2 == target_vit],
      Vitisnet_Network2 = vespucci_anno_table$V13[vespucci_anno_table$V2 == target_vit]
    )
    expansion_list <- rbind(gene_of_interest, expansion_list)
    nodes_list <- get_nodes(cpdag, expansion_list)
    check_deleted_nodes(nodes_list, colnames(pc_input_matrix), expansion_list)
    nodes_list <- nodes_list[, c(2,5,6,7,8,9,10,11,12,13, 14,15,1,4)]
    nodes_list[is.na(nodes_list)] <- ""
    cpdag <- get_pearson_corr (cpdag, pc_input_matrix) # INITIAL PEARSON CORRELATION (EDGES)
    nodes_dir <- paste(path, paste(target_vit, "_nodes.csv", sep = "", collapse = NULL) , sep = "/Vv/", collapse = NULL)
    write.csv(nodes_list, nodes_dir, row.names=FALSE, quote=FALSE)
    cpdag_dir <- paste(path, paste(target_vit, "_edges.csv", sep = "", collapse = NULL) , sep = "/Vv/", collapse = NULL)
    write.csv(cpdag, cpdag_dir, row.names=FALSE, quote=FALSE)
  } 
}




