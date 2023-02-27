set_pc_input <- function(input_matrix, target_gene, exp_list, cut, organism){ 
  
  if(organism == "Hs"){
    target_tid <- strsplit(target_gene, "-")[[1]][1]
    target_name <- strsplit(target_gene, "-")[[1]][2]
  
    if(cut == round(as.numeric(cut))){
      #cut the list according to number of genes
      exp_list <- exp_list[as.numeric(exp_list$rank) <= as.numeric(cut),]
      exp_list <- exp_list[!(exp_list$gene_name == target_name),]
      exp_list <- exp_list[!duplicated(exp_list$gene_name), ]
      tid_list <- exp_list$ID[as.numeric(exp_list$rank) <= as.numeric(cut)]
    
    } else {
      #cut the list according to relative frequency
      exp_list <- exp_list[as.numeric(exp_list$Frel) >= as.numeric(cut),]
      exp_list <- exp_list[!(exp_list$gene_name == target_name),] 
      exp_list <- exp_list[!duplicated(exp_list$gene_name), ]
      tid_list <- exp_list$ID[as.numeric(exp_list$Frel) >= as.numeric(cut)]
    }
    tid_list <- append(tid_list, target_tid)
  
  }else if (organism == "Vv") {
    if(cut == round(as.numeric(cut))){
      #cut the list according to number of genes
      tid_list <- exp_list$ID[as.numeric(exp_list$rank) <= as.numeric(cut)]
    } else {
      #cut the list according to relative frequency
      tid_list <- exp_list$ID[as.numeric(exp_list$Frel) >= as.numeric(cut)]
    }
    tid_list <- append(tid_list, target_gene) #, after=0
  }
  print("How many nodes are given as input to pc()")
  print(length(tid_list))
  subset_matrix <- data.frame()
  for (i in tid_list) {
    for (j in row.names(input_matrix)) {
      if (i == j) {
        subset_matrix <- rbind(subset_matrix, input_matrix[j, ])
      }
    }
  }
  print(subset_matrix[, c(1:2)])
  #subset_matrix <- subset_matrix[nrow(subset_matrix):1,]
  print(subset_matrix[, c(1:2)])
  t_subset_matrix <- t(subset_matrix)
  
  #n_tids <- paste("/", as.character(nrow(subset_matrix)),sep = "", collapse = NULL)#
  #tids_dir <- paste(path, n_tids, sep = "", collapse = NULL)#
  #tids_dir <- paste(tids_dir, "_tids.csv", sep = "", collapse = NULL)#
  #write.csv(tid_list, tids_dir)#
  
  return(t_subset_matrix)
}

load_matrix <- function(matrix_dir){
  input_matrix <- read.csv(matrix_dir, header = TRUE, sep = ",", dec = ".", row.names = 1)
  print(dim(input_matrix))
  return(input_matrix)
}

apply_pc_parallel <- function (pc_input_matrix, pc_output_dir){
  
  if(require(pcalg) && require(parallel) && require(ParallelPC)){
    n <- nrow(pc_input_matrix)
    v <- colnames(pc_input_matrix)
    
    start_time <- Sys.time()
    cat("\n pc() started at", format(start_time, "%a %b %d %X %Y"))
    
    pc_fit <-pc_parallel(
      suffStat = list(C = cor(pc_input_matrix), n = n),
      indepTest = gaussCItest, 
      alpha = 0.05,
      labels = v,
      u2pd = "relaxed",
      skel.method = "parallel",
      mem.efficient = TRUE, 
      conservative = FALSE,
      maj.rule = TRUE,
      solve.confl = TRUE,
      verbose = FALSE
    )
    
    end_time <- Sys.time()
    cat("\n pc() ended at", format(end_time, "%a %b %d %X %Y"))
    
    sink(pc_output_dir)
    showEdgeList(pc_fit, labels = NULL)
    sink()
  }
}

apply_pc <- function (pc_input_matrix, pc_output_dir){
  
  if(require(pcalg)){
    n <- nrow(pc_input_matrix)
    v <- colnames(pc_input_matrix)
    
    start_time <- Sys.time()
    cat("\n pc() started at", format(start_time, "%a %b %d %X %Y"))
    
    pc_fit <-pc(
      suffStat = list(C = cor(pc_input_matrix), n = n),
      indepTest = gaussCItest, 
      alpha = 0.05,
      labels = v,
      u2pd = "relaxed",
      skel.method = "stable",
      conservative = FALSE,
      maj.rule = TRUE,
      solve.confl = TRUE,
      verbose = FALSE
    )
    
    end_time <- Sys.time()
    cat("\n pc() ended at", format(end_time, "%a %b %d %X %Y"))
    
    sink(pc_output_dir)
    showEdgeList(pc_fit, labels = NULL)
    sink()
  }
}

txt_to_df <- function (txt_file){
  data <- readLines(txt_file)
  undir_df <- data.frame(source = character(),
                         interaction = character(),
                         target = character())
  i <- 5
  while (data[i] != "") {
    j <- unlist(strsplit(data[i], "  "))
    df_row <- data.frame(source = j[2], interaction = j[3], target = trimws(j[4]))
    undir_df <- rbind(undir_df, df_row)
    i <- i + 1
  }
  if(nrow(undir_df)==0){
    #print("undir_df is empty")
  } else{
  }
  i <- i + 2
  dir_df <- data.frame(source = character(),
                       interaction = character(),
                       target = character())
  while (!is.na(data[i])) {
    j <- unlist(strsplit(data[i], "  "))
    df_row <- data.frame(source = j[2], interaction = j[3], target = trimws(j[4]))
    dir_df <- rbind(dir_df, df_row)
    i <- i + 1
  }
  if(nrow(dir_df)==0){
    #print("dir_df is empty")
  } else{
  }
  cpdag <- rbind(undir_df, dir_df)
  row.names(cpdag) <- 1:nrow(cpdag)
  if(nrow(cpdag)==0){
    #print("cpdag is empty")
  } else{
    #print("dataframe: all edges: ")
    #print(cpdag)
  }
  
  return(cpdag)
}

get_nodes <- function (cpdag, expansion_list){
  lazy_nodes_list <-c(cpdag[["source"]], cpdag[["target"]])
  unique_nodes_list <- unique(lazy_nodes_list)
  nodes_list <- data.frame()
  for (node in unique_nodes_list) {
    nodes_list <- rbind(nodes_list, expansion_list[expansion_list$ID == node,])
  }
  row.names(nodes_list) <- 1:nrow(nodes_list)
  return(nodes_list)
}

get_pearson_corr <- function (cpdag, pc_input_matrix){
  for (i in 1:nrow(cpdag)) {
    x <- cpdag[i, "source"]
    y <- cpdag[i, "target"]
    pearson_corr <- cor(pc_input_matrix[,x],pc_input_matrix[,y])
    cpdag[i, "cor"] <- pearson_corr
  }
  cpdag[cpdag$cor > 0, "cor_sign"] <- "+"
  cpdag[cpdag$cor <= 0, "cor_sign"] <- "-"
  return(cpdag)
}

format_Hs <- function(nodes_list){
  i <- 0
  for (i in 1:nrow(nodes_list)){
    if(grepl(",", nodes_list$association_with_transcript[i])){
      nodes_list$association_with_transcript[i]<- gsub(",",";", nodes_list$association_with_transcript[i])
    }
    if(grepl(",", nodes_list$entrezgene_id[i])){
      nodes_list$entrezgene_id[i]<- gsub(",",";", nodes_list$entrezgene_id[i])
    }
    if(grepl(",", nodes_list$hgnc_id[i])){
      nodes_list$hgnc_id[i]<- gsub(",",";", nodes_list$hgnc_id[i])
    }
    if(grepl(",", nodes_list$uniprot_id[i])){
      nodes_list$uniprot_id[i]<- gsub(",",";", nodes_list$uniprot_id[i])
    }
    if(grepl(",", nodes_list$description[i])){
      nodes_list$description[i]<- gsub(",",";", nodes_list$description[i])
    }
    if(grepl(",", nodes_list$type[i])){
      nodes_list$type[i]<- gsub(",",";", nodes_list$type[i])
    }
  }
  return(nodes_list)
}

check_deleted_nodes <- function(nodes_list, input_nodes, expansion_list) {
  if (length(setdiff(input_nodes, nodes_list[,"ID"])) == 0){ 
    print("All input nodes belong to the output graph")
  } else {
    print("The following input nodes are isolated in the output graph: ")
    #print(setdiff(input_nodes, nodes_list[,"ID"]))
    del_nodes <- setdiff(input_nodes, nodes_list[,"ID"])
    del_nodes_list <- list()
    for (i in del_nodes) {
      for (j in expansion_list$ID) {
        if(i == j){
          del_nodes_list <- append(del_nodes_list, expansion_list$gene_name[expansion_list$ID==i])
        }
      }
    }
    print(unlist(del_nodes_list))
  }
}