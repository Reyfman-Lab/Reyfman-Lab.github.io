install_wflow_packages <- function(cran_pkg = NULL,
                                   bioconductor_pkg = NULL,
                                   github_pkg = NULL,
                                   update_bioc = FALSE) {
  if (!is.null(cran_pkg)) {
    cran_pkg_install <-
      cran_pkg[!cran_pkg %in% installed.packages()[, "Package"]]

    if (!identical(cran_pkg_install, "character(0)")) {
      install.packages(cran_pkg_install)
    }
  }

  if (!is.null(bioconductor_pkg)) {
    if((!"BiocManager" %in% installed.packages()[, "Package"])){
      install.packages("BiocManager")
    }
    bc_pkg_install <-
      bioconductor_pkg[!bioconductor_pkg %in% installed.packages()[, "Package"]]
    if (!identical(bc_pkg_install, "character(0)")) {

    }
    BiocManager::install(bc_pkg_install, update = update_bioc)
  }

  if (!is.null(github_pkg)) {
    if((!"devtools" %in% installed.packages()[, "Package"])){
      install.packages("devtools")
    }
    gh_pkg_install <-
      github_pkg[!github_pkg %in% installed.packages()[, "Package"]]
    if (!identical(gh_pkg_install, "character(0)")) {

    }
    devtools::install_github(gh_pkg_install)
  }

  invisible(sapply(c(cran_pkg, bioconductor_pkg, github_pkg), function(x)
    require(x, character.only = TRUE)))
}





plot_pca <- function(pilot_dge, genotype, title){

  pilot_pca <- pilot_dge %>%
    edgeR::cpm(log = TRUE) %>%
    t() %>%
    prcomp(center=TRUE, scale=FALSE) %>%
    extract2("x") %>%
    as_tibble(rownames = "Samples") %>%
    mutate(Genotype = genotype)

  pilot.pca.var <- pilot_pca %>%
    summarize(across(.cols = -c(Samples, Genotype), ~var(.x))) %>%
    mutate(Total = rowSums(.)) %>%
    mutate(across(.cols = -Total, ~(.x/Total) * 100)) %>%
    mutate(across(.cols = -Total, ~round(.x)))

  ggplot(pilot_pca, mapping = aes(x = PC1, y = PC2, color = Genotype)) +
    geom_point() +
    xlab(glue("PC1 ({pilot.pca.var$PC1}% of Variance)")) +
    ylab(glue("PC2 ({pilot.pca.var$PC2}% of Variance)")) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme_bw()
}

cpmGeneExpr <- function(counts, gene.col, goi, log = TRUE) {
  counts <- tibble::column_to_rownames(counts, var = gene.col)
  if (isTRUE(log)){
    counts <- edgeR::cpm(counts, log = TRUE)
  } else {
    counts <- edgeR::cpm(counts, log = FALSE)
  }
  counts <- tibble::as_tibble(counts, rownames = gene.col)
  counts <- dplyr::filter(counts, !!as.name(gene.col) %in% goi)
  counts <- tidyr::pivot_longer(counts,cols = -all_of(gene.col), names_to = "Samples",
                                values_to = "CPM")
  return(counts)
}


mrna_cpm_plot <- function(counts_df, gene) {
  counts_df <- dplyr::filter(counts_df, Gene_Name == gene)
  p1 <- ggplot(data = counts_df, mapping = aes(x = Samples, y = CPM, fill = Samples)) +
    geom_boxplot() +
    theme(legend.position = "none") +
    theme(panel.background = element_rect(fill = "NA")) +
    theme(panel.border = element_rect(fill = "NA", size = 2)) +
    scale_fill_manual(values = c("#00ba38","#f8766d")) +
    ylab("logCPM") +
    theme(text = element_text(size = 30)) +
    geom_point(size = 4) +
    ggtitle(gene)
}

plot_volcano <- function(top_tags) {

  top_tags <- top_tags %>%
    .$table %>%
    as_tibble(rownames = "Genes") %>%
    mutate(Significance = p.adjust(PValue, method = "BH") < .5,
           Significance = Significance * sign(logFC),
           Significance = factor(Significance),
           logP = -log10(PValue))

  ggplot(data = top_tags, aes(logFC, logP)) +
    geom_point(aes(col = Significance)) +
    scale_color_manual(name="FDR < 0.05",values=c("black","red","blue"),breaks=c(0,1,-1),
                       labels=c("non-DE","Up","Down")) +
    geom_label_repel(data = top_tags %>% slice_min(FDR, n = 10, with_ties = FALSE), size=7, aes(label = Genes))+
    theme_pubr()+
    theme(plot.title = element_text(hjust = 0.5))+
    ylab("-log(P-value)")+
    xlab("log2(Fold Change)") +
    ggtitle("Volcano Plot: F508del vs Control") +
    annotate(geom = "text",x = 2,
             y = 9,
             label = glue("Up: {top_tags %>% filter(FDR <= .05 & logFC > 0) %>% tally()}"),
             color = "red",
             size = 5) +
    annotate(geom = "text",
             x = -2,
             y = 9,
             label = glue("Down: {top_tags %>% filter(FDR <= .05 & logFC < 0) %>% tally()}"),
             color = "blue",
             size = 5)
}

plot_heatmap <- function(cpm_counts, group_annot, study){
  pheatmap(mat = cpm_counts,
           color = plasma(20),
           annotation_col = group_annot,
           cluster_cols = F,
           clustering_distance_rows = "correlation",
           clustering_distance_cols = "correlation",
           show_rownames = F,
           main = glue("Heatmap of {study} DEGS: F508del vs Control"),
           show_colnames = F,
           fontsize = 10,
           border_color = NA)
}

figure_save <- function(figure, name,
                        add_date = TRUE,
                        fig.height = 20,
                        fig.width = 20,
                        device = "png" ){

  if(isTRUE(add_date)){
    sys_date <- str_replace_all(Sys.Date(), "-","")
    name <- str_c(sys_date, name, sep = "_")
    name <- str_c(name, device, sep = ".")
  }

  ggsave(filename = name,
         plot = figure,
         width = fig.width,
         height = fig.height,
         device = device)
}

genesToEntrez <- function(data, ...){
  UseMethod("genesToEntrez")
}

genesToEntrez.tbl_df <- function(data,
                                 column,
                                 drop.na = FALSE) {


  ncbi_ids <- read_tsv("data/ncbi_genes.txt", col_names = TRUE, col_types = list(col_character(), col_character()))
  colnames(ncbi_ids[2]) <- column
  data <- dplyr::right_join(ncbi_ids, data, by = "Gene_Name")

  if (isTRUE(drop.na)){
    data <- tidyr::drop_na(data, NCBI_ID)
  }

  return(data)
}

genesToEntrez.DGEList<- function(data, drop.na = TRUE, return_DGE_Obj = FALSE){
  gene_names <- data %>% rownames(.) %>%
    as_tibble() %>%
    dplyr::rename(Gene_Name = value) %>%
    genesToEntrez(., column = "Gene_Name", drop.na = drop.na)

  if(isTRUE(return_DGE_Obj)) {
    data <- data[gene_names$Gene_Name, , keep.lib.sizes = TRUE]
    data$Entrez <- gene_names$NCBI_ID
    return(data)
  }

  return(gene_names)

}

genesToEntrez.DGELRT<- function(data, drop.na = TRUE, return_DGE_Obj = FALSE){
  gene_names <- data %>% rownames(.) %>%
    as_tibble() %>%
    dplyr::rename(Gene_Name = value) %>%
    genesToEntrez(., column = "Gene_Name", drop.na = drop.na)

  if(isTRUE(return_DGE_Obj)) {
    data <- data[gene_names$Gene_Name, ]
    data$Entrez <- gene_names$NCBI_ID
    return(data)
  }

  return(gene_names)

}

