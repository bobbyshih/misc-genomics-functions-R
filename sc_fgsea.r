# comparison_ident must be a character vector containing comparison groups
# control_ident must be a character vector containing a single control group

# fgsea_sets must be a named list containing gene set name and list of genes
# It is recommended that you create a df (msigdbr_df) containing gene sets of interest using the msigdbr package and the split function to assemble into a list
# fgsea_sets <- msigdbr_df %>% split(x = .$gene_symbol, f = .$gs_name)

# It is recommended that you create a separate ident class within the seurat object with easily understandable groups

flatten_list = function(x){
    if (typeof(x) != "list") {
        return(x)
    }
    sapply(x, function(y) paste(y, collapse = ", "))
}

sc_fgsea = function(seurat_obj, comparison_ident, control_ident, logfc.threshold = 0, min.pct = 0.1, slot = "data", test.use = "MAST", fgsea_sets = fgsea_sets, minSize = 15, maxSize = 500) {
    for (i in 1:length(comparison_ident)) {
      Seurat_DESeq <- FindMarkers(seurat_obj,
                                  ident.1 = comparison_ident[i],
                                  ident.2 = control_ident,
                                  slot = slot,
                                  test.use = test.use,
                                  logfc.threshold = logfc.threshold,
                                  min.pct = min.pct)

      write.csv(Seurat_DESeq,
                file = paste(comparison_ident[i], "_vs_", control_ident, "_DE", test.use, ".csv", sep = ""),
                row.names = TRUE)

      # Organize data for fgsea input
      ranks <- data.frame(matrix(ncol = 2, nrow = nrow(Seurat_DESeq)))
      ranks[,1] <- rownames(Seurat_DESeq)
      ranks[,2] <- Seurat_DESeq %>% arrange(desc(avg_log2FC)) %>% dplyr::select(avg_log2FC)
      ranks <- deframe(ranks)

      # Run fgsea
      fgseaRes <- fgsea(fgsea_sets, 
                        stats = ranks,
                        minSize = minSize,
                        maxSize = maxSize)
        
        fgseaRes %>% mutate(across("leadingEdge", flatten_list)) %>% write.csv(
            file = paste(comparison_ident[i], "_vs_", control_ident, "_fgseaRes.csv", sep = ""), 
                                                                           row.names = TRUE)
    
        saveRDS(fgseaRes, 
          file = paste(comparison_ident[i], "_vs_", control_ident, "_fgseaRes.rds", sep = ""))
    }
}