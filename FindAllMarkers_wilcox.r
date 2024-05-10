avg_log2 <- function(expr_array) {
    log2FC <- log(mean(expr_array), 2)
    return(log2FC)
}

FindAllMarkers_wilcox <- function(expr_matrix, annot_matrix, group_column = NULL, ident.1 = NULL, ident.2 = NULL) {
    # Perform a Wilcoxon rank-sum test as differential expression
    # expr_matrix should be a gene (row) x sample (column) data frame
    # annot_matrix should be a sample (row) x annotation (column) data frame
    # group_column specifies the annotation column in annot_matrix to be used
    # If only ident.1 is specified, will perform a comparison of ident.1 versus all other samples
    # If ident.1 and ident.2 is specified, will perform a comparison between the two

    suppressPackageStartupMessages({
        library(dplyr)
    })

    empty_vars <- c()
    if (is.null(expr_matrix)) {
    empty_vars <- c(empty_vars, "expr_matrix")
    }
    
    if (is.null(annot_matrix)) {
        empty_vars <- c(empty_vars, "annot_matrix")
    }
    
    if (is.null(group_column)) {
        empty_vars <- c(empty_vars, "group_column")
    }

    if (is.null(ident.1)) {
        empty_vars <- c(empty_vars, "ident.1")
    }

    if (length(empty_vars) > 0) {
        warning_message <- paste("Warning: The following variables are empty:", paste(empty_vars, collapse = ", "))
        return(warning_message)
    } else if (!(ident.1 %in% annot_matrix[, group_column])) {
        warning_message <- paste(ident.1, " is not found in ", group_column, sep = "")
        return(warning_message)
    } else if (!is.null(ident.2) & !(ident.2 %in% annot_matrix[, group_column])) {
        warning_message <- paste(ident.2, " is not found in ", group_column, sep = "")
        return(warning_message)
    } else if (!identical(colnames(expr_matrix), rownames(annot_matrix))) {
        warning_message <- "expr_matrix and annot_matrix sample indexes do not match"
        return(warning_message)
    }

    # Initialize a result data frame
    DE_result_df <- data.frame(gene = rownames(expr_matrix),
                                avg_log2FC = numeric(nrow(expr_matrix)),
                                p_value = numeric(nrow(expr_matrix)),
                                stringsAsFactors = TRUE)
    
    # Perform Wilcox Rank-Sum test over each gene
    if (is.null(ident.2)) { # Perform a comparison of ident.1 against all other samples
        for (i in 1:nrow(expr_matrix)) {
            gene_expr <- expr_matrix[i,]

            # Split gene expression values based on group labels
            mask <- annot_matrix[, group_column] == ident.1
            group_1 <- gene_expr[mask]
            group_2 <- gene_expr[!mask]

            # Perform Wilcoxon rank-sum test
            suppressWarnings(wilcox_result <- wilcox.test(group_1, group_2))
            log2FC <- avg_log2(group_1) - avg_log2(group_2)

            # Store log2FC and p-value
            DE_result_df$avg_log2FC[i] <- log2FC
            DE_result_df$p_value[i] <- wilcox_result$p.value
        }

        # Adjust the p-values for multiple testing (e.g., using Benjamini-Hochberg method
        DE_result_df$adj_p_value <- p.adjust(DE_result_df$p_value, method = "BH")

        # Return the result data frame
        DE_result_df <- DE_result_df[order(DE_result_df$avg_log2FC, decreasing = TRUE),]
        return(DE_result_df)

    } else if (!is.null(ident.1) & !is.null(ident.2)) { # Perform a comparison of ident.1 against ident.2        
        for (i in 1:nrow(expr_matrix)) {
            gene_expr <- expr_matrix[i,]

            # Split gene expression values based on group labels
            mask_1 <- annot_matrix[, group_column] == ident.1
            mask_2 <- annot_matrix[, group_column] == ident.2
            group_1 <- gene_expr[mask_1]
            group_2 <- gene_expr[mask_2]

            # Perform Wilcoxon rank-sum test
            suppressWarnings(wilcox_result <- wilcox.test(group_1, group_2))
            log2FC <- avg_log2(group_1) - avg_log2(group_2)

            # Store log2FC and p-value
            DE_result_df$avg_log2FC[i] <- log2FC
            DE_result_df$p_value[i] <- wilcox_result$p.value
        }

        # Adjust the p-values for multiple testing (e.g., using Benjamini-Hochberg method
        DE_result_df$adj_p_value <- p.adjust(DE_result_df$p_value, method = "BH")

        # Return the result data frame
        DE_result_df <- DE_result_df[order(DE_result_df$avg_log2FC, decreasing = TRUE),]
        return(DE_result_df)
    }
}