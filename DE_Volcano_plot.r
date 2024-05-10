# Define a function for plotting

DE_Volcano_plot <- function(DE_result_mat, log2FC_sig_thresh = 1.5, p_val_sig_thresh = 0.05, top_n_labels = 20, title = "", text_size = 20, box_padding = 0.5, max_overlaps = 30) {
    # Supply a DE result matrix with the column names "avg_log2FC", "adj_p_value", and "gene"
    # log2FC_sig_thresh specifies a log2FC value that is considered sig up- or down-regulated (default: 1.5)
    # p_val_sig_thresh specifies a adj p-val that is considered significant (default: 0.05)
    # top_n_labels specifies how many up- and down-regulated genes to label (default: 20)

    suppressPackageStartupMessages({
        library(ggplot2)
        library(ggrepel)
        library(dplyr)
    })

    required_columns <- c("avg_log2FC", "genes", "adj_p_value")

    if (all(!(required_columns %in% colnames(DE_result_mat)))) {
        missing_columns <- required_columns[!(required_columns %in% colnames(DE_result_mat))]
        missing_columns_text <- paste(missing_columns, collapse = ", ")
        message(missing_columns_text, "were not fouund in the input DE matrix")
        break
    } else {        
        if (min(DE_result_mat$adj_p_value) > p_val_sig_thresh) { # In case no gene has p-val < p_val_sig_thresh, set ylim to accomodate hline and annotate top 20 based on log2FC only
            top_20 <- subset(DE_result_mat, subset = avg_log2FC > log2FC_sig_thresh) %>% top_n(top_n_labels, avg_log2FC) %>% pull("gene") %>% as.character()
            bot_20 <- subset(DE_result_mat, subset = avg_log2FC < -log2FC_sig_thresh) %>% top_n(-top_n_labels, avg_log2FC) %>% pull("gene") %>% as.character()
            annot <- c(top_20, bot_20) 
            
            plot <- ggplot(data = DE_result_mat, aes(x = avg_log2FC, y = -log(adj_p_value, 10), label = gene)) +
                            geom_point(size = 1) +
                            geom_vline(xintercept = log2FC_sig_thresh, linetype = "dashed", color = "blue", linewidth = 1) +
                            geom_vline(xintercept = -log2FC_sig_thresh, linetype = "dashed", color = "blue", linewidth = 1) +
                            geom_hline(yintercept = -log(p_val_sig_thresh, 10), linetype = "dashed", color = "black", linewidth = 1) +
                            scale_x_continuous(breaks = seq(floor(min(DE_result_mat$avg_log2FC)), round(max(DE_result_mat$avg_log2FC)), 1)) + 
                            scale_y_continuous(expand = c(0, 0), limits = c(NA, -log(p_val_sig_thresh, 10) * 1.05)) + # This will add 5% padding to the top of the y-axis
                            labs(x = "log2FC", y = "-log10(adj_p_value)") +
                            ggtitle(title) +
                            theme_light() +
                            theme(text = element_text(size = text_size)) +
                            geom_label_repel(data = subset(DE_result_mat, gene %in% annot),
                                            box.padding = box_padding,
                                            max.overlaps = max_overlaps)
        } else {
            # Pull the top 20 up- and down-regulated genes for annotation
            top_20 <- subset(DE_result_mat, subset = adj_p_value < p_val_sig_thresh & avg_log2FC > log2FC_sig_thresh) %>% top_n(top_n_labels, avg_log2FC) %>% pull("gene") %>% as.character()
            bot_20 <- subset(DE_result_mat, subset = adj_p_value < p_val_sig_thresh & avg_log2FC < -log2FC_sig_thresh) %>% top_n(-top_n_labels, avg_log2FC) %>% pull("gene") %>% as.character()
            annot <- c(top_20, bot_20) 

            plot <- ggplot(data = DE_result_mat, aes(x = avg_log2FC, y = -log(adj_p_value, 10), label = gene)) +
                            geom_point(size = 1) +
                            geom_vline(xintercept = log2FC_sig_thresh, linetype = "dashed", color = "blue", linewidth = 1) +
                            geom_vline(xintercept = -log2FC_sig_thresh, linetype = "dashed", color = "blue", linewidth = 1) +
                            geom_hline(yintercept = -log(p_val_sig_thresh, 10), linetype = "dashed", color = "black", linewidth = 1) +
                            scale_x_continuous(breaks = seq(floor(min(DE_result_mat$avg_log2FC)), round(max(DE_result_mat$avg_log2FC)), 1)) + 
                            scale_y_continuous(expand = c(0, 0), limits = c(NA, -(DE_result_mat$adj_p_value %>% min() %>% log(10) * 1.05))) + # This will add 5% padding to the top of the y-axis
                            labs(x = "log2FC", y = "-log10(adj_p_value)") +
                            ggtitle(title) +
                            theme_light() +
                            theme(text = element_text(size = text_size)) +
                            geom_label_repel(data = subset(DE_result_mat, gene %in% annot),
                                            box.padding = box_padding,
                                            max.overlaps = max_overlaps)
        }

    }
    return(plot)
}