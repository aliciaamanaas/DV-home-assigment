library(shiny)
library(ggplot2)
library(ggfortify)
library(factoextra)
library(plotly)
library(dplyr)
library(ggpubr)
library(gtsummary)
library(DT)
library(Rtsne)
library(stringr)

# Adjust this path to where your rep1.csv … rep4.csv files are stored
folder_path <- "/Users/paulaartizduenas/Desktop/DV - Group Project/data"

rep_files <- list(
  rep1 = "rep1.csv",
  rep2 = "rep2.csv",
  rep3 = "rep3.csv",
  rep4 = "rep4.csv"
)

rep_list <- lapply(rep_files, function(f) {
  m <- read.csv(file.path(folder_path, f), row.names = 1, check.names = FALSE)
  as.matrix(m)
})

exprs_data <- t(do.call(cbind, rep_list))#rows = samples, cols = genes

sample_names <- rownames(exprs_data)

meta_data <- data.frame(
  row.names        = sample_names,
  wine             = str_match(sample_names, "^([AB])\\.")[, 2],
  time             = as.integer(str_match(sample_names, "\\.t(\\d+)_")[, 2]),
  replicate        = as.factor(as.integer(str_match(sample_names, "_rep(\\d+)$")[, 2])),
  stringsAsFactors = FALSE
)

exprs_data_clean <- exprs_data[, colSums(exprs_data) > 0]
exprs_data_norm  <- exprs_data_clean / rowSums(exprs_data_clean)

meta_data$Expression <- rowSums(exprs_data_clean)

#PCA
pca_not_scaled <- prcomp(exprs_data,       scale = FALSE, center = TRUE)
pca_scaled     <- prcomp(exprs_data_clean, scale = TRUE)
pca_norm       <- prcomp(exprs_data_norm,  scale = TRUE)

# Variance labels for normalised PCA
var_exp <- summary(pca_norm)$importance[2, ] * 100
pc1_lab <- paste0("PC1 (", round(var_exp[1], 2), "%)")
pc2_lab <- paste0("PC2 (", round(var_exp[2], 2), "%)")

# Batch effect
count_df    <- data.frame(row.names = rownames(exprs_data_clean),
                          total_exp = rowSums(exprs_data_clean))
meta_data_all <- merge(count_df, meta_data, by = "row.names")
meta_data_all$Row.names <- as.factor(meta_data_all$Row.names)


pca_points       <- as_tibble(pca_scaled$x) |> bind_cols(meta_data) |> as.data.frame()
pca_points$exprs <- rowSums(exprs_data_clean)
pc1_mod          <- lm(PC1 ~ exprs, pca_points)
sm               <- summary(pc1_mod)
r2               <- sm$r.squared
pval             <- sm$coefficients[, 4][2]
color4points     <- if (r2 > 0.6 & pval < 0.05) "red" else if (pval < 0.05) "orange" else "blue"

#tSNE
set.seed(42)
tsne_raw  <- Rtsne(exprs_data_clean, dims = 2, perplexity = 10, verbose = FALSE)
tsne_norm <- Rtsne(exprs_data_norm,  dims = 2, perplexity = 10, verbose = FALSE)

df_raw <- data.frame(
  tSNE1      = tsne_raw$Y[, 1],  tSNE2 = tsne_raw$Y[, 2],
  wine       = meta_data$wine,   replicate = meta_data$replicate,
  Expression = meta_data$Expression, time = as.factor(meta_data$time)
)
df_norm <- data.frame(
  tSNE1      = tsne_norm$Y[, 1], tSNE2 = tsne_norm$Y[, 2],
  wine       = meta_data$wine,   replicate = meta_data$replicate,
  Expression = meta_data$Expression, time = as.factor(meta_data$time)
)

set.seed(42)
tsne_res_1 <- Rtsne(exprs_data_norm, dims = 2, perplexity = 10)
tsne_res_2 <- Rtsne(exprs_data_norm, dims = 2, perplexity = 10)
tsne_res_3 <- Rtsne(exprs_data_norm, dims = 2, perplexity = 10)

set.seed(42)
tsne_perp5  <- Rtsne(exprs_data_norm, dims = 2, perplexity = 5)
tsne_perp10 <- Rtsne(exprs_data_norm, dims = 2, perplexity = 10)
tsne_perp15 <- Rtsne(exprs_data_norm, dims = 2, perplexity = 15)

set.seed(123)
tsne_iter100  <- Rtsne(exprs_data_norm, dims = 2, perplexity = 15, max_iter = 100)
tsne_iter500  <- Rtsne(exprs_data_norm, dims = 2, perplexity = 15, max_iter = 500)
tsne_iter1000 <- Rtsne(exprs_data_norm, dims = 2, perplexity = 15, max_iter = 1000)

pca_scores <- as.data.frame(pca_norm$x)
tukey_outlier <- function(x) {
  Q1 <- quantile(x, 0.25); Q3 <- quantile(x, 0.75); IQR_val <- Q3 - Q1
  x < (Q1 - 1.5 * IQR_val) | x > (Q3 + 1.5 * IQR_val)
}
is_outlier         <- tukey_outlier(pca_scores$PC1) | tukey_outlier(pca_scores$PC2)
meta_data$outlier  <- is_outlier
df_norm$outlier    <- is_outlier

exprs_norm_clean <- exprs_data_norm[!is_outlier, ]
meta_clean       <- meta_data[!is_outlier, ]
pca_no_out       <- prcomp(exprs_norm_clean, scale = TRUE)
var_no_out       <- summary(pca_no_out)$importance[2, ] * 100

set.seed(42)
tsne_final <- Rtsne(exprs_data_norm, dims = 2, perplexity = 15, max_iter = 1000)
df_final   <- data.frame(
  tSNE1      = tsne_final$Y[, 1], tSNE2 = tsne_final$Y[, 2],
  wine       = meta_data$wine,    replicate = meta_data$replicate,
  time       = meta_data$time,    Expression = meta_data$Expression,
  outlier    = meta_data$outlier
)

ui <- fluidPage(
  titlePanel(HTML("Home Assignment 1<br>
                   <font color='#F48225' size='4'>🍷 Wine production – RNAseq of the Flor yeast (Course 2025/2026)</font><br>
                   <font color='#A8A8A8' size='3'>Alicia Mañas, Lídia Sánchez and Paula Artiz</font>")),
  
  tags$head(tags$style(HTML("
    .answer-box  { background:#f0fff4; border-left:5px solid #28a745; padding:15px; margin-bottom:20px; border-radius:5px; }
    .code-box    { background:#f5f5f5; border:1px solid #ccc; padding:15px; margin-bottom:20px; border-radius:5px; font-family:monospace; overflow-x:auto; }
    .question-box{ background:#f4f8ff; border-left:5px solid #3366cc; padding:15px; margin-bottom:20px; border-radius:5px; }
    .alt-box     { background:#fff8f8; border-left:5px solid red; padding:15px; margin-bottom:20px; border-radius:5px; }
    .tab-content { padding:20px; }
    body { font-family:'Helvetica Neue', Arial, sans-serif; }
    h2 { border-bottom:2px solid #F48225; padding-bottom:4px; }
    .sidebar-btn { margin-bottom:8px; text-align:left; width:100%; }
  "))),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Navigation"),
      br(),
      actionButton("btn_q1", "Part 1 · Data Overview", class = "sidebar-btn"),
      actionButton("btn_q2", "Part 2 · PCA (Not Scaled)", class = "sidebar-btn"),
      actionButton("btn_q3", "Part 2 · PCA (Scaled)", class = "sidebar-btn"),
      actionButton("btn_q4", "Part 2 · Batch Effect", class = "sidebar-btn"),
      actionButton("btn_q5", "Part 2 · PC1 ~ Expression", class = "sidebar-btn"),
      actionButton("btn_q6", "Part 3 · Normalised PCA", class = "sidebar-btn"),
      actionButton("btn_q7", "Part 4 · tSNE Raw vs Norm", class = "sidebar-btn"),
      actionButton("btn_q8", "Part 4 · tSNE Parameters", class = "sidebar-btn"),
      actionButton("btn_q9", "Part 5 · Final Interpretation", class = "sidebar-btn")
    ),
    
    mainPanel(
      width = 9,
      
      div(
        h2("Introduction"),
        p("Flor yeasts (Saccharomyces cerevisiae) are the main microorganisms involved in the biological aging of sherry wines, where they form a velum on the surface of fortified must and enable survival under highly stressful winemaking conditions. This biological aging process contributes to key technological properties of sherry wines. In this study, we analyze two industrial flor yeast strains used for the production of sherry wines A and B, which differ in flavor, alcohol percentage, and color. The goal is to determine whether these differences are associated with genetic or transcriptional variation between strains, or with other biochemical processes occurring during wine production. RNAseq data were generated across six stages of velum development with four biological replicates per condition, resulting in 48 samples. After standard quality control, mapping, and gene quantification, the resulting count matrix is used to explore global expression patterns. The main objective is to apply dimensionality reduction techniques to assess whether the two yeast strains show distinct transcriptomic profiles that could explain differences between the wine products.")
      ),
      
      uiOutput("dynamic_content")
    )
  )
)

server <- function(input, output, session) {
  
  current_section <- reactiveVal("q1")
  
  observeEvent(input$btn_q1, current_section("q1"))
  observeEvent(input$btn_q2, current_section("q2"))
  observeEvent(input$btn_q3, current_section("q3"))
  observeEvent(input$btn_q4, current_section("q4"))
  observeEvent(input$btn_q5, current_section("q5"))
  observeEvent(input$btn_q6, current_section("q6"))
  observeEvent(input$btn_q7, current_section("q7"))
  observeEvent(input$btn_q8, current_section("q8"))
  observeEvent(input$btn_q9, current_section("q9"))
  
  output$q1_dt <- renderDT({
    df <- data.frame(
      Sample    = rownames(meta_data),
      Wine      = meta_data$wine,
      Time      = meta_data$time,
      Replicate = as.character(meta_data$replicate),
      Total_Exp = round(meta_data$Expression, 2)
    )
    datatable(df, options = list(pageLength = 15, scrollX = TRUE))
  })
  output$q1_summary <- renderPrint({
    cat("Samples per wine type:\n"); print(table(meta_data$wine))
    cat("\nSamples per replicate:\n"); print(table(meta_data$replicate))
    cat("\nSamples per time point:\n"); print(table(meta_data$time))
  })
  
  output$q2_fviz_wine  <- renderPlot({ fviz_pca_ind(pca_not_scaled, geom.ind="point", col.ind=meta_data$wine,      addEllipses=TRUE, legend.title="Wine") })
  output$q2_auto_wine  <- renderPlot({ autoplot(pca_not_scaled, data=meta_data, colour="wine")      + theme_classic() })
  output$q2_base_wine  <- renderPlot({
    plot(pca_not_scaled$x[,1], pca_not_scaled$x[,2], col=as.factor(meta_data$wine), pch=19, xlab="PC1", ylab="PC2", main="PCA Not Scaled — Wine")
    legend("topright", legend=levels(as.factor(meta_data$wine)), col=1:2, pch=19)
  })
  output$q2_fviz_rep   <- renderPlot({ fviz_pca_ind(pca_not_scaled, geom.ind="point", col.ind=meta_data$replicate, addEllipses=TRUE, legend.title="Replicate") })
  output$q2_auto_rep   <- renderPlot({ autoplot(pca_not_scaled, data=meta_data, colour="replicate") + theme_classic() })
  output$q2_base_rep   <- renderPlot({
    plot(pca_not_scaled$x[,1], pca_not_scaled$x[,2], col=as.integer(meta_data$replicate), pch=19, xlab="PC1", ylab="PC2", main="PCA Not Scaled — Replicate")
    legend("topright", legend=levels(meta_data$replicate), col=1:4, pch=19)
  })
  output$q2_scree_fviz <- renderPlot({ fviz_eig(pca_not_scaled, ncp=8, barfill="steelblue", barcolor="steelblue", main="Variance Explained") })
  output$q2_scree_base <- renderPlot({ plot(pca_not_scaled, main="Variance Explained (base)") })
  
  output$q3_fviz_wine  <- renderPlot({ fviz_pca_ind(pca_scaled, geom.ind="point", col.ind=meta_data$wine,      addEllipses=TRUE, legend.title="Wine") })
  output$q3_auto_wine  <- renderPlot({ autoplot(pca_scaled, data=meta_data, colour="wine")      + theme_classic() })
  output$q3_fviz_rep   <- renderPlot({ fviz_pca_ind(pca_scaled, geom.ind="point", col.ind=meta_data$replicate, addEllipses=TRUE, legend.title="Replicate") })
  output$q3_auto_rep   <- renderPlot({ autoplot(pca_scaled, data=meta_data, colour="replicate") + theme_classic() })
  output$q3_both_auto  <- renderPlot({ autoplot(pca_scaled, data=meta_data, colour="replicate", shape="wine") + theme_classic() + scale_shape_manual(values=c(3,17)) })
  output$q3_scree      <- renderPlot({ fviz_eig(pca_scaled, barfill="steelblue", barcolor="steelblue", main="Variance Explained — Scaled") })
  
  output$q4_bar_rep   <- renderPlot({
    ggplot(meta_data_all, aes(x=Row.names, y=total_exp, fill=replicate)) +
      geom_bar(stat="identity") + theme_classic() +
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      labs(title="Total Expression by Replicate", x="Sample", y="Total counts")
  })
  output$q4_bar_wine  <- renderPlot({
    ggplot(meta_data_all, aes(x=Row.names, y=total_exp, fill=wine)) +
      geom_bar(stat="identity") + theme_classic() +
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      labs(title="Total Expression by Wine type", x="Sample", y="Total counts")
  })
  output$q4_bar_grad  <- renderPlot({
    ggplot(meta_data_all, aes(x=Row.names, y=total_exp, fill=total_exp)) +
      geom_bar(stat="identity") + theme_classic() +
      scale_fill_gradient(low="green", high="red") +
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      labs(title="Total Expression Gradient", x="Sample", y="Total counts")
  })
  output$q4_box_rep   <- renderPlot({
    ggboxplot(meta_data_all, x="replicate", y="total_exp", fill="replicate", add="jitter", palette="jama") +
      stat_compare_means(method="wilcox.test") + labs(title="Expression by Replicate")
  })
  output$q4_box_wine  <- renderPlot({
    ggboxplot(meta_data_all, x="wine", y="total_exp", fill="wine", add="jitter", palette="jama") +
      stat_compare_means() + labs(title="Expression by Wine type")
  })
  output$q4_box_facet <- renderPlot({
    ggboxplot(meta_data_all, x="wine", y="total_exp", fill="wine", add="jitter", palette="jama") +
      stat_compare_means() + facet_grid(~replicate) + labs(title="Expression by Wine — faceted by Replicate")
  })
  output$q4_base_rep  <- renderPlot({
    plot(rowSums(exprs_data_clean), type="h", col=as.integer(meta_data$replicate), ylab="Total counts", xlab="Sample", main="Total Expression (Base R)")
    legend("topright", legend=levels(meta_data$replicate), col=1:4, pch=15)
  })
  
  output$q5_scatter <- renderPlot({
    ggscatter(pca_points, x="PC1", y="exprs", add="reg.line", conf.int=TRUE, cor.coef=TRUE,
              color=color4points, alpha=0.6, ggtheme=theme_bw(), xlab="PC1", ylab="Total expression") +
      labs(title="PC1 vs Total Expression")
  })
  output$q5_pca_expr <- renderPlot({
    autoplot(pca_scaled, data=meta_data, colour="Expression") +
      theme_classic() + scale_color_gradient(low="green", high="red") +
      labs(title="PCA Scaled — coloured by Total Expression")
  })
  output$q5_pca_shape <- renderPlot({
    autoplot(pca_scaled, data=meta_data, colour="Expression", shape="replicate") +
      theme_classic() + scale_color_gradient(low="green", high="red") +
      scale_shape_manual(values=c(3,17,15,18)) +
      labs(title="PCA Scaled — Expression + Replicate")
  })
  output$q5_lm <- renderPrint({ print(sm) })
  
  output$q6_rep   <- renderPlot({ autoplot(pca_norm, data=meta_data, colour="replicate") + theme_classic() + labs(title="Normalised PCA — Replicate") })
  output$q6_wine  <- renderPlot({ autoplot(pca_norm, data=meta_data, colour="wine")      + theme_classic() + labs(title="Normalised PCA — Wine type") })
  output$q6_both  <- renderPlot({ autoplot(pca_norm, data=meta_data, colour="replicate", shape="wine") + theme_classic() + scale_shape_manual(values=c(3,17)) + labs(title="Normalised PCA — Replicate + Wine") })
  output$q6_expr  <- renderPlot({ autoplot(pca_norm, data=meta_data, colour="Expression", shape="wine") + theme_classic() + scale_color_gradient(low="green", high="red") + scale_shape_manual(values=c(3,17)) + labs(title="Normalised PCA — Expression gradient") })
  output$q6_facet <- renderPlot({ autoplot(pca_norm, data=meta_data, colour="Expression", shape="wine") + theme_classic() + scale_color_gradient(low="green", high="red") + facet_grid(~replicate) + scale_shape_manual(values=c(3,17)) + labs(title="Normalised PCA — faceted by Replicate") })
  output$q6_base  <- renderPlot({
    par(mfrow=c(1,2))
    plot(pca_norm$x, pch=19, col=as.integer(meta_data$replicate), main="Replicate")
    legend("topright", legend=levels(meta_data$replicate), col=1:4, pch=19)
    plot(pca_norm$x, pch=19, col=as.factor(meta_data$wine),       main="Wine type")
    legend("topright", legend=levels(as.factor(meta_data$wine)), col=1:2, pch=19)
  })
  
  output$q7_raw_wine  <- renderPlot({ ggplot(df_raw, aes(tSNE1,tSNE2,colour=wine)) + geom_point(size=2.5,alpha=0.8) + scale_colour_manual(values=c("A"="#E74C3C","B"="#3498DB")) + theme_bw() + labs(title="tSNE Raw — Wine type",colour="Wine type") })
  output$q7_raw_rep   <- renderPlot({ ggplot(df_raw, aes(tSNE1,tSNE2,colour=replicate)) + geom_point(size=2.5,alpha=0.8) + scale_colour_brewer(palette="Set1") + theme_bw() + labs(title="tSNE Raw — Replicate") })
  output$q7_raw_expr  <- renderPlot({ ggplot(df_raw, aes(tSNE1,tSNE2,colour=Expression)) + geom_point(size=2.5,alpha=0.8) + scale_colour_viridis_c(option="plasma") + theme_bw() + labs(title="tSNE Raw — Total Expression") })
  output$q7_norm_wine <- renderPlot({ ggplot(df_norm, aes(tSNE1,tSNE2,colour=wine)) + geom_point(size=2.5,alpha=0.8) + scale_colour_manual(values=c("A"="#E74C3C","B"="#3498DB")) + theme_bw() + labs(title="tSNE Norm — Wine type",colour="Wine type") })
  output$q7_norm_rep  <- renderPlot({ ggplot(df_norm, aes(tSNE1,tSNE2,colour=replicate)) + geom_point(size=2.5,alpha=0.8) + scale_colour_brewer(palette="Set1") + theme_bw() + labs(title="tSNE Norm — Replicate") })
  output$q7_norm_expr <- renderPlot({ ggplot(df_norm, aes(tSNE1,tSNE2,colour=Expression)) + geom_point(size=2.5,alpha=0.8) + scale_colour_viridis_c(option="plasma") + theme_bw() + labs(title="tSNE Norm — Total Expression") })
  output$q7_norm_time <- renderPlot({ ggplot(df_norm, aes(tSNE1,tSNE2,colour=time)) + geom_point(size=2.5,alpha=0.8) + scale_colour_brewer(palette="RdYlBu",direction=-1) + theme_bw() + labs(title="tSNE Norm — Time point") })
  
  output$q8_seed <- renderPlot({
    make_df <- function(res,lbl) data.frame(x=res$Y[,1], y=res$Y[,2], run=lbl)
    df_runs <- rbind(make_df(tsne_res_1,"Run 1"), make_df(tsne_res_2,"Run 2"), make_df(tsne_res_3,"Run 3"))
    ggplot(df_runs, aes(x,y)) + geom_point(alpha=0.7) + facet_wrap(~run,nrow=1) + labs(title="tSNE — Effect of random seed",x="tSNE1",y="tSNE2") + theme_bw()
  })
  output$q8_perp <- renderPlot({
    mk <- function(res,lbl) data.frame(x=res$Y[,1], y=res$Y[,2], Perplexity=lbl, wine=meta_data$wine)
    df_perp <- rbind(mk(tsne_perp5,"Perp=5"), mk(tsne_perp10,"Perp=10"), mk(tsne_perp15,"Perp=15"))
    ggplot(df_perp, aes(x,y,colour=wine)) + geom_point(size=2,alpha=0.8) +
      scale_colour_manual(values=c("A"="#E74C3C","B"="#3498DB")) +
      facet_wrap(~Perplexity,nrow=1) + labs(title="tSNE — Effect of perplexity",x="tSNE1",y="tSNE2",colour="Wine") + theme_bw()
  })
  output$q8_iter <- renderPlot({
    mk <- function(res,lbl) data.frame(x=res$Y[,1], y=res$Y[,2], Iterations=lbl, wine=meta_data$wine)
    df_iter <- rbind(mk(tsne_iter100,"100 iter"), mk(tsne_iter500,"500 iter"), mk(tsne_iter1000,"1000 iter"))
    df_iter$Iterations <- factor(df_iter$Iterations, levels=c("100 iter","500 iter","1000 iter"))
    ggplot(df_iter, aes(x,y,colour=wine)) + geom_point(size=2,alpha=0.8) +
      scale_colour_manual(values=c("A"="#E74C3C","B"="#3498DB")) +
      facet_wrap(~Iterations,nrow=1) + labs(title="tSNE — Effect of iterations (perp=15)",x="tSNE1",y="tSNE2",colour="Wine") + theme_bw()
  })
  
  output$q9_pca_main  <- renderPlot({
    autoplot(pca_norm, data=meta_data, colour="wine", shape="replicate") +
      scale_colour_manual(values=c("A"="#E74C3C","B"="#3498DB"), labels=c("Wine A","Wine B")) +
      scale_shape_manual(values=c(16,17,15,3)) +
      labs(title="Final PCA — Normalised data", x=pc1_lab, y=pc2_lab, colour="Wine type", shape="Replicate") +
      theme_classic()
  })
  output$q9_pca_facet <- renderPlot({
    autoplot(pca_norm, data=meta_data, colour="wine", shape="replicate") +
      scale_colour_manual(values=c("A"="#E74C3C","B"="#3498DB")) +
      scale_shape_manual(values=c(16,17,15,3)) +
      facet_wrap(~replicate) +
      labs(title="PCA faceted by Replicate", x=pc1_lab, y=pc2_lab) + theme_bw()
  })
  output$q9_pca_noout <- renderPlot({
    pc1n <- paste0("PC1 (", round(var_no_out[1],2),"%)")
    pc2n <- paste0("PC2 (", round(var_no_out[2],2),"%)")
    autoplot(pca_no_out, data=meta_clean, colour="wine", shape="replicate") +
      scale_colour_manual(values=c("A"="#E74C3C","B"="#3498DB")) +
      scale_shape_manual(values=c(16,17,15,3)) +
      labs(title="PCA without Outlier (B.t1_rep3 removed)", x=pc1n, y=pc2n) + theme_classic()
  })
  output$q9_outlier_pca  <- renderPlot({
    autoplot(pca_norm, data=meta_data, colour="wine", shape="outlier") +
      scale_colour_manual(values=c("A"="#E74C3C","B"="#3498DB")) +
      scale_shape_manual(values=c(`FALSE`=16,`TRUE`=8), labels=c("Normal","Outlier")) +
      labs(title="PCA — Outlier highlighted (Tukey)", x=pc1_lab, y=pc2_lab) + theme_classic()
  })
  output$q9_outlier_tsne <- renderPlot({
    ggplot(df_final, aes(tSNE1,tSNE2,colour=wine,shape=outlier)) +
      geom_point(size=2.5,alpha=0.9) +
      scale_colour_manual(values=c("A"="#E74C3C","B"="#3498DB")) +
      scale_shape_manual(values=c(`FALSE`=16,`TRUE`=8), labels=c("Normal","Outlier")) +
      labs(title="tSNE — Outlier highlighted (Tukey)") + theme_bw()
  })
  output$q9_tsne_wine  <- renderPlot({ ggplot(df_final, aes(tSNE1,tSNE2,colour=wine)) + geom_point(size=2.5,alpha=0.8) + scale_colour_manual(values=c("A"="#E74C3C","B"="#3498DB"), labels=c("Wine A","Wine B")) + theme_bw() + labs(title="Final tSNE — Wine type",colour="Wine type") })
  output$q9_tsne_rep   <- renderPlot({ ggplot(df_final, aes(tSNE1,tSNE2,colour=replicate)) + geom_point(size=2.5,alpha=0.8) + scale_colour_brewer(palette="Set1") + theme_bw() + labs(title="Final tSNE — Replicate") })
  output$q9_tsne_time  <- renderPlot({ ggplot(df_final, aes(tSNE1,tSNE2,colour=as.factor(time))) + geom_point(size=2.5,alpha=0.8) + scale_colour_brewer(palette="RdYlBu",direction=-1) + theme_bw() + labs(title="Final tSNE — Time point",colour="Time point") })
  output$q9_outlier_id <- renderPrint({ cat("Outlier(s) detected:\n"); print(rownames(pca_scores)[is_outlier]) })
  
  output$dynamic_content <- renderUI({
    q <- current_section()
    
    if (q == "q1") {
      tagList(
        h2("Part 1 — Load Data & Metadata"),
        div(class="question-box", h4("Q1. How many samples belong to each wine type and to each replicate?")),
        div(class="answer-box",
            h4("ANSWER Q1:"),
            p("The dataset contains 48 samples in total: 24 from wine A and 24 from wine B."),
            p("Each wine was measured at 6 velum development time points (t1–t6) with 4 biological replicates per time point."),
            p("The wine variable plays the same role as treatment in the teacher's dataset, and replicate plays the role of experiment.")
        ),
        div(class="code-box", h4("CODE"), pre(
          'rep_list <- lapply(rep_files, function(f) as.matrix(read.csv(f, row.names=1)))
exprs_data <- t(do.call(cbind, rep_list))  # rows=samples, cols=genes

meta_data <- data.frame(
  wine      = str_match(sample_names, "^([AB])\\.")[, 2],
  time      = as.integer(str_match(sample_names, "\\.t(\\d+)_")[, 2]),
  replicate = as.factor(as.integer(str_match(sample_names, "_rep(\\d+)$")[, 2]))
)

table(rownames(meta_data) == rownames(exprs_data))  # verify order'
        )),
      h4("Sample Overview Table"),
      DTOutput("q1_dt"),
      br(),
      h4("Distribution Summary"),
      verbatimTextOutput("q1_summary")
      )
    }
    
    else if (q == "q2") {
      tagList(
        h2("Part 2 — PCA (Not Scaled)"),
        div(class="question-box", h4("Q2. How are samples distributed? Is there a batch effect associated with replicates?")),
        div(class="answer-box",
            h4("ANSWER Q2:"),
            p("Most of the variance does not correlate with wine type. However, there is a clear clustering by replicate: PC2 clearly separates samples from different replicates."),
            p("There is a noticeable batch effect caused by the replicates. Differences caused by technical aspects between replicates hide the important biological difference (wine A vs wine B) and need to be corrected through normalisation.")
        ),
        div(class="code-box", h4("CODE"), pre(
          'pca_res <- stats::prcomp(exprs_data, scale = FALSE)

# ggfortify
autoplot(pca_res, data = meta_data, colour = "wine") + theme_classic()
autoplot(pca_res, data = meta_data, colour = "replicate") + theme_classic()

# factoextra
fviz_pca_ind(pca_res, geom.ind = "point", col.ind = meta_data$wine, addEllipses = TRUE)
fviz_eig(pca_res, main = "Principal components: Variance explained")'
        )),
        h4("PCA Not Scaled — Coloured by Wine type"),
        tabsetPanel(
          tabPanel("factoextra", plotOutput("q2_fviz_wine", height="450px")),
          tabPanel("ggfortify",  plotOutput("q2_auto_wine", height="450px")),
          tabPanel("Base R",     plotOutput("q2_base_wine", height="450px"))
        ),
        h4("PCA Not Scaled — Coloured by Replicate"),
        tabsetPanel(
          tabPanel("factoextra", plotOutput("q2_fviz_rep",  height="450px")),
          tabPanel("ggfortify",  plotOutput("q2_auto_rep",  height="450px")),
          tabPanel("Base R",     plotOutput("q2_base_rep",  height="450px"))
        ),
        h4("Variance Explained"),
        tabsetPanel(
          tabPanel("factoextra", plotOutput("q2_scree_fviz", height="400px")),
          tabPanel("Base R",     plotOutput("q2_scree_base", height="400px"))
        )
      )
    }
    
    else if (q == "q3") {
      tagList(
        h2("Part 2 — PCA (Scaled)"),
        div(class="question-box", h4("Q3. Has the batch effect been corrected after scaling?")),
        div(class="answer-box",
            h4("ANSWER Q3:"),
            p("Scaling by itself does not eliminate the batch effect. Replicates continue to group apart, and now PC1 somewhat reflects the batch effect."),
            p("Scaling balances the gene inputs but cannot fix variations in total library size across samples — direct normalisation is required.")
        ),
        div(class="code-box", h4("CODE"), pre(
          'exprs_data_clean <- exprs_data[, colSums(exprs_data) > 0]
pca_res_scaled <- stats::prcomp(exprs_data_clean, scale = TRUE)

autoplot(pca_res_scaled, data = meta_data, colour = "wine")      + theme_classic()
autoplot(pca_res_scaled, data = meta_data, colour = "replicate") + theme_classic()
autoplot(pca_res_scaled, data = meta_data, colour = "replicate", shape = "wine") + theme_classic()

fviz_pca_ind(pca_res_scaled, col.ind = meta_data$wine,      addEllipses = TRUE)
fviz_pca_ind(pca_res_scaled, col.ind = meta_data$replicate, addEllipses = TRUE)'
        )),
        h4("Scaled PCA — Wine type"),
        tabsetPanel(
          tabPanel("factoextra", plotOutput("q3_fviz_wine", height="450px")),
          tabPanel("ggfortify",  plotOutput("q3_auto_wine", height="450px"))
        ),
        h4("Scaled PCA — Replicate"),
        tabsetPanel(
          tabPanel("factoextra", plotOutput("q3_fviz_rep",  height="450px")),
          tabPanel("ggfortify",  plotOutput("q3_auto_rep",  height="450px"))
        ),
        h4("Scaled PCA — Replicate (colour) + Wine (shape)"),
        plotOutput("q3_both_auto", height="450px"),
        h4("Variance Explained — Scaled"),
        plotOutput("q3_scree", height="400px")
      )
    }
    
    else if (q == "q4") {
      tagList(
        h2("Part 2 — Understanding the Batch Effect"),
        div(class="question-box", h4("Q4. Is there any difference between replicates in total expression?")),
        div(class="answer-box",
            h4("ANSWER Q4:"),
            p("Replicates show significant differences in overall expression, which confirms that variations in library size are causing the batch effect."),
            p("There is no significant difference in overall expression when comparing wine A and wine B within the same replicate."),
            p("Samples with high total expression are found on the left side of PC1 (negative values), while those with low expression are on the right side.")
        ),
        div(class="code-box", h4("CODE"), pre(
          'count_df <- data.frame(total_exp = rowSums(exprs_data_clean))
meta_data_all <- merge(count_df, meta_data, by = "row.names")

# Barplots
ggplot(meta_data_all, aes(x=Row.names, y=total_exp, fill=replicate)) + geom_bar(stat="identity")
ggplot(meta_data_all, aes(x=Row.names, y=total_exp, fill=wine))      + geom_bar(stat="identity")

# Boxplots with Wilcoxon
ggboxplot(meta_data_all, x="replicate", y="total_exp", fill="replicate", add="jitter") +
  stat_compare_means(method = "wilcox.test")'
        )),
        h4("Barplot by Replicate"),
        plotOutput("q4_bar_rep", height="380px"),
        h4("Barplot by Wine type"),
        plotOutput("q4_bar_wine", height="380px"),
        h4("Expression Gradient Barplot"),
        plotOutput("q4_bar_grad", height="380px"),
        h4("Boxplot by Replicate"),
        plotOutput("q4_box_rep", height="420px"),
        h4("Boxplot by Wine type"),
        plotOutput("q4_box_wine", height="380px"),
        h4("Boxplot faceted by Replicate"),
        plotOutput("q4_box_facet", height="400px"),
        div(class="alt-box", h4("CODE ALTERNATIVE — Base R"), plotOutput("q4_base_rep", height="350px"))
      )
    }
    
    else if (q == "q5") {
      tagList(
        h2("Part 2 — PC1 Correlation with Total Expression"),
        div(class="question-box", h4("Q4 (cont.). Do samples with high expression have positive or negative PC1? Is PC1 correlated with total expression?")),
        div(class="answer-box",
            h4("ANSWER:"),
            p("PC1 is strongly negatively correlated with total expression (R ≈ −0.99 in analogous teacher's data)."),
            p("PC1 is driven mainly by library size rather than genuine biological differences. Normalisation is therefore essential to reveal the real signal.")
        ),
        div(class="code-box", h4("CODE"), pre(
          'pca_points <- as_tibble(pca_scaled$x) |> bind_cols(meta_data) |> as.data.frame()
pca_points$exprs <- rowSums(exprs_data_clean)

pc1_mod <- lm(PC1 ~ exprs, pca_points)
summary(pc1_mod)

ggscatter(pca_points, x = "PC1", y = "exprs", add = "reg.line",
          conf.int = TRUE, cor.coef = TRUE, ggtheme = theme_bw())'
        )),
        h4("PC1 vs Total Expression — Scatter & Linear Fit"),
        plotOutput("q5_scatter", height="500px"),
        h4("PCA Scaled — coloured by Expression"),
        plotOutput("q5_pca_expr", height="450px"),
        h4("PCA Scaled — Expression + Replicate shape"),
        plotOutput("q5_pca_shape", height="450px"),
        h4("Linear Model Summary"),
        verbatimTextOutput("q5_lm")
      )
    }
    
    else if (q == "q6") {
      tagList(
        h2("Part 3 — Normalised PCA"),
        div(class="question-box", h4("Q5. After normalisation, what groups can be distinguished? What do PC1 and PC2 represent?")),
        div(class="answer-box",
            h4("ANSWER Q5:"),
            p("Once the data is properly normalized, the batch effect completely vanishes, allowing samples from all four replicates to blend in PCA space."),
            p("PC1 (69.14%) distinguishes between the two wine types (A and B), representing the primary transcriptional variation between yeast strains."),
            p("PC2 (24.71%) captures variation associated with velum development stages (t1–t6), reflecting temporal changes in gene expression during the aging process."),
            p("Together, PC1 and PC2 account for ~93.85% of total variance.")
        ),
        div(class="code-box", h4("CODE"), pre(
          'exprs_data_norm <- exprs_data_clean / rowSums(exprs_data_clean)
pca_norm <- stats::prcomp(exprs_data_norm, scale = TRUE)

autoplot(pca_norm, data = meta_data, colour = "replicate") + theme_classic()
autoplot(pca_norm, data = meta_data, colour = "wine")      + theme_classic()
autoplot(pca_norm, data = meta_data, colour = "Expression", shape = "wine") +
  scale_color_gradient(low = "green", high = "red") + theme_classic()'
        )),
        tabsetPanel(
          tabPanel("Replicate",          plotOutput("q6_rep",   height="480px")),
          tabPanel("Wine type",          plotOutput("q6_wine",  height="480px")),
          tabPanel("Replicate + Wine",   plotOutput("q6_both",  height="480px")),
          tabPanel("Expression gradient",plotOutput("q6_expr",  height="480px")),
          tabPanel("Faceted",            plotOutput("q6_facet", height="480px")),
          tabPanel("Base R",             plotOutput("q6_base",  height="480px"))
        )
      )
    }
    
    else if (q == "q7") {
      tagList(
        h2("Part 4 — tSNE: Raw vs Normalised Data"),
        div(class="question-box", h4("Q3a/Q3b. How are samples distributed in tSNE? Is the batch effect visible?")),
        div(class="answer-box",
            h4("ANSWER Q3a (Raw):"),
            p("In the raw tSNE plot, samples form distinct clusters primarily according to replicate rather than biological condition. This is clear evidence of a batch effect dominating the data structure.")
        ),
        div(class="answer-box",
            h4("ANSWER Q3b (Normalised):"),
            p("Once normalised, samples mainly group by wine type. The batch-related effects no longer dominate the clustering. The clear separation between Wine A and Wine B shows that normalisation maintained important biological differences while minimising technical variation.")
        ),
        div(class="code-box", h4("CODE"), pre(
          'library(Rtsne)
exprs_clean <- exprs_data[, colSums(exprs_data) > 0]
exprs_norm  <- exprs_clean / rowSums(exprs_clean)

set.seed(42)
tsne_raw  <- Rtsne(exprs_clean, dims = 2, perplexity = 10)
tsne_norm <- Rtsne(exprs_norm,  dims = 2, perplexity = 10)

ggplot(df_raw,  aes(tSNE1, tSNE2, colour = wine))      + geom_point() + theme_bw()
ggplot(df_norm, aes(tSNE1, tSNE2, colour = wine))      + geom_point() + theme_bw()
ggplot(df_norm, aes(tSNE1, tSNE2, colour = replicate)) + geom_point() + theme_bw()'
        )),
        h3("Raw Data"),
        tabsetPanel(
          tabPanel("Wine type",       plotOutput("q7_raw_wine",  height="450px")),
          tabPanel("Replicate",       plotOutput("q7_raw_rep",   height="450px")),
          tabPanel("Total Expression",plotOutput("q7_raw_expr",  height="450px"))
        ),
        h3("Normalised Data"),
        tabsetPanel(
          tabPanel("Wine type",       plotOutput("q7_norm_wine", height="450px")),
          tabPanel("Replicate",       plotOutput("q7_norm_rep",  height="450px")),
          tabPanel("Total Expression",plotOutput("q7_norm_expr", height="450px")),
          tabPanel("Time point",      plotOutput("q7_norm_time", height="450px"))
        )
      )
    }
    
    else if (q == "q8") {
      tagList(
        h2("Part 4 — tSNE Parameters"),
        div(class="question-box", h4("Q4a/Q4b/Q4c. Effect of seed, perplexity and number of iterations.")),
        div(class="answer-box",
            h4("ANSWER Q4a — Seed:"),
            p("Without a fixed seed, clusters remain broadly consistent but their positions, orientations and relative distances change each time. tSNE is a stochastic method that initialises coordinates randomly, converging to different (though equally valid) local minima.")
        ),
        div(class="answer-box",
            h4("ANSWER Q4b — Perplexity:"),
            p("Low perplexity (5): many small, fragmented clusters. Intermediate (10–15): better balance between local and global structure. The most interpretable result uses perplexity = 15. Too high a value risks merging distinct groups.")
        ),
        div(class="answer-box",
            h4("ANSWER Q4c — Iterations:"),
            p("100 iterations: poor convergence; messy layout. 500 iterations: well-defined clusters already visible. 1000 iterations: slightly more refined, marginal improvement over 500. A reasonable minimum is 500 iterations.")
        ),
        div(class="code-box", h4("CODE"), pre(
          '# Random seed effect
tsne_1 <- Rtsne(exprs_norm, dims=2, perplexity=10)
tsne_2 <- Rtsne(exprs_norm, dims=2, perplexity=10)
tsne_3 <- Rtsne(exprs_norm, dims=2, perplexity=10)

# Perplexity effect
set.seed(42)
tsne_p5  <- Rtsne(exprs_norm, dims=2, perplexity=5)
tsne_p10 <- Rtsne(exprs_norm, dims=2, perplexity=10)
tsne_p15 <- Rtsne(exprs_norm, dims=2, perplexity=15)

# Iterations effect
set.seed(123)
tsne_i100  <- Rtsne(exprs_norm, dims=2, perplexity=15, max_iter=100)
tsne_i500  <- Rtsne(exprs_norm, dims=2, perplexity=15, max_iter=500)
tsne_i1000 <- Rtsne(exprs_norm, dims=2, perplexity=15, max_iter=1000)'
        )),
        h4("Effect of Random Seed (no set.seed)"),
        plotOutput("q8_seed", height="400px"),
        h4("Effect of Perplexity"),
        plotOutput("q8_perp", height="400px"),
        h4("Effect of Number of Iterations (perplexity = 15)"),
        plotOutput("q8_iter", height="400px")
      )
    }
    
    else if (q == "q9") {
      tagList(
        h2("Part 5 — Final Interpretation"),
        div(class="answer-box",
            h4("FINAL SUMMARY:"),
            p("Both PCA and t-SNE consistently show that the primary driver of variation in this dataset is wine type: samples from wine A and wine B form two clearly separated groups."),
            p("PC1 (69.14%) captures the inter-wine difference, while PC2 (24.71%) reflects temporal changes in gene expression across velum development stages (t1–t6)."),
            p("Replicates are well distributed within each group, indicating good technical reproducibility and no strong batch effect after normalisation."),
            p("One outlier was identified — B.t1_rep3 — which disproportionately distorts the PCA but can be removed on the basis of the Tukey criterion without altering the overall biological interpretation."),
            p("The differences between wine A and wine B are driven by the yeast strains used, and this signal is robust across both dimensionality reduction methods.")
        ),
        div(class="question-box", h4("Q5a/Q5b/Q5c. What do PC1 and PC2 represent? Compare PCA vs tSNE. Was an outlier detected?")),
        div(class="code-box", h4("CODE — Outlier detection"), pre(
          'pca_scores <- as.data.frame(pca_norm$x)
tukey_outlier <- function(x) {
  Q1 <- quantile(x, 0.25); Q3 <- quantile(x, 0.75); IQR_val <- Q3 - Q1
  x < (Q1 - 1.5*IQR_val) | x > (Q3 + 1.5*IQR_val)
}
is_outlier <- tukey_outlier(pca_scores$PC1) | tukey_outlier(pca_scores$PC2)
rownames(pca_scores)[is_outlier]  # → "B.t1_rep3"'
        )),
        h4("Outlier identified"),
        verbatimTextOutput("q9_outlier_id"),
        h4("Final PCA (normalised)"),
        tabsetPanel(
          tabPanel("Wine + Replicate",     plotOutput("q9_pca_main",  height="500px")),
          tabPanel("Faceted by Replicate", plotOutput("q9_pca_facet", height="500px")),
          tabPanel("Outlier highlighted",  plotOutput("q9_outlier_pca", height="500px")),
          tabPanel("Without Outlier",      plotOutput("q9_pca_noout", height="500px"))
        ),
        h4("Final tSNE (normalised, perplexity = 15, 1000 iter)"),
        tabsetPanel(
          tabPanel("Wine type",          plotOutput("q9_tsne_wine", height="480px")),
          tabPanel("Replicate",          plotOutput("q9_tsne_rep",  height="480px")),
          tabPanel("Time point",         plotOutput("q9_tsne_time", height="480px")),
          tabPanel("Outlier highlighted",plotOutput("q9_outlier_tsne", height="480px"))
        )
      )
    }
  })
}

shinyApp(ui = ui, server = server)