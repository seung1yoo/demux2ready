

### Load_libraries

library(dplyr)
library(ggplot2)
library(PCAtools)
library(EnhancedVolcano)
library(vidger)
library(RColorBrewer)
library(psych)
library(corrplot)


### Pair scatter plots
setwd("/workspace/project/GLC23009-RNAS_230817/rstudio")
file_name <- "../MultiParser.Cufflinks.genes.exp.fpkm.tsv"
data <- read.table(file_name, header = TRUE, sep = "\t", row.names = 1)
sample_count <- length(colnames(data))-5
data_fpkm <- data[,c(1:sample_count)]
colnames(data_fpkm)
data_fpkm_sub <- data_fpkm[which(apply(data_fpkm, 1, max) > 0.3),]
data_fpkm_sub <- data_fpkm_sub[which(apply(data_fpkm_sub, 1, max) < 3000),]
pairs(data_fpkm_sub, 
      row1attop = TRUE, gap = 0.3, font.labels = 2, 
      lower.panel = panel.cor, upper.panel = upper.panel)


corPlot(data_fpkm_sub, cex = 0.6, zlim = c(0.8, 1),
        main = "", xlas=0, ylas=2, ysrt=0, xsrt=45,
        show.legend=FALSE)

corrplot.mixed(cor(data_fpkm_sub, method = "pearson"), 
               tl.pos = 'lt', diag = 'l', 
               order="hclust", lower="color",upper="number")








### PCA
setwd("/workspace/project/GLC23009-RNAS_230817/rstudio")
file_name <- "../MultiParser.Cufflinks.genes.exp.fpkm.tsv"
data <- read.table(file_name, header = TRUE, sep = "\t", row.names = 1)
sample_count <- length(colnames(data))-5
data_fpkm <- data[,c(1:sample_count)]
colnames(data_fpkm)
group_names <- colnames(data_fpkm)
group_names <- gsub("_Lot1", "", group_names)
group_names <- gsub("_Lot2", "", group_names)
group_names <- gsub("_Lot3", "", group_names)
group_names
data_meta <- data.frame(row.names = colnames(data_fpkm),
                        sample = colnames(data_fpkm),
                        group = group_names)

p <- pca(data_fpkm, metadata = data_meta, removeVar = 0.1)
p$n <- length(p$components)
pcaPlotting_screeplot(p)
pcaPlotting_biplot(p)
pcaPlotting_pairsplot(p)
pcaPlotting_loadingsplot(p)

### Volcano plot
data <- read.table("../MultiParser.Cuffdiff.genes.deg.tsv", sep = "\t", header = TRUE, row.names = 1)
colnames(data)
data_anno <- select(data, !starts_with("DEG"))

comp_list <- list(
  list(comp_id = "DEG001", control = "YiP3_P15", case = "EBF_P3"),
  list(comp_id = "DEG002", control = "YiP3_P15", case = "EBF_P4"),
  list(comp_id = "DEG003", control = "YiP3_P15", case = "EBF_P5"),
  list(comp_id = "DEG004", control = "YiP3_P15", case = "EBF_P7"),
  list(comp_id = "DEG005", control = "YiP3_P15", case = "EBF_P10"),
  list(comp_id = "DEG006", control = "YiP3_P15", case = "hMSC_P4"),
  list(comp_id = "DEG007", control = "YiP3_P15", case = "hMSC_P7"),
  list(comp_id = "DEG008", control = "EBF_P3", case = "EBF_P4"),
  list(comp_id = "DEG009", control = "EBF_P3", case = "EBF_P5"),
  list(comp_id = "DEG010", control = "EBF_P3", case = "EBF_P7"),
  list(comp_id = "DEG011", control = "EBF_P3", case = "EBF_P10"),
  list(comp_id = "DEG012", control = "EBF_P4", case = "hMSC_P4"),
  list(comp_id = "DEG013", control = "EBF_P7", case = "hMSC_P7"),
  list(comp_id = "DEG014", control = "hMSC_P4", case = "hMSC_P7")
)

for (comp_info in comp_list){
  print(comp_info$comp_id)
  data_deg <- select(data, starts_with(comp_info$comp_id))
  figure <- EnhancedVolcano(data_deg,
                  lab = data_anno$GeneName,
                  x = paste(comp_info$comp_id, "Log2FC", sep = "."),
                  y = paste(comp_info$comp_id, "p", sep = "."),
                  title = "Volcano plot", subtitle = paste(comp_info$control, comp_info$case, sep = " vs. "),
                  border = "full", col=c('grey', 'darkgrey', 'green', 'skyblue'),
                  legendPosition = "right",
                  pCutoff = 0.05, FCcutoff = 0.584,
                  ylim = c(0, 5.5), drawConnectors = TRUE)
  ggsave(paste("VolcanoPlot", comp_info$comp_id, "png", sep = "."),
         plot=figure, width = 11, height = 9, dpi = 350,  units = "in")
}
figure

# MA plot
View(df.cuff)
summary(df.cuff)

comp_list <- list(
  list(comp_id = "DEG001", control = "YiP3_P15", case = "EBF_P3", df_cuffdiff = read.table("../cuffdiff/DEG001/gene_exp.diff", header = TRUE)),
  list(comp_id = "DEG002", control = "YiP3_P15", case = "EBF_P4", df_cuffdiff = read.table("../cuffdiff/DEG002/gene_exp.diff", header = TRUE)),
  list(comp_id = "DEG003", control = "YiP3_P15", case = "EBF_P5", df_cuffdiff = read.table("../cuffdiff/DEG003/gene_exp.diff", header = TRUE)),
  list(comp_id = "DEG004", control = "YiP3_P15", case = "EBF_P7", df_cuffdiff = read.table("../cuffdiff/DEG004/gene_exp.diff", header = TRUE)),
  list(comp_id = "DEG005", control = "YiP3_P15", case = "EBF_P10", df_cuffdiff = read.table("../cuffdiff/DEG005/gene_exp.diff", header = TRUE)),
  list(comp_id = "DEG006", control = "YiP3_P15", case = "hMSC_P4", df_cuffdiff = read.table("../cuffdiff/DEG006/gene_exp.diff", header = TRUE)),
  list(comp_id = "DEG007", control = "YiP3_P15", case = "hMSC_P7", df_cuffdiff = read.table("../cuffdiff/DEG007/gene_exp.diff", header = TRUE)),
  list(comp_id = "DEG008", control = "EBF_P3", case = "EBF_P4", df_cuffdiff = read.table("../cuffdiff/DEG008/gene_exp.diff", header = TRUE)),
  list(comp_id = "DEG009", control = "EBF_P3", case = "EBF_P5", df_cuffdiff = read.table("../cuffdiff/DEG009/gene_exp.diff", header = TRUE)),
  list(comp_id = "DEG010", control = "EBF_P3", case = "EBF_P7", df_cuffdiff = read.table("../cuffdiff/DEG010/gene_exp.diff", header = TRUE)),
  list(comp_id = "DEG011", control = "EBF_P3", case = "EBF_P10", df_cuffdiff = read.table("../cuffdiff/DEG011/gene_exp.diff", header = TRUE)),
  list(comp_id = "DEG012", control = "EBF_P4", case = "hMSC_P4", df_cuffdiff = read.table("../cuffdiff/DEG012/gene_exp.diff", header = TRUE)),
  list(comp_id = "DEG013", control = "EBF_P7", case = "hMSC_P7", df_cuffdiff = read.table("../cuffdiff/DEG013/gene_exp.diff", header = TRUE)),
  list(comp_id = "DEG014", control = "hMSC_P4", case = "hMSC_P7", df_cuffdiff = read.table("../cuffdiff/DEG014/gene_exp.diff", header = TRUE))
)

idx = 0
for (comp_info in comp_list){
  idx <- idx + 1
  print(idx)
  print(comp_info$comp_id)
  
  comp_info$comp_id
  hl_df <- filter(comp_info$df_cuffdiff, abs(log2.fold_change.)>=1.0, abs(log2.fold_change.)<=5.0, 
                  p_value < 0.05, q_value < 0.1, 
                  value_1 > 10, value_2 > 10, 
                  gene != c("-"))
  
  hl_id <- hl_df[order(-abs(hl_df$log2.fold_change.), hl_df$q_value),]$gene
  comp_list[[idx]]$hl_id <- hl_id
  
  comp_info$df_cuffdiff$log2.fold_change. <- log(comp_info$df_cuffdiff$value_2+1, 2) - log(comp_info$df_cuffdiff$value_1+1, 2)
  comp_list[[idx]]$df_cuffdiff$log2.fold_change. <- comp_info$df_cuffdiff$log2.fold_change.
  
  if (idx == 1) {
    df_rbind_cuffdiff <- comp_info$df_cuffdiff
  } else {
    df_rbind_cuffdiff <- rbind(df_rbind_cuffdiff, comp_info$df_cuffdiff)
  }
}

comp_list[[14]]$hl_id
comp_list[[1]]$df_cuffdiff$q_value
comp_info$df_cuffdiff$log2.fold_change.
df_rbind_cuffdiff$log2.fold_change.
dim(df_rbind_cuffdiff)
View(df_rbind_cuffdiff)
summary(df_rbind_cuffdiff)

### Figure. Boxplot 
figure <- vsBoxPlot(
  data = df_rbind_cuffdiff, d.factor = NULL, type = 'cuffdiff', title = TRUE, 
  legend = TRUE, grid = TRUE, fill.color = "Paired", aes = "viosumm", data.return = TRUE
)
figure
ggsave("Boxplot.png", plot=figure$plot, width = 9, height = 6, dpi = 350, units = "in")

### Figure. DEG Matrix 
figure <- vsDEGMatrix(data = df_rbind_cuffdiff, padj = 0.05, d.factor = NULL, type = 'cuffdiff', 
                      title = TRUE, legend = TRUE, grid = TRUE, data.return = TRUE
)
ggsave(
  "DEG.Matrix.png", 
  plot = figure$plot, width = 6, height = 6, dpi = 350, units = "in"
)

### Figure. ScatterPlots
for (comp_info in comp_list) {
  print(comp_info$comp_id)
  figure <- vsScatterPlot(x = comp_info$control, y = comp_info$case, 
                          data = df_rbind_cuffdiff, type = "cuffdiff",
                          d.factor = NULL, title = TRUE, grid = TRUE, data.return = TRUE)
  ggsave(paste("ScatterPlot", comp_info$comp_id, "png", sep = "."),
         plot = figure$plot, width = 6, height = 6, dpi = 350, units = "in")
}

#figure <- vsScatterMatrix(data = df_rbind_cuffdiff, d.factor = NULL, type = 'cuffdiff', 
#                          comp = NULL, title = TRUE, grid = TRUE, man.title = NULL, data.return = TRUE)
#ggsave("ScatterPlot.Matrix.png", 
#       plot = figure$plot, width = 12, height = 12, dpi = 350, units = "in")


### MA plot
for (comp_info in comp_list) {
  print(comp_info$comp_id)
  figure <- vsMAPlot(x = comp_info$control, y = comp_info$case, 
                     data = comp_info$df_cuffdiff, d.factor = NULL, type = 'cuffdiff',
                     padj = 0.05, lfc = 0.584, y.lim = NULL, 
                     title = TRUE, legend = TRUE, grid = TRUE, data.return = TRUE,
                     xaxis.title.size = 13, yaxis.title.size = 13, main.title.size = 13,
                     highlight = comp_info$hl_id)
  ggsave(paste("MAPlot", comp_info$comp_id, "png", sep = "."),
         plot = figure$plot, width = 9, height = 6, dpi = 350, units = "in")
}

#vsMAMatrix(
#  data = df_rbind_cuffdiff, d.factor = NULL, type = 'cuffdiff', 
#  padj = 0.05, y.lim = NULL, lfc = 0.583, title = TRUE, 
#  grid = TRUE, counts = TRUE, data.return = FALSE
#)

### Volcano plot
for (comp_info in comp_list) {
  print(comp_info$comp_id)
  figure <- vsVolcano(
    x = comp_info$case, y = comp_info$control, 
    data = comp_info$df_cuffdiff, d.factor = NULL, type = 'cuffdiff', 
    padj = 0.05, lfc = 0.584, x.lim = c(-10, 10),
    title = TRUE, legend = TRUE, grid = TRUE, data.return = TRUE
  )
  ggsave(paste("VolcanoPlot2", comp_info$comp_id, "png", sep = "."),
         plot = figure$plot, width = 9, height = 6, dpi = 350, units = "in")
}
#vsVolcanoMatrix(
#  data = df_rbind_cuffdiff, d.factor = NULL, type = 'cuffdiff', 
#  padj = 0.05, x.lim = NULL, lfc = 0.584, title = TRUE, 
#  legend = TRUE, grid = TRUE, counts = TRUE
#)


### FourWay (It is usefull with venn-diagram)
figure <- vsFourWay(
  x = "EBF_P3", y = "EBF_P4", control = "YiP3_P15",
  data = df_rbind_cuffdiff, d.factor = NULL, type = 'cuffdiff', 
  padj = 0.05, lfc = 0.584, x.lim = c(-10,10), y.lim = c(-10,10),
  title = TRUE, legend = TRUE, grid = TRUE, data.return = TRUE,
)
figure




### Def()

# Correlation panel
panel.cor <- function(x, y){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y, method = "pearson"), digits=3)
  txt <- paste0(r)
  cex.cor <- 0.5/strwidth(txt)
  cex.value <- cex.cor * r
  if (cex.value < 1.0){
    cex.value <- 1.0
  }
  text(0.5, 0.5, txt, cex = cex.value, col = "red", font = 2)
}
# Customize upper panel
reg <- function(x, y, col) abline(lm(y~x), col=col) 
upper.panel <- function(x, y, col = par("col"), bg = NA, pch = par("pch"), cex = 1, col.smooth = "red", span = 2/3, iter = 3, ...)  {
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) reg(x[ok], y[ok], col.smooth)
}
# put histograms on the diagonal
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "skyblue", ...)
}
# PCA
pcaPlotting_screeplot <- function(p) {
  figure <- screeplot(p, components = getComponents(p, 1:p$n), 
                      title = '',
                      hline = 80, axisLabSize = 13, returnPlot = FALSE) +
            geom_text(aes(p$n-1, 80, label = '80% explained variation', vjust = -1), size=4)
  return(figure)
}

pcaPlotting_biplot <- function(p) {
  figure <- biplot(p, 
                   lab = p$metadata$sample, 
                   colby = 'group', 
                   hline = 0, vline = c(0),
                   ylabhjust = 0.5, ylabvjust = 0.5,
                   axisLabSize = 13,
                   labSize = 4,
                   legendPosition = 'right')
  return(figure)
}

pcaPlotting_pairsplot <- function(p) {
  max_pc <- p$n
  if (p$n > 5){
    max_pc <- 5
  }
  figure <- pairsplot(p, components = getComponents(p, c(1:max_pc)),
                      triangle = TRUE, trianglelabSize = 13,
                      hline = 0, vline = 0, hlineWidth = 0.1, vlineWidth = 0.1,
                      pointSize = 3, gridlines.major = FALSE, gridlines.minor = FALSE,
                      colby = 'group',
                      title = '', titleLabSize = 16, plotaxes = FALSE,
                      margingaps = unit(c(0.1, 0.1, 0.1, 0.1), 'cm'),
                      returnPlot = FALSE)
  return(figure)
}  

pcaPlotting_loadingsplot <- function(p) {
  if (p$n < 4){
    max_n_components <- p$n
  } else {
    max_n_components <- 4
  }
  figure <- plotloadings(p,
                         components = getComponents(p, c(1:max_n_components)),
                         rangeRetain = 0.1,
                         axisLabSize = 15,
                         labSize = 3.0,
                         absolute = FALSE,
                         title = '',
                         caption = 'Top 10% variables', captionLabSize = 13,
                         legendPosition = 'bottom', legendLabSize = 15, legendIconSize = 2.0,
                         shape = 23, shapeSizeRange = c(1, 16),
                         col = c('white', 'pink'),
                         drawConnectors = FALSE)
  return(figure)
}


vsMAPlot <- function(
    x, y, data, d.factor = NULL, type = c("cuffdiff", "deseq", "edger"),
    padj = 0.05, y.lim = NULL, lfc = NULL, title = TRUE, legend = TRUE,
    grid = TRUE, highlight = NULL, data.return = FALSE, xaxis.text.size = 10,
    yaxis.text.size = 10, xaxis.title.size = 10, yaxis.title.size = 10,
    main.title.size = 15, legend.text.size = 9
) {
  if (missing(type) || !type %in% c("cuffdiff", "deseq", "edger")) {
    stop(
      paste(
        "Please specify analysis type",
        "(\"cuffdiff\", \"deseq\", or \"edger\")"
      )
    )
  }
  
  type <- match.arg(type)
  if (type == "cuffdiff") {
    dat <- .getCuffMA(x, y, data)
  } else if (type == "deseq") {
    dat <- .getDeseqMA(x, y, data, d.factor)
  } else if (type == "edger") {
    dat <- .getEdgeMA(x, y, data)
  }
  
  if (!isTRUE(title)) {
    m.lab <- NULL
  } else {
    m.lab  <- ggtitle("MA plot", subtitle = paste(x, "vs.", y))
              
  }
  
  if (!isTRUE(legend)) {
    leg <- theme(legend.position = "none")
  } else {
    leg <- guides(
      colour = guide_legend(
        override.aes = list(size = 3)
      ),
      shape = guide_legend(
        override.aes = list(size = 3)
      )
    )
  }
  
  if (!isTRUE(grid)) {
    grid <- theme_classic()
  } else {
    grid <- theme_bw()
  }
  
  dat$isDE <- ifelse(dat$padj <= padj, TRUE, FALSE)
  py <- dat$M
  
  if (is.null(y.lim)) {
    y.lim = c(-1.5, 1.5) * quantile(abs(py[is.finite(py)]), probs = 0.99)
  }
  if (is.null(lfc)) {
    lfc = 1
  }
  
  dat <- .ma.ranker(dat, padj, lfc, y.lim)
  
  tmp.size <- .ma.out.ranker(py, y.lim[2])
  tmp.col <- .ma.col.ranker(dat$isDE, py, lfc)
  tmp.shp <- .ma.shp.ranker(py, y.lim)
  
  tmp.cnt <- .ma.col.counter(dat, lfc)
  b <- tmp.cnt[[1]]
  g <- tmp.cnt[[2]]
  
  comp1 <- .ma.comp1(y.lim, padj, lfc, b, g)
  point <- geom_point(
    alpha = 0.7,
    aes(color = tmp.col, shape = tmp.shp, size = tmp.size)
  )
  comp2 <- .ma.comp2(
    comp1[[4]], comp1[[6]], comp1[[5]], comp1[[1]], comp1[[2]], comp1[[3]]
  )
  
  text.size <- theme(
    axis.text.x = element_text(size = xaxis.text.size),
    axis.text.y = element_text(size = yaxis.text.size),
    axis.title.x = element_text(size = xaxis.title.size),
    axis.title.y = element_text(size = yaxis.title.size),
    plot.title = element_text(size = main.title.size, face = 'bold'),
    plot.subtitle = element_text(size = 10),
    legend.text = element_text(size = legend.text.size)
  )
  
  A <- id <- M <- NULL
  if (is.null(highlight)) {
    tmp.plot <- ggplot(
      dat, aes(x = A, y = pmax(y.lim[1], pmin(y.lim[2], py)))
    ) +
      point +
      comp2$color + comp2$shape + comp1$hline1 + comp1$hline2 +
      comp1$hline3 + comp1$x.lab + comp1$y.lab + m.lab + ylim(y.lim) +
      comp2$size + grid + leg + text.size
  } else {
    tl <- length(setdiff(highlight, dat$id))
    if (!is.atomic(highlight)) {
      stop("\"highlight\" must be vector.")
    } else if (all(highlight %in% dat$id)) {
      hl <- highlight
    } else if (tl > 0 && tl < length(highlight)) {
      remove <- setdiff(highlight, dat$id)
      message("Some IDs not found in data frame:")
      print(remove)
      message("Plotting the remaining samples...")
      hl <- highlight[!highlight %in% remove]
    } else if (!all(highlight %in% dat$id)) {
      stop("No IDs in highlight vector are present in data frame.")
    }
    
    A <- NULL
    tmp.plot <- ggplot(
      dat, aes(x = A, y = pmax(y.lim[1], pmin(y.lim[2], py)))
    ) +
      geom_point(
        alpha = 0.4,
        aes(color = tmp.col, shape = tmp.shp, size = tmp.size)
      ) +
      comp2$color + comp2$shape + comp1$hline1 + comp1$hline2 +
      comp1$hline3 + comp1$x.lab + comp1$y.lab + m.lab + ylim(y.lim) +
      comp2$size +
      ggrepel::geom_label_repel(
        data = dat[which(dat$id %in% hl), ],
        aes(label = id, x = A, y = M),
        segment.size = 1,
        segment.color = "gray10",
        box.padding = unit(0.4, "lines"),
        point.padding = unit(0.4, "lines")
      ) +
      geom_point(
        data = dat[which(dat$id %in% hl), ],
        aes(x = A, y = M),
        color = "red",
        size = 1
      ) +
      grid + leg + text.size
  }
  
  if (isTRUE(data.return)) {
    dat2 <- dat[, -ncol(dat)]
    plot.l <- list(data = dat2, plot = tmp.plot)
  } else {
    print(tmp.plot)
  }
}

.ma.ranker <- function(data, padj, lfc, y.lim) {
  dat <- data
  dat$color <- 'grey'
  dat$color[dat$isDE == TRUE & abs(dat$M) >= lfc] <- 'blue'
  dat$color[dat$isDE == TRUE & abs(dat$M) < lfc] <- 'green'
  dat$size <- .ma.out.ranker(dat$A, y.lim[2])
  dat$shape <- 'circle'
  dat$shape[dat$logFC < y.lim[1]] <- 'down.triangle'
  dat$shape[dat$logFC > y.lim[2]] <- 'up.triangle'
  return(dat)
}
.getCuffMA <- function(x, y, data) {
  sample_1 <- sample_2 <- NULL
  deg <- data
  deg <- subset(deg, (sample_1 == x & sample_2 == y) | 
                  (sample_1 == y & sample_2 == x))
  dat <- data.frame(id = deg$gene)
  
  if (x %in% deg$sample_1 && y %in% deg$sample_2) {
    dat$x <- deg$value_1
    dat$y <- deg$value_2
  } else if (y %in% deg$sample_1 && x %in% deg$sample_2) {
    dat$x <- deg$value_2
    dat$y <- deg$value_1
  }
  dat$x <- log10(dat$x + 1)
  dat$y <- log10(dat$y + 1)
  # dat$A <- 0.5 * log2(dat$x * dat$y)
  dat$A <- 0.5 * (dat$x + dat$y)
  dat$M <- deg[, "log2.fold_change."]
  # dat$M <- log2(dat$x / dat$y)
  dat$pval <- deg$p_value
  dat$padj <- deg$q_value
  dat <- do.call(
    data.frame, lapply(dat, function(x) replace(x, is.infinite(x), NA))
  )
  
  dat <- dat[complete.cases(dat), ]
  return(dat)
}
.ma.out.ranker <- function(log2fc, lim) {
  vec2 <- log2fc[which(abs(log2fc) >= lim)]
  tmp <- quantile(abs(vec2))
  ifelse(
    abs(log2fc) < lim, 'sub',
    ifelse(
      abs(vec2) >= tmp[[4]], 't4',
      ifelse(
        abs(vec2) >= tmp[[3]] & abs(vec2) < tmp[[4]], 't3',
        ifelse(
          abs(vec2) >= tmp[[2]] & abs(vec2) < tmp[[3]], 't2','t1'
        )
      )
    )
  ) 
}
.ma.col.ranker <- function(isDE, log2fc, lfc) {
  ifelse(
    isDE == TRUE & abs(log2fc) < lfc, 'grn',
    ifelse(isDE == TRUE & abs(log2fc) > lfc, 'blu', 'gry')
  )
}
.ma.shp.ranker <- function(log2fc, lim){
  ifelse(log2fc < lim[1], 'tri1', ifelse(log2fc > lim[2], 'tri2', 'circ'))
}
.ma.col.counter <- function(dat, lfc) {
  de <- dat$isDE
  py <- abs(dat$M)
  
  blu.c <- nrow(dat[which(py >= lfc & de == TRUE), ])
  grn.c <- nrow(dat[which(py < lfc & de == TRUE), ])
  
  l.count <- list(blu.c, grn.c)
  return(l.count)
}
.ma.comp1 <- function(y.lim, padj, lfc, b, g) {
  list(
    sh1 = paste('lfc < ', round(y.lim[1], 2)),
    sh2 = paste('lfc > ', round(y.lim[2], 2)),
    sh3 = paste(round(y.lim[1], 2), '< lfc <', round(y.lim[2], 2)),
    col1 = paste('padj >', padj),
    col2 = paste('padj <', padj, '; lfc >' , lfc, ' (', b ,')'),
    col3 = paste('padj <', padj, '; lfc <', lfc, ' (', g ,')'),
    hline1 = geom_hline(
      yintercept = 0, color = 'red3', size = 0.5,
      alpha = 0.8, linetype = 'longdash'
    ), 
    hline2 = geom_hline(
      yintercept = -lfc, 
      color = 'grey32', 
      size = 1, 
      alpha = 0.8, 
      linetype = 'dashed'
    ),
    hline3 = geom_hline(
      yintercept = lfc, 
      color = 'grey32', 
      size = 1,
      alpha = 0.8, 
      linetype = 'dashed'
    ),
    x.lab = xlab(expression(paste('Log'['10'], ' average expression'))),
    y.lab = ylab(expression(paste('Log'['2'], ' fold change')))
    
  )
}
.ma.comp2 <- function(a, b, c, d, e, f) {
  list(
    color = scale_color_manual(
      name = '', 
      values = c('gry' = 'grey73', 'grn' ='green', 'blu' ='royalblue1'), 
      labels = c('gry' = a, 'grn' = b, 'blu' = c)
    ),
    shape = scale_shape_manual(
      name = '', 
      values = c('tri1' = 6, 'tri2' = 2, 'circ' = 16),
      labels = c('tri1' = d, 'tri2' = e, 'circ' = f)
    ),
    size = scale_size_manual(
      name = '', 
      values = c(
        'sub' = 1, 
        't1' = 1.5, 
        't2' = 2, 
        't3' = 3, 
        't4' = 4
      ),
      labels = c(
        'sub' = 'SUB', 
        't1' = 'T-1', 
        't2' = 'T-2', 
        't3' = 'T-3', 
        't4' = 'T-4'
      )
    )
  )
}
