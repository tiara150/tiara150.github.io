install.packages("BiocManager")
BiocManager::install("iSEE")
library(tidyverse)
library(Biostrings)
library(biomaRt)
library(iSEE)

## The following list of commands will generate the plots created in iSEE
## Copy them into a script or an R session containing your SingleCellExperiment.
## All commands below refer to your SingleCellExperiment object as `se`.

se <- sce
colData(se)[,"sizeFactors(se)"] <- sizeFactors(se)
colormap <- ExperimentColorMap()
colormap <- synchronizeAssays(colormap, se)
all_coordinates <- list()
custom_data_fun <- NULL
custom_stat_fun <- NULL

################################################################################
## Reduced dimension plot 1
################################################################################

red.dim <- reducedDim(se, 1);
plot.data <- data.frame(X = red.dim[, 1], Y = red.dim[, 2], row.names=colnames(se));
plot.data <- subset(plot.data, !is.na(X) & !is.na(Y));

# Saving data for transmission
all_coordinates[['redDimPlot1']] <- plot.data

# Creating the plot
ggplot() +
  geom_point(aes(x = X, y = Y), alpha = 1, plot.data, color='#000000', size=1) +
  labs(x = "Dimension 1", y = "Dimension 2", title = "(1) PCA") +
  coord_cartesian(xlim = range(plot.data$X, na.rm = TRUE),
                  ylim = range(plot.data$Y, na.rm = TRUE), expand = TRUE) +
  theme_bw() +
  theme(legend.position = 'bottom', legend.box = 'vertical', legend.text=element_text(size=9), legend.title=element_text(size=11),
        axis.text=element_text(size=10), axis.title=element_text(size=12), title=element_text(size=12))

################################################################################
## Column data plot 1
################################################################################

plot.data <- data.frame(Y = colData(se)[,"NREADS"], row.names=colnames(se));
plot.data$X <- factor(character(ncol(se)))
plot.data <- subset(plot.data, !is.na(X) & !is.na(Y));

# Saving data for transmission
all_coordinates[['colDataPlot1']] <- plot.data

# Setting up plot coordinates
plot.data$GroupBy <- plot.data$X;
set.seed(100);
plot.data$jitteredX <- iSEE::jitterViolinPoints(plot.data$X, plot.data$Y, 
                                                width=0.4, varwidth=FALSE, adjust=1,
                                                method='quasirandom', nbins=NULL);

# Creating the plot
ggplot() +
  geom_violin(aes(x = X, y = Y, group = GroupBy), alpha = 0.2, data=plot.data, scale = 'width', width = 0.8) +
  geom_point(aes(y = Y, x = jitteredX), alpha = 1, plot.data, color='#000000', size=1) +
  labs(x = "", y = "NREADS", title = "NREADS ") +
  coord_cartesian(ylim = range(plot.data$Y, na.rm=TRUE), expand = TRUE) +
  scale_x_discrete(drop = FALSE) +
  theme_bw() +
  theme(legend.position = 'bottom', legend.text=element_text(size=9),
        legend.title=element_text(size=11), legend.box = 'vertical',
        axis.text.x = element_text(angle=90, size=10, hjust=1, vjust=0.5),
        axis.text.y=element_text(size=10),
        axis.title=element_text(size=12), title=element_text(size=12))

################################################################################
## Feature assay plot 1
################################################################################

plot.data <- data.frame(Y=assay(se, 6, withDimnames=FALSE)[1,], row.names = colnames(se))
plot.data$X <- factor(character(ncol(se)))
plot.data <- subset(plot.data, !is.na(X) & !is.na(Y));

# Saving data for transmission
all_coordinates[['featAssayPlot1']] <- plot.data

# Setting up plot coordinates
plot.data$GroupBy <- plot.data$X;
set.seed(100);
plot.data$jitteredX <- iSEE::jitterViolinPoints(plot.data$X, plot.data$Y, 
                                                width=0.4, varwidth=FALSE, adjust=1,
                                                method='quasirandom', nbins=NULL);

# Creating the plot
ggplot() +
  geom_violin(aes(x = X, y = Y, group = GroupBy), alpha = 0.2, data=plot.data, scale = 'width', width = 0.8) +
  geom_point(aes(y = Y, x = jitteredX), alpha = 1, plot.data, color='#000000', size=1) +
  labs(x = "", y = "0610007P14Rik (logcounts)", title = "0610007P14Rik") +
  coord_cartesian(ylim = range(plot.data$Y, na.rm=TRUE), expand = TRUE) +
  scale_x_discrete(drop = FALSE) +
  theme_bw() +
  theme(legend.position = 'bottom', legend.text=element_text(size=9),
        legend.title=element_text(size=11), legend.box = 'vertical',
        axis.text.x = element_text(angle=90, size=10, hjust=1, vjust=0.5),
        axis.text.y=element_text(size=10),
        axis.title=element_text(size=12), title=element_text(size=12))

################################################################################
## Row data plot 1
################################################################################

plot.data <- data.frame(Y = rowData(se)[,"ave_count"], row.names=rownames(se));
plot.data$X <- factor(character(nrow(se)))
plot.data <- subset(plot.data, !is.na(X) & !is.na(Y));

# Saving data for transmission
all_coordinates[['rowDataPlot1']] <- plot.data

# Setting up plot coordinates
plot.data$GroupBy <- plot.data$X;
set.seed(100);
plot.data$jitteredX <- iSEE::jitterViolinPoints(plot.data$X, plot.data$Y, 
                                                width=0.4, varwidth=FALSE, adjust=1,
                                                method='quasirandom', nbins=NULL);

# Creating the plot
ggplot() +
  geom_violin(aes(x = X, y = Y, group = GroupBy), alpha = 0.2, data=plot.data, scale = 'width', width = 0.8) +
  geom_point(aes(y = Y, x = jitteredX), alpha = 1, plot.data, color='#000000', size=1) +
  labs(x = "", y = "ave_count", title = "ave_count ") +
  coord_cartesian(ylim = range(plot.data$Y, na.rm=TRUE), expand = TRUE) +
  scale_x_discrete(drop = FALSE) +
  theme_bw() +
  theme(legend.position = 'bottom', legend.text=element_text(size=9),
        legend.title=element_text(size=11), legend.box = 'vertical',
        axis.text.x = element_text(angle=90, size=10, hjust=1, vjust=0.5),
        axis.text.y=element_text(size=10),
        axis.title=element_text(size=12), title=element_text(size=12))

################################################################################
## Sample assay plot 1
################################################################################

plot.data <- data.frame(Y=assay(se, 6, withDimnames=FALSE)[,1], row.names = rownames(se));
plot.data$X <- factor(character(nrow(se)));
plot.data <- subset(plot.data, !is.na(X) & !is.na(Y));

# Saving data for transmission
all_coordinates[['sampAssayPlot1']] <- plot.data

# Setting up plot coordinates
plot.data$GroupBy <- plot.data$X;
set.seed(100);
plot.data$jitteredX <- iSEE::jitterViolinPoints(plot.data$X, plot.data$Y, 
                                                width=0.4, varwidth=FALSE, adjust=1,
                                                method='quasirandom', nbins=NULL);

# Creating the plot
ggplot() +
  geom_violin(aes(x = X, y = Y, group = GroupBy), alpha = 0.2, data=plot.data, scale = 'width', width = 0.8) +
  geom_point(aes(y = Y, x = jitteredX), alpha = 1, plot.data, color='#000000', size=1) +
  labs(x = "", y = "SRR2140028 (logcounts)", title = "SRR2140028") +
  coord_cartesian(ylim = range(plot.data$Y, na.rm=TRUE), expand = TRUE) +
  scale_x_discrete(drop = FALSE) +
  theme_bw() +
  theme(legend.position = 'bottom', legend.text=element_text(size=9),
        legend.title=element_text(size=11), legend.box = 'vertical',
        axis.text.x = element_text(angle=90, size=10, hjust=1, vjust=0.5),
        axis.text.y=element_text(size=10),
        axis.title=element_text(size=12), title=element_text(size=12))

################################################################################
## Heat map 1
################################################################################

value.mat <- as.matrix(assay(se, 6)[1L, , drop=FALSE]);
plot.data <- reshape2::melt(value.mat, varnames = c('Y', 'X'));

plot.data[['OrderBy1']] <- colData(se)[['NREADS']][match(plot.data$X, rownames(colData(se)))];
plot.data <- dplyr::arrange(plot.data, OrderBy1);
plot.data$X <- factor(plot.data$X, levels = unique(plot.data$X));

# Centering and scaling
plot.data$value <- plot.data$value - ave(plot.data$value, plot.data$Y);

# Creating the heat map
p0 <- ggplot(plot.data, aes(x = X, y = Y)) +
  geom_raster(aes(fill = value)) +
  labs(x='', y='') +
  scale_fill_gradientn(colors=c('purple','black','yellow'),
                       values=c(0,0.70731046945827,1),
                       limits=c(-7.94193346326628,3.2864221262516), na.value='grey50') +
  scale_y_discrete(expand=c(0, 0)) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.line=element_blank());
heatlegend <- cowplot::get_legend(p0 + theme(legend.position='bottom'));

# Adding annotations
legends <- list()

p1 <- ggplot(plot.data, aes(x = X, y = 1)) +
  geom_raster(aes(fill = OrderBy1)) +
  labs(x='', y='') +
  scale_y_continuous(breaks=1, labels='NREADS') +
  scale_fill_gradientn(colors=colDataColorMap(colormap, 'NREADS', discrete=FALSE)(21L), na.value='grey50', name='NREADS') +
  theme(axis.text.x=element_blank(), axis.ticks=element_blank(), axis.title.x=element_blank(),
        rect=element_blank(), line=element_blank(), axis.title.y=element_blank(),
        plot.margin = unit(c(0,0,-0.5,0), 'lines'));
legends[[1]] <- cowplot::get_legend(p1 + theme(legend.position='bottom', plot.margin = unit(c(0,0,0,0), 'lines')));

# Laying out the grid
cowplot::plot_grid(
  cowplot::plot_grid(
    p1 + theme(legend.position='none'),
    p0 + theme(legend.position='none'),
    ncol=1, align='v', rel_heights=c(0.1, 1)),
  heatlegend, ncol=1, rel_heights=c(0.9, 0.1))

################################################################################
## To guarantee the reproducibility of your code, you should also
## record the output of sessionInfo()
sessionInfo()