library(ggplot2)
library(dplyr)
library(splines)
library(ggformula)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
	stop("Usage: Rscript depth_plot.R <input_directory> <output_file1 (depth plot)> <output_file2 (zscore plot)  <windowsize>")
}
inputdir <- args[1]
out1 <- args[2]
out2 <- args[3]
windowsize <- args[4]

pattern <- paste0('_windowed_depth\\.')
depth_files <- list.files(inputdir, pattern = pattern, full.names = TRUE)

depth_data <- lapply(depth_files, function(file) {
	
	name <- strsplit(basename(file), '_window')[[1]][1]
	data <- read.table(file, header = FALSE, col.names = c('chr', 'start', 'end', 'depth'))
	data$position <- (data$start+ data$end) / 2
	data$pos_label <- paste(data$start/1000, 'kbp', '-', (data$end+1)/1000, 'kbp', sep='')
	data$sample <- name

	data$zscore <- (data$depth - mean(data$depth, na.rm = TRUE)) / sd(data$depth, na.rm = TRUE)

	return(data)
})

depth_data <- bind_rows(depth_data)
head(depth_data)
depth_data <- depth_data %>% filter(chr %in% c(as.character(1:19), 'X', 'Y'))

# Detph plot
p <- ggplot(depth_data, aes(x = position, y = depth)) +
	geom_point(aes(fill = factor(sample)), size = 0.4, alpha = 0.1, shape = 21, stroke = NA)+
	scale_y_continuous(limits = c(0,30))+
	labs(title = "Genome-wide Depth by Sample",
		x = paste0("Genomic Position (window size:", windowsize, ")"),
		y = "Depth",
		color = "Sample"
	)+
	scale_x_continuous(labels = label_number(suffix = " k", scale = 1e-6))+
	facet_grid(~ factor(chr, levels = unique(depth_data$chr)), scales = 'free', space = 'free_x')+
	theme_bw()+
	geom_spline(spar = 0.2, aes(color = sample), linewidth = 0.2)+
	guides(fill = 'none')+
	theme(axis.text.x = element_text(size = 15, angle = 45, hjust = 1))

ggsave(out1, plot = p, width = 45, height = 4, dpi = 300, device = 'png')

# Z-score plot
p <- ggplot(depth_data, aes(x = position, y = zscore)) +
	geom_point(aes(fill = factor(sample)), size = 0.4, alpha = 0.1, shape = 21, stroke = NA)+
	scale_y_continuous(limits = c(-1.5, 1.5))+
	labs(title = "Genome-wide Depth by Sample",
		x = paste0("Genomic Position (window size:", windowsize, ")"),
		y = 'Z-Score',
		color = 'Sample'
	)+
	scale_x_continuous(labels = label_number(suffix = " k", scale = 1e-6))+
	facet_grid(~ factor(chr, levels = unique(depth_data$chr)), scales = 'free', space = 'free_x',)+
	theme_bw()+
	geom_spline(spar = 0.2, aes(color = sample), linewidth = 0.2)+
	guides(fill = 'none')+
	theme(axis.text.x = element_text(size = 15), angle = 45, hjust = 1)

ggsave(out2, plot = p, width = 45, height = 4, dpi = 300, device = 'png')

