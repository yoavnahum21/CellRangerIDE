# library(deMULTIplex2)
library(dplyr)
library(ggplot2)
# f"log_file <- '{os.path.join(deMULTIplex2_dir, 'deMULTIplex_log.txt')}'"
# "sink(log_file, append = TRUE)" # Redirect output to log file

cat('Starting demultiplexing process...\n')
mat <- read.csv(csv_path)
rownames(mat) <- mat$X
mat$X <- NULL
colnames(mat) <- gsub('X', '', colnames(mat))
res <- demultiplexTags(mat, plot.path = '{deMULTIplex2_dir}', plot.name = 'demultiplex_plot', plot.diagnostics = TRUE)
prob_mtx_df <- as.data.frame(res$prob_mtx)
res_mtx_df <- as.data.frame(res$res_mtx)
final_assign_df <- as.data.frame(res$final_assign)
write.csv(prob_mtx_df, '{os.path.join(deMULTIplex2_dir, 'prob_mtx_df_deMULTIplex2.csv')}')
write.csv(res_mtx_df, '{os.path.join(deMULTIplex2_dir, 'res_mtx_df_deMULTIplex2.csv')}')
write.csv(final_assign_df, '{os.path.join(deMULTIplex2_dir, 'final_assign_deMULTIplex2.csv')}')