suppressMessages({library(ggplot2)})
plot_corr <- function(dt, out_path) { p <- ggplot(dt, aes(x=Beta,y=TPM)) + geom_point(alpha=0.7, size=1.4) + geom_smooth(method="lm", se=FALSE, size=0.8); ggsave(out_path, p, width=6, height=4.5, dpi=300, bg="white") }
