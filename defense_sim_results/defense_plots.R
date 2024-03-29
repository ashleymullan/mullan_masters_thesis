
library(tidyverse)
library(latex2exp)
setwd("/Users/ashleymullan/Documents/Grad School/Wake Forest/M.S. Coursework/Research/Food-Access/masters_thesis/defense_sim_results/")

#Read in Results
vn <- read_csv("vnnz.csv")
vq <- read_csv("vqnz.csv")
vb0 <- read_csv("vb0nz.csv")
vb1 <- read_csv("vb1nz.csv")
vpr <- read_csv("vprnz.csv")

#Pivot results for plotting
#write pivot function
pivot <- function(df){
  df_long <- df |>
    pivot_longer(cols = ends_with("beta1"), 
                 names_to = "method_type", 
                 values_to = "betahat1")
  return(df_long)
}
#pivot the data
vn_long <- pivot(vn)
vq_long <- pivot(vq)
vb0_long <- pivot(vb0)
vb1_long <- pivot(vb1)
vpr_long <- pivot(vpr)

#add the error level for that one study
vpr_long <- vpr_long |> 
  mutate(error = case_when(
    tpr == 0.9 ~ "S",
    tpr == 0.75 ~ "M",
    tpr == 0.5 ~ "L"
  )) |>
  mutate(error = factor(error, levels = c("S", "M", "L"))) 

#define base plot function
plot_base <- function(df_long, x_name, xlab, truth = FALSE){
  df <- data.frame(x = df_long |> pull(x_name) |> as.factor(), 
                   y = df_long |> pull(betahat1),
                   n = df_long |> pull("n") |> as.character(),
                   method_type = case_when(
                     df_long$method_type == "gs_beta1" ~ "Gold Standard", 
                     df_long$method_type == "cc_beta1" ~ "Complete Case",
                     df_long$method_type == "our_beta1" ~ "MLE",
                     df_long$method_type == "naive_beta1" ~ "Naive"))
  df <- df |> 
    filter(abs(y) < 3) |>
    mutate(method_type = factor(method_type, 
                                levels = c("Naive", "Gold Standard",
                                           "MLE", "Complete Case")),
           n = factor(n, levels = c("390", "2200"),
                      labels = c("N = 390", "N = 2200")))
  p <- df |> ggplot(aes(x = x, y = y)) +
    geom_boxplot(aes(fill = method_type)) +
    facet_wrap(vars(n)) +
    theme_minimal() +
    labs(x = xlab,
         y = TeX("$\\hat{\\beta_1}$"),
         title = "",
         fill = "Method") +
    scale_fill_colorblind() +
    theme(legend.position = "top")
  #add a reference line for all necessary cases (fixed prevalence ratio)
  if(truth) {p <- p + geom_hline(yintercept = 0.155,
                                 color = "gray",
                                 linetype = "dashed")}
  
  #little tweaks for specific axis labeling
  if(x_name == "q") {p <- p + scale_x_discrete(labels = c('25%','50%','75%','90%'))}
  if(x_name == "b0") {p <- p + scale_x_discrete(labels = c('5%', '10%', '35%'))}
  if(x_name == "b1") {p <- p + scale_x_discrete(labels = c('1.1', '1.17', '1.25'))}
  return(p)
}

#run the plots
vq_plot <- plot_base(vq_long, "q", "Query Percentage", TRUE) #49 of 64K dropped
vpr_plot <- plot_base(vpr_long, "error", "Error Setting", TRUE) #1 of 96K dropped
vb0_plot <- plot_base(vb0_long, "b0", "Prevalence", TRUE) #2 of 12K dropped
vb1_plot <- plot_base(vb1_long, "b1", "Prevalence Ratio") #none dropped

#save the plots
ggsave(filename = "vq_plot.pdf", 
       plot = vq_plot,
       width = 5,
       height = 3.5,
       units = "in")
ggsave(filename = "vpr_plot.pdf", 
       plot = vpr_plot,
       width = 5,
       height = 3.5,
       units = "in")
ggsave(filename = "vb0_plot.pdf", 
       plot = vb0_plot,
       width = 5,
       height = 3.5,
       units = "in")
ggsave(filename = "vb1_plot.pdf", 
       plot = vb1_plot,
       width = 5,
       height = 3.5,
       units = "in")

