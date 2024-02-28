library(tidyverse)
library(latex2exp)
vn <- read_csv("/Users/ashleymullan/Documents/Grad School/Wake Forest/M.S. Coursework/Research/Food-Access/masters_thesis/varied_n_full_sims.csv")
vq <- read_csv("/Users/ashleymullan/Documents/Grad School/Wake Forest/M.S. Coursework/Research/Food-Access/masters_thesis/varied_q_full_sims.csv")


vn_long <- vn |>
  pivot_longer(cols = ends_with("beta1"), 
               names_to = "method_type", 
               values_to = "betahat1")

vq_long <- vq |>
  pivot_longer(cols = ends_with("beta1"), 
               names_to = "method_type", 
               values_to = "betahat1")
vn_plot <- 
  vn_long |>
  mutate(method_type = case_when(
    method_type == "gs_beta1" ~ "Gold Standard", 
    method_type == "cc_beta1" ~ "Complete Case",
    method_type == "our_beta1" ~ "Ours",
    method_type == "naive_beta1" ~ "Naive")) |>
  #filter(method_type != "Naive") |> #temporarily axe naive
  ggplot(aes(x = as.factor(n), y = betahat1)) +
  geom_boxplot(aes(fill = method_type)) +
  theme_minimal() +
  scale_x_discrete(labels = c('100','1,000','10,000')) +
  theme(plot.title = element_text(size = 9.5),
        plot.subtitle = element_text(size = 9.5),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  geom_hline(yintercept = 3, color = "gray", linetype = "dashed") +
  labs(y = TeX("$\\hat{\\beta_1}$"),
       fill = "Method",
       x = "(Original) Sample Size",
       title = "Method Comparison with Fixed Data Quality",
       subtitle = "25% Missingness")
vn_plot_no_naive <- 
  vn_long |>
  mutate(method_type = case_when(
    method_type == "gs_beta1" ~ "Gold Standard", 
    method_type == "cc_beta1" ~ "Complete Case",
    method_type == "our_beta1" ~ "Ours",
    method_type == "naive_beta1" ~ "Naive")) |>
  filter(method_type != "Naive") |> #temporarily axe naive
  ggplot(aes(x = as.factor(n), y = betahat1)) +
  geom_boxplot(aes(fill = method_type)) +
  theme_minimal() +
  scale_x_discrete(labels = c('100','1,000','10,000')) +
  theme(plot.title = element_text(size = 9.5),
        plot.subtitle = element_text(size = 9.5),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  geom_hline(yintercept = 3, color = "gray", linetype = "dashed") +
  labs(y = TeX("$\\hat{\\beta_1}$"),
       fill = "Method",
       x = "(Original) Sample Size",
       title = "Method Comparison with Fixed Data Quality",
       subtitle = "25% Missingness")
ggsave(filename = "vn_plot.pdf", 
       plot = vn_plot,
       path = "/Users/ashleymullan/Documents/Grad School/Wake Forest/M.S. Coursework/Research/Food-Access/masters_thesis",
       width = 5,
       height = 3.5,
       units = "in")
ggsave(filename = "vn_plot_no_naive.pdf", 
       plot = vn_plot_no_naive,
       path = "/Users/ashleymullan/Documents/Grad School/Wake Forest/M.S. Coursework/Research/Food-Access/masters_thesis",
       width = 5,
       height = 3.5,
       units = "in")


vq_plot <- 
  vq_long |>
  mutate(method_type = case_when(
    method_type == "gs_beta1" ~ "Gold Standard", 
    method_type == "cc_beta1" ~ "Complete Case",
    method_type == "our_beta1" ~ "Ours",
    method_type == "naive_beta1" ~ "Naive")) |>
  #filter(method_type != "Naive") |> #temporarily axe naive
  ggplot(aes(x = as.factor(q), y = betahat1)) +
  geom_boxplot(aes(fill = method_type)) +
  theme_minimal() +
  scale_x_discrete(labels = c('50%','75%','90%')) +
  theme(plot.title = element_text(size = 9.5),
        plot.subtitle = element_text(size = 9.5),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  geom_hline(yintercept = 3, color = "gray", linetype = "dashed") +
  labs(y = TeX("$\\hat{\\beta_1}$"),
       fill = "Method",
       x = "Percent Queried",
       title = "Method Comparison with Fixed Data Quality",
       subtitle = "N = 1000")
vq_plot_no_naive <- 
  vq_long |>
  mutate(method_type = case_when(
    method_type == "gs_beta1" ~ "Gold Standard", 
    method_type == "cc_beta1" ~ "Complete Case",
    method_type == "our_beta1" ~ "Ours",
    method_type == "naive_beta1" ~ "Naive")) |>
  filter(method_type != "Naive") |> #temporarily axe naive
  ggplot(aes(x = as.factor(q), y = betahat1)) +
  geom_boxplot(aes(fill = method_type)) +
  theme_minimal() +
  scale_x_discrete(labels = c('50%','75%','90%')) +
  theme(plot.title = element_text(size = 9.5),
        plot.subtitle = element_text(size = 9.5),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  geom_hline(yintercept = 3, color = "gray", linetype = "dashed") +
  labs(y = TeX("$\\hat{\\beta_1}$"),
       fill = "Method",
       x = "Percent Queried",
       title = "Method Comparison with Fixed Data Quality",
       subtitle = "N = 1000")
ggsave(filename = "vq_plot.pdf", 
       plot = vq_plot,
       path = "/Users/ashleymullan/Documents/Grad School/Wake Forest/M.S. Coursework/Research/Food-Access/masters_thesis",
       width = 5,
       height = 3.5,
       units = "in")
ggsave(filename = "vq_plot_no_naive.pdf", 
       plot = vq_plot_no_naive,
       path = "/Users/ashleymullan/Documents/Grad School/Wake Forest/M.S. Coursework/Research/Food-Access/masters_thesis",
       width = 5,
       height = 3.5,
       units = "in")




