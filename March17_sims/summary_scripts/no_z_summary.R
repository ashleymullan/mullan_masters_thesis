library(dplyr)
library(readr)
library(kableExtra)

setwd("/deac/sta/lotspeichGrp/mullae22/thesis/Sim_Results")
vnnz <- read_csv("vnnz.csv") 
vnnz <- vnnz |>
  cbind(fpr = rep(0.1, times = nrow(vnnz)),
        tpr = rep(0.9, times = nrow(vnnz)),
        b0 = rep(-2.246, times = nrow(vnnz)),
        b1 = rep(0.155, times = nrow(vnnz)))

vqnz <- read_csv("vqnz.csv")
vqnz <- vqnz |>
  cbind(fpr = rep(0.1, times = nrow(vqnz)),
        tpr = rep(0.9, times = nrow(vqnz)),
        b0 = rep(-2.246, times = nrow(vqnz)),
        b1 = rep(0.155, times = nrow(vqnz)))

vprnz <- read_csv("vprnz.csv") 
vprnz <- vprnz |>
  cbind(b0 = rep(-2.246, times = nrow(vprnz)),
        b1 = rep(0.155, times = nrow(vprnz)))

vb0nz <- read_csv("vb0nz.csv") 
vb0nz <- vb0nz |>
  cbind(fpr = rep(0.1, times = nrow(vb0nz)),
        tpr = rep(0.9, times = nrow(vb0nz)),
        b1 = rep(0.155, times = nrow(vb0nz)))

vb1nz <- read_csv("vb1nz.csv") 
vb1nz <- vb1nz |>
  cbind(fpr = rep(0.1, times = nrow(vb1nz)),
        tpr = rep(0.9, times = nrow(vb1nz)),
        b0 = rep(-2.246, times = nrow(vb1nz)))

full_results <- rbind(vnnz, vqnz, vprnz, vb0nz, vb1nz)
full_results |> nrow() #94000 reps
full_results |> drop_na() |> nrow() #89890 reps, 4110 dropped

#compute coverage proportion
cp = function(est, se, truth) {
  mean((est - 1.96 * se) <= truth & truth <= (est + 1.96 * se))
}

full_result_summary <- full_results |> 
  drop_na() |>
  group_by(n,q,fpr,tpr,b0,b1) |> 
  summarize(bias_gs = mean((gs_beta1 - b1) / b1), ese_gs = sd(gs_beta1), 
            bias_n = mean((naive_beta1 - b1) / b1), ese_n = sd(naive_beta1), 
            bias_cc = mean((cc_beta1 - b1) / b1), ese_cc = sd(cc_beta1), 
            bias_mle = mean((our_beta1 - b1) / b1), ese_mle = sd(our_beta1), 
            ase_mle = mean(our_beta1_se), 
            cp_mle = cp(est = our_beta1, se = our_beta1_se, truth = b1)
  ) |> 
  dplyr::mutate(re_cc = (ese_gs ^ 2) / (ese_cc ^ 2), 
                re_mle = (ese_gs ^ 2) / (ese_mle ^ 2)) 

#function to format numbers for LaTeX
format_num = function(num, digits = 3) {
  paste0("$", format(round(num, 3), nsmall = digits), "$")
}

## Format for LaTeX
full_result_summary = full_result_summary |> 
  mutate_at(.vars = 5:18, .funs = format_num, digits = 3) |>
  mutate_at(.vars = 2:4, .funs = format_num, digits = 2)

#change col names
colnames(full_result_summary) = c("$\\pmb{N}$", "$\\pmb{Q}$", 
                                  "FPR", "TPR", 
                                  "\\pmb{\\beta_0}", "\\pmb{\\beta_1}",
                       rep(c("Bias", "ESE"), times = 4),
                       "ASE","CP", "RE", "RE")
#GGNNCCMM MMCM

full_result_summary |> 
  kable(format = "latex", 
        booktabs = TRUE, 
        escape = FALSE, 
        align = "cccccccccccccccccc") |> 
  kable_styling() |> 
  add_header_above(header = c(" " = 6, "Gold Standard" = 2, 
                              "Naive" = 2, "Complete Case" = 2, 
                              "MLE" = 4, "Complete Case" = 1, "MLE" = 1), 
                   bold = TRUE) |> 
  row_spec(row = 0, bold = TRUE)
## And a \multicolumn used to separate the three parameters

