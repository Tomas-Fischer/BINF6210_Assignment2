##***************************
##
## Assignment 1
##
## Group 5
##
## 2026-11-02
##
##***************************

## _ Packages used -------
library(tidyverse)
library(ggplot2)
conflicted::conflicts_prefer(dplyr::filter())
conflicted::conflicts_prefer(dplyr::rename())
library(viridis)
theme_set(theme_light())

#####################
##### Setup ----#####
#####################

#Import the csv files
df_spiralis_te = read.csv("../BINF6210_Assignment2/2026_11_02_Assignment1_Group5/data/GCA_001447595.1.csv")
df_nelsoni_te = read.csv("../BINF6210_Assignment2/2026_11_02_Assignment1_Group5/data/GCA_001447455.1.csv")
## edited the directory

#Filter out the subsets of Retrotransposons from transposons
ls_retro_te = c("LINE/CR1", "LINE/RTE-BovB", "LTR/DIRS", "LTR/Gypsy", "LTR/Ngaro", "LTR/Pao", "PLE", "PLE/Naiad", "SINE/tRNA-Deu-L2", "SINE/tRNA-RTE")
ls_transpo_te = c("DNA/CMC-Chapaev", "DNA/hAT", "DNA/hAT-Ac", "DNA/MULE-MuDR", "DNA/PIF-Harbinger", "DNA/TcMar-Mariner", "DNA/TcMar-Pogo", "DNA/TcMar-Tc1", "DNA/TcMar-Tc4", "RC/Helitron")
ls_retro_subset = c("PLE","PLE/Naiad","LTR/Pao", "LINE/CR1", "LTR/Gypsy")

#############################
#### Data Wrangling ----#####
#############################
# EDIT 1 ---- 
## including a data frame or tibble with all of the information for both species
## having the summary is useful, but having a large df can be useful if you need to do further analysis
spiralis <- df_spiralis_te |>
  mutate(species = "spiralis")
nelsoni <- df_nelsoni_te |>
  mutate(species = "nelsoni")
df_te_all <- rbind(spiralis, nelsoni)
df_te_all <- df_te_all |> 
  mutate(
    te_type = if_else(
      matching_class %in% ls_transpo_te, "DNA transposon",
      if_else(matching_class %in% ls_retro_te, "Retrotransposon","Other")))
View(df_te_all)
# large df that can be used for analyses on other variables, and includes columns with species and te_type

#Summary of percentages for fig 1
df_spiralis_summary <- df_spiralis_te |>
  mutate(
    te_type = case_when(
      matching_class %in% ls_transpo_te ~ "DNA transposon",
      matching_class %in% ls_retro_te ~ "Retrotransposon",
      TRUE ~ "Other"
    )
  ) |>
  filter(te_type != "Other") |>
  group_by(te_type) |>
  summarise(n_copies = n(), total_bp = sum(qry_width)) |>
  mutate(
    pct_copies = 100 * n_copies / sum(n_copies),
    pct_bp = 100 * total_bp / sum(total_bp),
    species = "spiralis"
  )
df_spiralis_summary
# useful to have these here to check if your piping is working correct

# Repeating above coding steps for T. nelsoni

df_nelsoni_summary <- df_nelsoni_te |>
  mutate(
    te_type = case_when(
      matching_class %in% ls_transpo_te ~ "DNA transposon",
      matching_class %in% ls_retro_te ~ "Retrotransposon",
      TRUE ~ "Other"
    )
  ) |>
  filter(te_type != "Other") |>
  group_by(te_type) |>
  summarise(n_copies = n(), total_bp = sum(qry_width)) |>
  mutate(
    pct_copies = 100 * n_copies / sum(n_copies),
    pct_bp = 100 * total_bp / sum(total_bp),
    species = "nelsoni"
  )
df_nelsoni_summary
# useful to have as a check

df_te_summary <- rbind(df_spiralis_summary, df_nelsoni_summary)

# df_te_summary

#Filter the retrotransposon data for Figure 2

df_spiralis_retro_te = df_spiralis_te |> 
  filter(matching_class %in% ls_retro_te) |> 
  group_by(matching_class) |> 
  count() |> 
  arrange(desc(n))


df_nelsoni_retro_te = df_nelsoni_te |> 
  filter(matching_class %in% ls_retro_te) |> 
  group_by(matching_class) |> 
  count() |> 
  arrange(desc(n)) 

df_retros_long =
  left_join(df_nelsoni_retro_te, df_spiralis_retro_te, by="matching_class") |>
  replace_na(list (n.y = 0)) |>
  rename("T. nelsoni" = "n.x", "T. spiralis" = "n.y") |> 
  pivot_longer(cols = 2:3, names_to = "Species", values_to = "Count")

# _Example of EDIT 1 ----
# for example, using the large df to do this step
retros_long <- df_te_all |>
  filter(matching_class %in% ls_retro_te) |>
  group_by(matching_class, species) |>
  count() |>
  arrange(desc(n))
retros_long

#Figure 3 data wrangling

df_spiralis_div <- df_spiralis_te |>
  filter(matching_class %in% ls_retro_subset) |>
  mutate(
    species = "T. spiralis"
  ) 
df_nelsoni_div <- df_nelsoni_te |>
  filter(matching_class %in% ls_retro_subset) |>
  mutate(
    species = "T. nelsoni"
  ) 

df_trich_div <- rbind(df_spiralis_div, df_nelsoni_div)

count(df_trich_div, species, matching_class)

####################
##### Plots ----####
####################

#Figure 1 - DNA Transposon vs Retrotransposon makeup between T. nelsoni and T. spiralis

ggplot(df_te_summary, aes(x = species, y = pct_copies, fill = te_type)) +
  geom_col() +
  labs(x = NULL, y = "% of TE copies", fill = NULL,
       title = "DNA transposons vs. Retrotransposons in Trichinella species") +
  scale_fill_viridis_d(option = "viridis") +
  theme_minimal()
ggsave("../output/2026-11-01_Fig1_retro_transpo_comparison.png", width = 6, dpi=600) 

#Figure 2 2 - Retrotransposon class composition between T. nelsoni and T. spiralis
ggplot(df_retros_long, aes(x = matching_class, y = Count, fill = `Species`)) +
  geom_col(position = "dodge") +
  scale_fill_viridis_d(option = "viridis") +
  labs(title = "Retrotransposon composition between T. spiralis and T. nelsoni", x = "Retrotransposon class", y = "Number of Retrotransposons") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
ggsave("../output/2026-11-01_Fig2_retro_comparison.png", width = 8, dpi=600) 

#Figure 3 - Percent divergence between certain retrotransposons between T. nelsoni and T. spiralis
ggplot(df_trich_div, aes(x = perc_div, fill = species)) +
  geom_histogram(binwidth = 1, position = "identity", alpha = 0.5) +
  facet_wrap(~matching_class) +
  labs(x = "% divergence from consensus", y = "No. of copies", fill = NULL,
       title = "Repeat landscape for select retrotransposons in Trichinella species") +
  scale_fill_viridis_d(option = "viridis") +
  theme_dark()
ggsave("../output/2026-11-01_Fig3_retro_percdel_comparison.png", width = 8, dpi=600) 

#Mac users had to change permissions to Read and Write instead of Read only. Some of our group members had this error, but you probably already know this, but wanted to mention this just in case!

# EDIT 2 ----
# adding a figure relating to the hypothesis that spiralis will have lower divergence
ggplot(df_trich_div, aes(x = matching_class, y = perc_div, fill = species)) +
  geom_boxplot(colour = "black", linewidth = 0.5) +
  scale_fill_viridis_d(option = "viridis") +
  theme_minimal() +
  labs(x = "Retrotransposon class", y = "% divergence from consensus", fill = NULL,
       title = "Average percentage of divergence for select
       retrotransposons in Trichinella species")
# we can see that there does not appear to be any significant differences between the two species, but that the percentage divergence is higher for LTR/Gypsy for both species

#####################
##### Stats-----#####
#####################

# Chi squared test for Fig 1
counts <- matrix(df_te_summary$n_copies, nrow = 2, byrow = TRUE)
#counts
chisq.test(counts)

# X-squared = 0.00018247, df = 1,
# p-value = 0.9892

# EDIT 3 ----
## a statistical test relating to the figure I added
## testing whether average % divergence is different across species
div_lm <- lm(perc_div ~ 1 + species + te_type, data = df_te_all)
div_lm
div_anova <- aov(div_lm)
summary(div_anova)
# tells us that species and te_type individually have a significant effect on perc_div
# however, this is not useful for testing the hypothesis whether spiralis would have less diverged retrotransposons

# install.packages("emmeans")
library(emmeans)

div_emm <- emmeans(div_lm, ~ species + te_type)
div_emm
pairs(div_emm)
# just a further test from above
# gives us the pairs, and tells us p-values for species and te_type

# now testing if retrotransposon have lower divergence across species
retro_df <- df_te_all |>
  filter(te_type == "Retrotransposon")

div_lm_retro <- lm(perc_div ~ 1 + species, data = retro_df)
div_lm_retro
div_anova_retro <- aov(div_lm_retro)
summary(div_anova_retro)
# now we are getting a p-value > 0.05, indicating there is no significant difference in retrotransposon percentage divergence across species

div_emm_retro <- emmeans(div_lm_retro, ~ species)
div_emm_retro
pairs(div_emm_retro)
# this is confirming that there is no significant difference between species for retrotranposons in general

div_lm_retro2 <- lm(perc_div ~ 1 + species + matching_class, data = retro_df)
div_lm_retro2
div_anova_retro2 <- aov(div_lm_retro2)
summary(div_anova_retro2)
# now we are getting a p-value < 0.05 across matching class, significantly different
# across species, there is no significant effect on perc_div

div_emm_retro2 <- emmeans(div_lm_retro2, ~ species + matching_class)
div_emm_retro2
pairs(div_emm_retro2)
# tells us there is no significant differences in any of the same retrotransposon class comparing the two species
# there are some significant differences, but they are no meaningful for this analysis
# e.g. in nelsoni, LTR/Gypsy is significantly different from LTR/Pao, PLE, PLE/Naiad, SINE/tRNA-Deu-L2, but this does not relate to the given research question

#KS Test for Fig 3
spiralis_cr1 <- df_spiralis_div |>
  filter(matching_class == "LINE/CR1")

nelsoni_cr1 <- df_nelsoni_div |>
  filter(matching_class == "LINE/CR1")

ks_cr1 <- ks.test(spiralis_cr1$perc_div, nelsoni_cr1$perc_div)
rm(spiralis_cr1, nelsoni_cr1)

# LTR/Pao class

spiralis_pao <- df_spiralis_div |>
  filter(matching_class == "LTR/Pao")

nelsoni_pao <- df_nelsoni_div |>
  filter(matching_class == "LTR/Pao")

ks_pao <- ks.test(spiralis_pao$perc_div, nelsoni_pao$perc_div)
rm(spiralis_pao, nelsoni_pao)

# PLE class

spiralis_ple <- df_spiralis_div |>
  filter(matching_class == "PLE")

nelsoni_ple <- df_nelsoni_div |>
  filter(matching_class == "PLE")

ks_ple <- ks.test(spiralis_ple$perc_div, nelsoni_ple$perc_div)
rm(spiralis_ple, nelsoni_ple)

# PLE/Naiad class

spiralis_nai <- df_spiralis_div |>
  filter(matching_class == "PLE/Naiad")

nelsoni_nai <- df_nelsoni_div |>
  filter(matching_class == "PLE/Naiad")

ks_nai <- ks.test(spiralis_nai$perc_div, nelsoni_nai$perc_div)
rm(spiralis_nai, nelsoni_nai)

# LTR/Gypsy class

spiralis_gypsy <- df_spiralis_div |>
  filter(matching_class == "LTR/Gypsy")

nelsoni_gypsy <- df_nelsoni_div |>
  filter(matching_class == "LTR/Gypsy")

ks_gypsy <- ks.test(spiralis_gypsy$perc_div, nelsoni_gypsy$perc_div)
ks_gypsy
rm(spiralis_gypsy, nelsoni_gypsy)

# Adjusting p-values to ensure there are no false positives by multiplying by the number of KS tests performed (5)

ls_ks_pvalues <- c("LINE/CR1" = ks_cr1$p.value, "LTR/Pao" = ks_pao$p.value, "PLE" = ks_ple$p.value, "PLE/Naiad" = ks_nai$p.value, "LTR/Gypsy" = ks_gypsy$p.value)

p.adjust(ls_ks_pvalues, method = "bonferroni")

# P values:
# LINE/CR1   LTR/Pao       PLE PLE/Naiad 
# 1.0000000 0.2088786 1.0000000 1.0000000 
# LTR/Gypsy 
# 1.0000000 

rm(ks_cr1, ks_nai, ks_pao, ks_ple, ks_gypsy)