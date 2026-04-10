# Neat analysis script 

# --- Load Essential Libraries ----
library(dplyr)
library(tidyverse)
library(ggplot2)
library(lme4)
library(glmmTMB)
library(DHARMa)
library(ggpubr)
library(ggtext)
library(ggeffects)
library(RColorBrewer)
library(stringr)
library(emmeans)
library(summarytools)

# --- Load Data ----
main_data <- read.csv(
  "https://docs.google.com/spreadsheets/d/e/2PACX-1vSOE-j1tXxme4Nxmf3lw8yp7Z03wIpUOw3B51WHd0ZQ2BMNUE890hcHBaKsYw1-ESdr9Yn8jomxIh_X/pub?gid=609329912&single=true&output=csv"
)

# --- Select Relevant Columns ---
main_data <- main_data[, c(5:10, 12)]  

# --- Data Wrangling ----
data_mod <- main_data %>%
  mutate(
    Gnathopod_ratio = Gnathopod_2_length / Body_length,
    Rugosity = case_when(                      # Proxy for structural complexity
      Habitat == "B" ~ 1,                      # Black coral forests
      Habitat == "S" ~ 2,                      # Seagrass meadows
      Habitat == "M" ~ 3,                      # Macroalgae
      Habitat == "R" ~ 4                       # Rhodolith beds
    ),
    Sex = case_when(
      Female == 1 ~ 0,
      Male == 1 ~ 1
    ),
    Sex_label = case_when(
      Sex == 1 ~ "Male",
      Sex == 0 ~ "Female"
    )
  )

# --- Tidy Habitat Variable ---
data_mod$Habitat <- str_trim(data_mod$Habitat)
data_mod$Habitat <- factor(data_mod$Habitat, levels = c("R", "M", "S", "B"))

# --- Exploring dataset ---

# --- 1. Species and Sample Sizes ---
species_summary <- data_mod %>%
  group_by(Species) %>%
  summarise(
    N_total = n(),
    N_male = sum(Sex == 1, na.rm = TRUE),
    N_female = sum(Sex == 0, na.rm = TRUE)
  ) %>%
  arrange(desc(N_total))

print(species_summary)

data_mod %>%
  filter(Species == "Ampithoe ramondi") %>%
  group_by(Habitat, Sex_label) %>%
  summarise(N = n())

data_mod %>%
  filter(Species == "Caprella acanthifera") %>%
  group_by(Habitat, Sex_label) %>%
  summarise(N = n(), .groups = "drop")

# --- 2. Variables Present ---
# Check column names and types
str(data_mod)

# Quick summary stats
summary(data_mod)

# --- 3. Missing Data ---
missing_summary <- data_mod %>%
  summarise_all(~sum(is.na(.))) %>%
  gather(key = "Variable", value = "Missing_Count") %>%
  arrange(desc(Missing_Count))

print(missing_summary)

# --- 4. Sample Size Threshold ---
# Minimum N per species to include
min_N <- 5
species_include <- species_summary %>%
  filter(N_total >= min_N) %>%
  pull(Species)

cat("Species meeting minimum sample size (N >=", min_N, "):\n")
print(species_include)

# --- 5. Quick overview of sexual dimorphism across all species ---
# Means by species and sex
sex_dimorphism <- data_mod %>%
  filter(Species %in% species_include) %>%
  group_by(Species, Sex_label) %>%
  summarise(
    mean_Body_length = mean(Body_length, na.rm = TRUE),
    sd_Body_length = sd(Body_length, na.rm = TRUE),
    mean_Gnathopod_2_length = mean(Gnathopod_2_length, na.rm = TRUE),
    sd_Gnathopod_2_length = sd(Gnathopod_2_length, na.rm = TRUE),
    mean_Gnathopod_ratio = mean(Gnathopod_ratio, na.rm = TRUE),
    sd_Gnathopod_ratio = sd(Gnathopod_ratio, na.rm = TRUE),
    N = n()
  ) %>%
  arrange(Species, Sex_label)

print(sex_dimorphism)

# --- Analaysis without considering habitats ---

# --- Select species with minimum 5 individuals per sex ---
min_N_per_sex <- 5

species_include <- species_summary %>%
  filter(N_male >= min_N_per_sex & N_female >= min_N_per_sex) %>%
  pull(Species)

cat("Species meeting minimum sample size per sex (>= ", min_N_per_sex, "):\n")
print(species_include)

data_analysis <- data_mod %>%
  filter(Species %in% species_include)

# --- Step 2: Linear model for Body_length ---
model_body <- lm(Body_length ~ Species * Sex_label, data = data_analysis)
summary(model_body)

# Model diagnostics:

              # Step 1: Fit the body length model (with species and sex)
model_body <- lm(Body_length ~ Species * Sex_label, data = data_analysis)

              # Step 2: Extract residuals and fitted values
residuals_body <- residuals(model_body)
fitted_body <- fitted(model_body)

              # Step 3: Residuals vs Fitted plot
plot(fitted_body, residuals_body,
     xlab = "Fitted values",
     ylab = "Residuals",
     main = "Residuals vs Fitted - Body Length",
     pch = 19, col = "#56B4E9")
abline(h = 0, lty = 2, col = "red")

              # Step 4: Q-Q plot for residuals
qqnorm(residuals_body, main = "Q-Q Plot - Body Length Residuals")
qqline(residuals_body, col = "red")

              # Histogram of residuals
hist(residuals_body,
     breaks = 30,
     col = "#56B4E9",
     border = "white",
     main = "Histogram of Body Length Residuals",
     xlab = "Residuals",
     ylab = "Frequency")


# --- Step 3: Linear model for Gnathopod_ratio ---
model_ratio <- lm(Gnathopod_ratio ~ Species * Sex_label, data = data_analysis)
summary(model_ratio)

# Model diagnostics:

      # Step 2: Extract residuals and fitted values
residuals_ratio <- residuals(model_ratio)
fitted_ratio <- fitted(model_ratio)

      # Step 3: Residuals vs Fitted plot
plot(fitted_ratio, residuals_ratio,
     xlab = "Fitted values",
     ylab = "Residuals",
     main = "Residuals vs Fitted - Gnathopod Ratio",
     pch = 19, col = "#56B4E9")
abline(h = 0, lty = 2, col = "red")

      # Step 4: Q-Q plot for residuals
qqnorm(residuals_ratio, main = "Q-Q Plot - Body Length Residuals")
qqline(residuals_ratio, col = "red")

      # Histogram of residuals
hist(residuals_ratio,
     breaks = 30,
     col = "#56B4E9",
     border = "white",
     main = "Histogram of Gnathopod Ratio Residuals",
     xlab = "Residuals",
     ylab = "Frequency")

# PLOTS

# Define sex colors
sex_colors <- c("Female" = "#e377c2", "Male" = "#1f77b4")

# Create a new column with abbreviated species names
data_analysis$Species_abbrev <- recode(data_analysis$Species,
                                       "Gammaropsis ostromowi" = "G. ostromowi",
                                       "Pseudoprotella phasma" = "P. phasma",
                                       "Caprella acanthifera" = "C. acanthifera",
                                       "Erichtonius punctatus" = "E. punctatus",
                                       "Aora spinicornis" = "A. spinicornis",
                                       "Ampithoe ramondi" = "A.ramondi"
)

# Update abbreviated species names in your data
data_analysis$Species_abbrev <- recode(data_analysis$Species,
                                       "Gammaropsis ostromowi"   = "G. ostromowi",
                                       "Pseudoprotella phasma"   = "P. phasma",
                                       "Caprella acanthifera"    = "C. acanthifera",
                                       "Erichtonius punctatus"   = "E. punctatus",
                                       "Aora spinicornis"        = "A. spinicornis",
                                       "Ampithoe ramondi"        = "A. ramondi"
)

# Map abbreviated species names to bold + italic expressions
species_labels <- c(
  "A. spinicornis"   = "bolditalic('A. spinicornis')",
  "C. acanthifera"   = "bolditalic('C. acanthifera')",
  "E. punctatus"     = "bolditalic('E. punctatus')",
  "G. ostromowi"     = "bolditalic('G. ostromowi')",
  "P. phasma"        = "bolditalic('P. phasma')",
  "A. ramondi"       = "bolditalic('A. ramondi')"
)

# Body length plot with ultra-bold, large text
plot_body_length <- ggplot(data_analysis, aes(x = Species_abbrev, y = Body_length, fill = Sex_label)) +
  geom_boxplot(alpha = 0.6, position = position_dodge(width = 0.8), outlier.shape = NA, width = 0.6) +
  geom_jitter(aes(color = Sex_label),
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
              alpha = 0.6, size = 2.5) +
  scale_fill_manual(values = sex_colors, name = "Sex") +
  scale_color_manual(values = sex_colors, name = "Sex") +
  scale_x_discrete(labels = function(x) parse(text = species_labels[x])) +
  labs(x = NULL, y = "Body Length (mm)") +
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 36),
    legend.text  = element_text(face = "bold", size = 32),
    axis.title   = element_text(face = "bold", size = 36),
    axis.text    = element_text(face = "bold", size = 30),
    axis.text.x  = element_text(angle = 45, hjust = 1, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 1.2)
  )

# Display the plot
print(plot_body_length)

# Save high-resolution, publication-ready plot
ggsave(
  filename = "BodyLength_by_Species_Sex.png",
  plot = plot_body_length,
  width = 16,
  height = 10,
  dpi = 600,
  units = "in",
  bg = "white"
)


# Gnathopod ratio plot with bold, large text
plot_gnathopod_ratio <- ggplot(data_analysis, aes(x = Species_abbrev, y = Gnathopod_ratio, fill = Sex_label)) +
  geom_boxplot(alpha = 0.6, position = position_dodge(width = 0.8), outlier.shape = NA, width = 0.6) +
  geom_jitter(aes(color = Sex_label),
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
              alpha = 0.6, size = 2.5) +
  scale_fill_manual(values = sex_colors, name = "Sex") +
  scale_color_manual(values = sex_colors, name = "Sex") +
  scale_x_discrete(labels = function(x) parse(text = species_labels[x])) +
  labs(x = NULL, y = "Gnathopod Ratio") +
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 36),
    legend.text  = element_text(face = "bold", size = 32),
    axis.title   = element_text(face = "bold", size = 36),
    axis.text    = element_text(face = "bold", size = 30),
    axis.text.x  = element_text(angle = 45, hjust = 1, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 1.2)
  )

# Display plot
print(plot_gnathopod_ratio)

# Save high-quality, publication-ready figure
ggsave(
  filename = "GnathopodRatio_by_Species_Sex.png",
  plot = plot_gnathopod_ratio,
  width = 16,
  height = 10,
  dpi = 600,
  units = "in",
  bg = "white"
)


# ---- Sex-specific Allometry Plot with larger font and bold-italic species names ----

species_for_manova <- c(
  "Ampithoe ramondi", "Caprella acanthifera", 
  "Gammaropsis ostromowi", "Pseudoprotella phasma", 
  "Aora spinicornis", "Erichtonius punctatus"
)

# Subset and log-transform
data_allometry <- data_mod %>%
  filter(Species %in% species_for_manova) %>%
  mutate(
    log_body = log(Body_length),
    log_gnath = log(Gnathopod_2_length)
  )

# Abbreviated species labels in bold italics
species_labels <- c(
  "Ampithoe ramondi"     = "**_A. ramondi_**",
  "Caprella acanthifera" = "**_C. acanthifera_**",
  "Gammaropsis ostromowi"= "**_G. ostromowi_**",
  "Pseudoprotella phasma"= "**_P. phasma_**",
  "Aora spinicornis"     = "**_A. spinicornis_**",
  "Erichtonius punctatus"= "**_E. punctatus_**"
)

# Allometry plot - ultra-big bold
plot_allometry <- ggplot(data_allometry, aes(x = log_body, y = log_gnath, color = Sex_label, fill = Sex_label)) +
  geom_point(alpha = 0.6, size = 3) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.5, alpha = 0.3) +
  facet_wrap(~ Species, labeller = labeller(Species = species_labels)) +
  scale_color_manual(values = c("Female" = "#e377c2", "Male" = "#1f77b4"), name = "Sex") +
  scale_fill_manual(values = c("Female" = "#e377c2", "Male" = "#1f77b4"), guide = "none") +
  labs(
    x = "log(Body Length)",
    y = "log(Gnathopod 2 Length)"
  ) +
  theme_minimal(base_size = 24) +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(face = "bold", size = 40),
    axis.text = element_text(face = "bold", size = 32),
    legend.title = element_text(face = "bold", size = 36),
    legend.text = element_text(face = "bold", size = 32),
    strip.text = ggtext::element_markdown(size = 36),  # bold-italic facet labels
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 1.5)
  )

# Display
print(plot_allometry)

# Save high-resolution publication-ready plot
ggsave(
  filename = "Sex_specific_allometry.png",
  plot = plot_allometry,
  width = 20,
  height = 12,
  dpi = 600,
  units = "in",
  bg = "white"
)


six_species_sex_speicific_allometry <- lm(log(Gnathopod_2_length) ~ log(Body_length) * Sex_label * Species, data = data_allometry)
summary(six_species_sex_speicific_allometry)

# For A.ramondi and C.acanthifera ----

# ---- Sex-specific Allometry Plot for Ampithoe ramondi and Caprella acanthifera ----

species_for_manova <- c("Ampithoe ramondi", "Caprella acanthifera")

# Subset and log-transform
data_allometry <- data_mod %>%
  filter(Species %in% species_for_manova) %>%
  mutate(
    log_body = log(Body_length),
    log_gnath = log(Gnathopod_2_length)
  )

# Abbreviated species labels in bold italics
species_labels <- c(
  "Ampithoe ramondi"     = "**_A. ramondi_**",
  "Caprella acanthifera" = "**_C. acanthifera_**"
)

# Allometry plot - ultra-big bold
plot_allometry <- ggplot(data_allometry, aes(x = log_body, y = log_gnath, color = Sex_label, fill = Sex_label)) +
  geom_point(alpha = 0.6, size = 3) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.5, alpha = 0.3) +
  facet_wrap(~ Species, labeller = labeller(Species = species_labels)) +
  scale_color_manual(values = c("Female" = "#e377c2", "Male" = "#1f77b4"), name = "Sex") +
  scale_fill_manual(values = c("Female" = "#e377c2", "Male" = "#1f77b4"), guide = "none") +
  labs(
    x = "log(Body Length)",
    y = "log(Gnathopod 2 Length)"
  ) +
  theme_minimal(base_size = 24) +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(face = "bold", size = 40),
    axis.text = element_text(face = "bold", size = 32),
    legend.title = element_text(face = "bold", size = 36),
    legend.text = element_text(face = "bold", size = 32),
    strip.text = ggtext::element_markdown(size = 36),  # bold-italic facet labels
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 1.5)
  )

# Display
print(plot_allometry)

# Save high-resolution publication-ready plot
ggsave(
  filename = "Sex_specific_allometry_Ar_Ca.png",
  plot = plot_allometry,
  width = 20,
  height = 12,
  dpi = 600,
  units = "in",
  bg = "white"
)


# --- Subsets for Selected Species ---
data_ramondi <- subset(data_mod, Species == "Ampithoe ramondi")
data_acanthifera <- subset(data_mod, Species == "Caprella acanthifera")

# --- Optional: Clean Habitat Levels for Subsets Too ---
data_acanthifera$Habitat <- factor(str_trim(data_acanthifera$Habitat), levels = c("R", "M", "S", "B"))

# --- Sample Size Overview (Descriptive Table) ---
sample_summary <- data_mod %>%
  filter(Species %in% c("Caprella acanthifera", "Ampithoe ramondi")) %>%
  group_by(Species, Habitat, Sex_label) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(Species, Habitat, Sex_label)

print(sample_summary)

# --- Pivoted Summary (Optional for Supplementary Info Table) ---
sample_summary_wide <- sample_summary %>%
  tidyr::pivot_wider(names_from = Sex_label, values_from = n, values_fill = 0)

print(sample_summary_wide)

# --- Sexual Dimorphism Index (SDI) ----

# Focused species (where sample sizes are sufficient for both sexes)
species_focus <- c("Ampithoe ramondi", "Caprella acanthifera")

# Function to compute SDI for a trait
compute_sdi <- function(df, trait) {
  male_mean   <- mean(df[[trait]][df$Sex_label == "Male"], na.rm = TRUE)
  female_mean <- mean(df[[trait]][df$Sex_label == "Female"], na.rm = TRUE)
  
  # Return NA if invalid (zero or missing denominator)
  if (is.na(male_mean) || is.na(female_mean) || female_mean == 0) return(NA)
  
  # SDI formula: (Male - Female) / Female
  return((male_mean - female_mean) / female_mean)
}

# Filter dataset to the species of interest
data_filtered <- data_mod %>%
  filter(Species %in% species_focus)

# Count samples per Species × Habitat × Sex
sample_counts <- data_filtered %>%
  group_by(Species, Habitat, Sex_label) %>%
  summarise(n = n(), .groups = "drop")

# Identify valid groups (at least 5 individuals per sex)
valid_habitats <- sample_counts %>%
  pivot_wider(names_from = Sex_label, values_from = n, values_fill = 0) %>%
  filter(Male >= 5 & Female >= 5)

# Retain only data from valid groups
valid_data <- data_filtered %>%
  semi_join(valid_habitats, by = c("Species", "Habitat"))

# Compute SDI for each trait across valid groups
sdi_results <- valid_data %>%
  group_by(Species, Habitat) %>%
  summarise(
    SDI_Body_length        = compute_sdi(cur_data(), "Body_length"),
    SDI_Gnathopod_2_length = compute_sdi(cur_data(), "Gnathopod_2_length"),
    SDI_Gnathopod_ratio    = compute_sdi(cur_data(), "Gnathopod_ratio"),
    .groups = "drop"
  )

# View SDI table
print(sdi_results)

# Optional: Save SDI results to CSV
write.csv(sdi_results, "sexual_dimorphism_indices.csv", row.names = FALSE)

# ---- SDI Plot by Species and Habitat ----

# Clean and recode habitat labels
sdi_results <- sdi_results %>%
  mutate(
    Habitat = str_trim(gsub('"', '', Habitat)),
    Habitat = case_when(
      Habitat == "R" ~ "Rhodolith beds",
      Habitat %in% c("M", "M ") ~ "Macroalgae",
      Habitat == "S" ~ "Seagrass meadows",
      Habitat == "B" ~ "Black coral forests",
      TRUE ~ Habitat
    ),
    Habitat = factor(Habitat, levels = c("Rhodolith beds", "Macroalgae", "Seagrass meadows", "Black coral forests"))
  )

# Pivot longer for plotting
sdi_long <- sdi_results %>%
  pivot_longer(
    cols = starts_with("SDI_"), 
    names_to = "Trait", 
    values_to = "SDI"
  ) %>%
  mutate(Trait = case_when(
    Trait == "SDI_Body_length" ~ "Body Length",
    Trait == "SDI_Gnathopod_ratio" ~ "Gnathopod Ratio"
  )) %>%
  filter(!is.na(Trait))

# Ensure trait factor levels are ordered
sdi_long$Trait <- factor(sdi_long$Trait, levels = c("Body Length", "Gnathopod Ratio"))

# Function to italicize species labels
italicize_species <- function(species) {
  paste0("italic('", species, "')")
}

# Define custom colors
trait_colors <- c(
  "Body Length" = "#4E79A7",
  "Gnathopod Ratio" = "#76B7B2"
)

# Define trait labels
trait_labels <- c(
  "Body Length" = "Body Length",
  "Gnathopod Ratio" = "Gnathopod 2 Length / Body Length"
)

# Filter out irrelevant combinations and clean SDI values
sdi_filtered <- sdi_long %>%
  filter(
    !(Species == "A. ramondi" & Habitat == "Black coral forests"),
    !(Species == "C. acanthifera" & Habitat %in% c("Rhodolith beds", "Seagrass meadows")),
    !is.na(SDI),
    SDI != 0
  )

# Drop unused factor levels
sdi_filtered$Habitat <- droplevels(sdi_filtered$Habitat)

# ---- SDI Plot by Species and Habitat ----

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

# Clean and recode habitat labels
sdi_results <- sdi_results %>%
  mutate(Habitat = factor(Habitat, levels = c("R", "M", "S", "B"))
  )

# Pivot longer for plotting
sdi_long <- sdi_results %>%
  pivot_longer(
    cols = starts_with("SDI_"), 
    names_to = "Trait", 
    values_to = "SDI"
  ) %>%
  mutate(Trait = case_when(
    Trait == "SDI_Body_length" ~ "Body Length",
    Trait == "SDI_Gnathopod_ratio" ~ "Gnathopod Ratio",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Trait))

# Ensure trait factor levels are ordered
sdi_long$Trait <- factor(sdi_long$Trait, levels = c("Body Length", "Gnathopod Ratio"))

# Function to italicize species labels
italicize_species <- function(species) {
  paste0("italic('", species, "')")
}

# Define trait labels for legend
trait_labels <- c(
  "Body Length" = "Body Length",
  "Gnathopod Ratio" = "Gnathopod Ratio"
)

# Filter out irrelevant combinations and clean SDI values
sdi_filtered <- sdi_long %>%
  filter(
    !(Species == "A. ramondi" & Habitat == "Black coral forests"),
    !(Species == "C. acanthifera" & Habitat %in% c("Rhodolith beds", "Seagrass meadows")),
    !is.na(SDI),
    SDI != 0
  )

# Drop unused factor levels
sdi_filtered$Habitat <- droplevels(sdi_filtered$Habitat)

# ---- Plot SDI ----
# Only two colors: turquoise and orange
trait_colors <- c("Body Length" = "#40E0D0",  # Turquoise shade
                  "Gnathopod Ratio" = "#FFA500") # Orange shade

plot_sdi <- ggplot(sdi_filtered, aes(x = Habitat, y = SDI, fill = Trait)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black") +
  facet_wrap(~ Species, scales = "free_x", labeller = labeller(
    Species = as_labeller(setNames(
      sapply(unique(sdi_filtered$Species), function(sp) parse(text = italicize_species(sp))),
      unique(sdi_filtered$Species)
    ))
  )) +
  scale_fill_manual(values = trait_colors, labels = trait_labels, name = "Traits") +
  labs(
    x = NULL,
    y = "Sexual Dimorphism Index (SDI)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(face = "italic", size = 20),
    axis.title.x = element_text(face = "bold", size = 25),
    axis.title.y = element_text(face = "bold", size = 30),   
    axis.text.x = element_text(size = 30, face = "bold"),   # Bigger X-axis labels
    axis.text.y = element_text(size = 22, face = "bold"),   
    legend.title = element_text(face = "bold", size = 26),  # Bigger legend title
    legend.text = element_text(face = "bold", size = 24),   # Bigger legend text
    panel.grid.major.x = element_blank()
  )

# Show the plot
plot_sdi

# Save the figure
ggsave(
  filename = "SDI_by_species_and_habitat.png",
  plot = plot_sdi,
  width = 18,
  height = 10,
  dpi = 600,
  units = "in",
  bg = "white"
)

# ---- Sex-specific Allometry Analysis ----

# Filter for the two species
combined_species <- data_mod %>% 
  filter(Species %in% c("Ampithoe ramondi", "Caprella acanthifera"))

# Create combined plot
plot_sex_allometry_combined <- ggplot(combined_species, 
                                      aes(x = log(Body_length), y = log(Gnathopod_2_length), color = Sex_label)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.2) +
  facet_wrap(~ Species, labeller = labeller(
    Species = as_labeller(setNames(
      lapply(unique(combined_species$Species), function(sp) bquote(italic(.(sp)))),
      unique(combined_species$Species)
    ))
  )) +
  scale_color_manual(values = c("Female" = "#e377c2", "Male" = "#56B4E9"), name = "Sex") +
  labs(
    x = "log(Body length)",
    y = "log(Gnathopod 2 length)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(face = "italic", size = 20),
    axis.title = element_text(face = "bold", size = 20),
    axis.text = element_text(face = "bold", size = 16),
    legend.title = element_text(face = "bold", size = 18),
    legend.text = element_text(face = "bold", size = 16),
    panel.grid = element_blank(),          # remove grid
    axis.line = element_line(size = 1.2),  # draw axes
    panel.border = element_blank()         # no border around panel
  )

# Display the plot
print(plot_sex_allometry_combined)

library(broom)

# Get stats for every line in the plot
allometry_stats <- combined_species %>%
  group_by(Species, Sex_label) %>%
  do({
    mod <- lm(log(Gnathopod_2_length) ~ log(Body_length), data = .)
    stats <- tidy(mod, conf.int = TRUE)
    summ <- summary(mod)
    data.frame(
      term = stats$term,
      estimate = stats$estimate,
      conf.low = stats$conf.low,
      conf.high = stats$conf.high,
      p.value = stats$p.value,
      r = sqrt(summ$r.squared),
      n = nrow(.)
    )
  })

print(allometry_stats)

# Save the figure
ggsave("amphipod_paper_figure_1.png", plot = plot_sex_allometry_combined, 
       width = 18, height = 10, units = "in", dpi = 600, bg = "white")

# For supplementary information----

# Combined allometry plot for both species
plot_sex_allometry_combined <- function(data, species_list = c("Ampithoe ramondi", "Caprella acanthifera")) {
  
  # Filter only the species of interest
  species_data <- data %>% filter(Species %in% species_list)
  
  # Fit linear models separately for each species and print summaries
  for (sp in species_list) {
    model <- lm(log(Gnathopod_2_length) ~ log(Body_length) * Sex_label, data = species_data %>% filter(Species == sp))
    cat("\nModel summary for", sp, ":\n")
    print(summary(model))
    
    # Residual diagnostics
    sim_res <- simulateResiduals(model)
    plot(sim_res)
  }
  
  # Create combined allometry plot
  plot <- ggplot(species_data, aes(x = log(Body_length), y = log(Gnathopod_2_length), color = Sex_label)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 1.2) +
    facet_wrap(~ Species, labeller = labeller(Species = as_labeller(setNames(
      lapply(species_list, function(sp) bquote(italic(.(sp)))),
      species_list
    )))) +
    labs(
      x = "log(Body length)",
      y = "log(Gnathopod 2 length)",
      color = "Sex"
    ) +
    scale_color_manual(values = c("Female" = "#e377c2", "Male" = "#56B4E9")) +
    theme_minimal(base_size = 14) +
    theme(
      axis.title = element_text(face = "bold", size = 20),
      axis.text = element_text(face = "bold", size = 16),
      legend.title = element_text(face = "bold", size = 18),
      legend.text = element_text(face = "bold", size = 16),
      strip.text = element_text(face = "italic", size = 18)
    )
  
  print(plot)
  
  # Save combined plot
  ggsave("allometry_combined.png", plot = plot, width = 18, height = 10, units = "in", dpi = 600, bg = "white")
}

# Run combined plot
plot_sex_allometry_combined(data_mod)

# ---- Sex-specific Allometry Plot ----

# Filter data for target species and log-transform traits

species_for_manova <- c("Ampithoe ramondi", "Caprella acanthifera", 
                        "Gammaropsis ostromowi", "Pseudoprotella phasma", 
                        "Aora spinicornis", "Erichtonius punctatus")
data_manova <- subset(data_mod, Species %in% species_for_manova)

data_allometry <- data_manova %>%
  filter(Species %in% c("Ampithoe ramondi", "Caprella acanthifera")) %>%
  mutate(
    log_body = log(Body_length),
    log_gnath = log(Gnathopod_2_length)
  )

# Markdown-formatted species labels for italics in facets
species_labels <- c(
  "Ampithoe ramondi" = "*Ampithoe ramondi*",
  "Caprella acanthifera" = "*Caprella acanthifera*"
)

# Create allometry plot
plot_allometry <- ggplot(data_allometry, aes(x = log_body, y = log_gnath, color = Sex_label)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1.2) +
  facet_wrap(~ Species, labeller = labeller(Species = species_labels)) +
  scale_color_manual(
    values = c("Female" = "#e377c2", "Male" = "#56B4E9"),
    name = "Sex"
  ) +
  labs(
    title = "Sex-specific Allometry of Gnathopod 2 Length",
    x = "log(Body Length)",
    y = "log(Gnathopod 2 Length)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    strip.text = ggtext::element_markdown(size = 14)  # Requires ggtext package
  )
print(plot_allometry)

# Save the plot
ggsave(
  filename = "Sex_specific_allometry.png",
  plot = plot_allometry,
  width = 18,
  height = 10,
  dpi = 600,
  units = "in",
  bg = "white"
)

# ---- Pooled Trait Distributions by Habitat and Sex ----

# Convert wide to long format
df_long <- pivot_longer(
  data_manova,
  cols = c(Body_length, Gnathopod_2_length, Gnathopod_ratio),
  names_to = "Trait",
  values_to = "Value"
)

# Set descriptive habitat labels
df_long$Habitat <- factor(df_long$Habitat,
                          levels = c("R", "M", "S", "B"),
                          labels = c("Rhodolith beds", "Macroalgae", "Seagrass meadows", "Black coral forests"))

# Define color palette for sexes
sex_colors <- c("Female" = "#e377c2", "Male" = "#1f77b4")

# Plot pooled traits
plot_traits <- ggplot(df_long, aes(x = Habitat, y = Value, fill = Sex_label)) +
  geom_boxplot(alpha = 0.6, position = position_dodge(width = 0.8), outlier.shape = NA) +
  geom_jitter(aes(color = Sex_label),
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
              alpha = 0.5, size = 1) +
  facet_wrap(~ Trait, scales = "free_y",
             labeller = as_labeller(c(
               "Body_length" = "Body length",
               "Gnathopod_2_length" = "Gnathopod 2 length",
               "Gnathopod_ratio" = "Gnathopod 2 / Body length"
             ))) +
  scale_fill_manual(values = sex_colors, name = "Sex") +
  scale_color_manual(values = sex_colors, name = "Sex") +
  labs(
    x = "Habitats",
    y = "Trait values",
   # title = "Pooled Trait Distributions by Habitat and Sex"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18)
  )

# Display plot
print(plot_traits)

# Save high-resolution plot
ggsave(
  filename = "Trait_distributions_by_habitat_and_sex.png",
  plot = plot_traits,
  width = 18,
  height = 10,
  dpi = 600,
  units = "in",
  bg = "white"
)

# MANOVA----

manova_model <- manova(cbind(Body_length, Gnathopod_2_length, Gnathopod_ratio) ~ 
                         Sex_label * Habitat + Species, data = data_manova)
summary(manova_model, test = "Pillai")

lmer_body <- lmer(Body_length ~ Sex_label * Habitat + (1|Species), data = data_manova)
summary(lmer_body)
unique(data_manova$Species)
# Diagnostic plots
par(mfrow = c(1, 2))
plot(fitted(lmer_body), resid(lmer_body), main = "Residuals vs Fitted")
abline(h = 0, col = "red")
qqnorm(resid(model_ratio), main = "QQ Plot of Residuals")
qqline(resid(model_ratio), col = "red")

# Add small constant to avoid log(0)
data_manova$log_Gnathopod2 <- log(data_manova$Gnathopod_2_length + 0.001)

# Run models on transformed variables
lmer_log_gna <- lmer(log_Gnathopod2 ~ Sex_label * Habitat + (1 | Species), data = data_manova)
summary(lmer_log_gna)

# Plot residuals again
par(mfrow = c(1, 2))
plot(fitted(lmer_log_gna), resid(lmer_log_gna), main="Residuals vs Fitted: log Gnathopod 2 length", xlab="Fitted", ylab="Residuals")
abline(h=0, col="red", lty=2)

qqnorm(resid(lmer_log_gna)); qqline(resid(lmer_log_gna), col="red")

# Model gnathopod ratio ~ Sex * Habitat + random species effect
model_ratio <- lmer(Gnathopod_ratio ~ Sex_label * Habitat + (1 | Species), data = data_manova)

# Summary to check results
summary(model_ratio)

# Diagnostic plots
par(mfrow = c(1, 2))
plot(fitted(model_ratio), resid(model_ratio), main = "Residuals vs Fitted")
abline(h = 0, col = "red")
qqnorm(resid(model_ratio), main = "QQ Plot of Residuals")
qqline(resid(model_ratio), col = "red")

# ---- Univariate models for Ampithoe ramondi and Caprella acanthifera ----

# Log transform for A. ramondi
data_ramondi$log_Gnathopod2 <- log(data_ramondi$Gnathopod_2_length + 0.001)

# ---- Ampithoe ramondi ----

# Body length model
lm_body_ramondi <- lm(Body_length ~ Sex_label * Habitat, data = data_ramondi)
summary(lm_body_ramondi)

# Diagnostic plots
par(mfrow = c(1, 2))
plot(fitted(lm_body_ramondi), resid(lm_body_ramondi), 
     main = "Residuals vs Fitted: Body length", 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red", lty = 2)
qqnorm(resid(lm_body_ramondi)); qqline(resid(lm_body_ramondi), col = "red")

library(emmeans)
library(multcomp)

# This will show you exactly which pairs are significantly different
pairs(em_interaction, adjust = "tukey")

# # Log Gnathopod 2 model
# lm_log_gna_ramondi <- lm(log_Gnathopod2 ~ Sex_label * Habitat, data = data_ramondi)
# summary(lm_log_gna_ramondi)
# 
# # Diagnostics
# par(mfrow = c(1, 2))
# plot(fitted(lm_log_gna_ramondi), resid(lm_log_gna_ramondi),
#      main = "Residuals vs Fitted: log(Gnathopod 2 length) (A. ramondi)",
#      xlab = "Fitted values", ylab = "Residuals")
# abline(h = 0, col = "red", lty = 2)
# qqnorm(resid(lm_log_gna_ramondi)); qqline(resid(lm_log_gna_ramondi), col = "red")

# Beta regression: Gnathopod ratio
beta_ratio_ramondi <- glmmTMB(Gnathopod_ratio ~ Sex_label * Habitat,
                              family = beta_family(link = "logit"),
                              data = data_ramondi)
summary(beta_ratio_ramondi)

# Get the pairwise comparisons for the interaction
pairs(emmeans(beta_ratio_ramondi, ~ Sex_label * Habitat, type = "response"), adjust = "tukey")

# DHARMa residual diagnostics
sim_resid_ramondi <- simulateResiduals(beta_ratio_ramondi)
plot(sim_resid_ramondi)

# ---- ggplots for A. ramondi ----

# Shared color palette
sex_colors <- c("Female" = "#e377c2", "Male" = "#1f77b4")

# Body length plot for Ampithoe ramondi with ultra-big bold text
plot_body_ramondi <- ggplot(data_ramondi, aes(x = Habitat, y = Body_length, fill = Sex_label)) +
  geom_boxplot(position = position_dodge(0.8), alpha = 0.6, outlier.shape = NA) +
  geom_jitter(aes(color = Sex_label), position = position_jitterdodge(0.15, 0.8), alpha = 0.6, size = 3) +
  scale_fill_manual(values = sex_colors) +
  scale_color_manual(values = sex_colors) +
  scale_x_discrete(labels = c("R" = "Rhodolith beds", "M" = "Macroalgae-dominated reefs", "S" = "Seagrass meadows")) +
  labs(
    x = NULL,
    y = "Body length (mm)",
    fill = "Sex",
    color = "Sex"
  ) +
  theme_minimal(base_size = 24) +  # base font really big
  theme(
    axis.title = element_text(face = "bold", size = 40),   # massive axis titles
    axis.text = element_text(face = "bold", size = 32),    # huge tick labels
    axis.text.x = element_text(face = "bold", size = 25),  # make x-axis labels bold
    legend.title = element_text(face = "bold", size = 36), # huge legend title
    legend.text = element_text(face = "bold", size = 32),  # huge legend text
    plot.title = element_text(face = "bold", hjust = 0.5, size = 40),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 1.5)  # thicker axes
  )

# Print plot
print(plot_body_ramondi)

# Save high-resolution, publication-ready figure
ggsave(
  filename = "amphipod_paper_figure_2.png",
  plot = plot_body_ramondi,
  width = 18,    # increase width for clarity
  height = 12,   # increase height for clarity
  units = "in",
  dpi = 600,     # high resolution for journals
  bg = "white"
)


# Gnathopod 2 length plot
plot_gnath2_ramondi <- ggplot(data_ramondi, aes(x = Habitat, y = Gnathopod_2_length, fill = Sex_label)) +
  geom_boxplot(position = position_dodge(0.8), alpha = 0.6, outlier.shape = NA) +
  geom_jitter(aes(color = Sex_label), position = position_jitterdodge(0.15, 0.8), alpha = 0.5, size = 1) +
  scale_fill_manual(values = sex_colors) + scale_color_manual(values = sex_colors) +
  scale_x_discrete(labels = c("R" = "Rhodolith beds", "M" = "Macroalgae", "S" = "Seagrass meadows")) +
  labs(x = "Habitats", y = "Gnathopod 2 length (mm)", title = expression("Gnathopod 2 length in "*italic("A. ramondi")*" by Habitat and Sex")) +
  theme_minimal(base_size = 16) +
  theme(legend.position = "top", legend.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"), plot.title = element_text(hjust = 0.5, face = "bold"))

# Gnathopod ratio plot for Ampithoe ramondi with ultra-big bold text
plot_ratio_ramondi <- ggplot(data_ramondi, aes(x = Habitat, y = Gnathopod_ratio, fill = Sex_label)) +
  geom_boxplot(position = position_dodge(0.8), alpha = 0.6, outlier.shape = NA) +
  geom_jitter(aes(color = Sex_label), position = position_jitterdodge(0.15, 0.8), alpha = 0.6, size = 3) +
  scale_fill_manual(values = sex_colors) +
  scale_color_manual(values = sex_colors) +
  scale_x_discrete(labels = c("R" = "Rhodolith beds", "M" = "Macroalgae-dominated reefs", "S" = "Seagrass meadows")) +
  labs(
    x = NULL,
    y = "Gnathopod 2 length / Body length",
    fill = "Sex",
    color = "Sex"
  ) +
  theme_minimal(base_size = 24) +  # base font very large
  theme(
    axis.title = element_text(face = "bold", size = 40),   # massive axis titles
    axis.text = element_text(face = "bold", size = 32),    # tick labels
    axis.text.x = element_text(face = "bold", size = 24),  # x-axis labels
    legend.title = element_text(face = "bold", size = 36), # legend title
    legend.text = element_text(face = "bold", size = 32),  # legend text
    plot.title = element_text(face = "bold", hjust = 0.5, size = 40),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 1.5)  # thicker axes
  )

# Print the plot
print(plot_ratio_ramondi)

# Save high-resolution, publication-ready figure
ggsave(
  filename = "amphipod_paper_figure_3_ultrabig.png",
  plot = plot_ratio_ramondi,
  width = 18,   # width for clarity
  height = 12,  # height for clarity
  units = "in",
  dpi = 600,
  bg = "white"
)

# # ---- Save A. ramondi plots ----
# ggsave("Body_length_A_ramondi.png", plot_body_ramondi, width = 18, height = 10, units = "in", dpi = 600, bg = "white")
# ggsave("Gnathopod2_A_ramondi.png", plot_gnath2_ramondi, width = 18, height = 10, units = "in", dpi = 600, bg = "white")
# ggsave("Gnathopod_ratio_A_ramondi.png", plot_ratio_ramondi, width = 18, height = 10, units = "in", dpi = 600, bg = "white")

# ---- Repeat similarly for C. acanthifera ---

# ---- UNIVARIATE MODELS FOR Caprella acanthifera ----

# ------------------- Body length model 

lm_body_acanthifera <- lm(Body_length ~ Sex_label * Habitat, data = data_acanthifera)
summary(lm_body_acanthifera)
# Get the pairwise comparisons for the interaction
pairs(emmeans(lm_body_acanthifera, ~ Sex_label * Habitat), adjust = "tukey")

# Diagnostic plots for Body length
par(mfrow = c(1, 2))
plot(fitted(lm_body_acanthifera), resid(lm_body_acanthifera),
     main = "Residuals vs Fitted: Body length",
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red", lty = 2)
qqnorm(resid(lm_body_acanthifera))
qqline(resid(lm_body_acanthifera), col = "red")

# # ------------------- Log-transformed Gnathopod 2 length 
# 
# # Small offset to avoid log(0)
# data_acanthifera$log_Gnathopod2 <- log(data_acanthifera$Gnathopod_2_length + 0.001)
# 
# lm_log_gna_acanthifera <- lm(log_Gnathopod2 ~ Sex_label * Habitat, data = data_acanthifera)
# summary(lm_log_gna_acanthifera)
# 
# # Diagnostic plots for log-transformed Gnathopod 2 length
# par(mfrow = c(1, 2))
# plot(fitted(lm_log_gna_acanthifera), resid(lm_log_gna_acanthifera),
#      main = "Residuals vs Fitted: log Gnathopod 2 length (C. acanthifera)",
#      xlab = "Fitted values", ylab = "Residuals")
# abline(h = 0, col = "red", lty = 2)
# qqnorm(resid(lm_log_gna_acanthifera))
# qqline(resid(lm_log_gna_acanthifera), col = "red")

# ------------------- Beta regression for Gnathopod ratio 

beta_ratio_acanthifera <- glmmTMB(Gnathopod_ratio ~ Sex_label * Habitat,
                                  family = beta_family(link = "logit"),
                                  data = data_acanthifera)
summary(beta_ratio_acanthifera)
# Pairwise comparisons for Gnathopod ratio in C. acanthifera
pairs(emmeans(beta_ratio_acanthifera, ~ Sex_label * Habitat, type = "response"), adjust = "tukey")

# Residual diagnostics for beta regression
sim_resid_acanthifera <- simulateResiduals(beta_ratio_acanthifera)
plot(sim_resid_acanthifera)

# ------------------- ggplot2 Visualization 

sex_colors <- c("Female" = "#e377c2", "Male" = "#1f77b4")

# Body length plot for Caprella acanthifera - big bold
plot_body_acanthifera <- ggplot(data_acanthifera, aes(x = Habitat, y = Body_length, fill = Sex_label)) +
  geom_boxplot(alpha = 0.6, position = position_dodge(0.8), outlier.shape = NA, width = 0.6) +
  geom_jitter(aes(color = Sex_label),
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
              alpha = 0.6, size = 3) +
  scale_fill_manual(values = sex_colors, name = "Sex") +
  scale_color_manual(values = sex_colors, name = "Sex") +
  scale_x_discrete(labels = c("R" = "Rhodolith beds", 
                              "M" = "Macroalgae-dominated reefs", 
                              "S" = "Seagrass meadows", 
                              "B" = "Black coral forests")) +
  labs(
    x = NULL,
    y = "Body length (mm)"
  ) +
  theme_minimal(base_size = 24) +
  theme(
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 36),
    legend.text  = element_text(face = "bold", size = 32),
    axis.title   = element_text(face = "bold", size = 40),
    axis.text    = element_text(face = "bold", size = 32),
    axis.text.x  = element_text(face = "bold", size = 32),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 1.5)
  )

# Gnathopod ratio plot for Caprella acanthifera - big bold
plot_ratio_acanthifera <- ggplot(data_acanthifera, aes(x = Habitat, y = Gnathopod_ratio, fill = Sex_label)) +
  geom_boxplot(alpha = 0.6, position = position_dodge(0.8), outlier.shape = NA, width = 0.6) +
  geom_jitter(aes(color = Sex_label),
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
              alpha = 0.6, size = 3) +
  scale_fill_manual(values = sex_colors, name = "Sex") +
  scale_color_manual(values = sex_colors, name = "Sex") +
  scale_x_discrete(labels = c("R" = "Rhodolith beds", 
                              "M" = "Macroalgae-dominated reefs", 
                              "S" = "Seagrass meadows", 
                              "B" = "Black coral forests")) +
  labs(
    x = NULL,
    y = "Gnathopod 2 length / Body length"
  ) +
  theme_minimal(base_size = 24) +
  theme(
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 36),
    legend.text  = element_text(face = "bold", size = 32),
    axis.title   = element_text(face = "bold", size = 40),
    axis.text    = element_text(face = "bold", size = 32),
    axis.text.x  = element_text(face = "bold", size = 32),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 1.5)
  )

# Print plots
print(plot_body_acanthifera)
print(plot_ratio_acanthifera)

# Save high-resolution publication-ready figures
ggsave("amphipod_paper_figure_4.png", plot = plot_body_acanthifera,
       width = 18, height = 12, units = "in", dpi = 600, bg = "white")

ggsave("amphipod_paper_figure_5.png", plot = plot_ratio_acanthifera,
       width = 18, height = 12, units = "in", dpi = 600, bg = "white")

# # Save plots
# 
# ggsave("Body_length_C_acanthifera.png", plot_body_acanthifera, width = 18, height = 10, units = "in", dpi = 600, bg = "white")
# ggsave("Gnathopod2_C_acanthifera.png", plot_gnath2_acanthifera, width = 18, height = 10, units = "in", dpi = 600, bg = "white")
# ggsave("Gnathopod_ratio_C_acanthifera.png", plot_ratio_acanthifera, width = 18, height = 10, units = "in", dpi = 600, bg = "white")

# ---- Estimated Marginal Means Interaction Plots (Sex × Habitat) ----

# Define habitat order and labels
habitat_labels <- c("Rhodolith beds", "Macroalgae", "Seagrass meadows", "Black coral forests")
habitat_codes <- c("R", "M", "S", "B")

# Estimated Marginal Means: Body Length ---------------------------------------

# Get EMMs from model
emm_body <- emmeans(lmer_body, ~ Sex_label * Habitat)
emm_df_body <- as.data.frame(emm_body)

# Recode and order habitat
emm_df_body$Habitat <- factor(emm_df_body$Habitat, 
                              levels = habitat_codes,
                              labels = habitat_labels)

# Plot
# Interaction plot for Body Length with improved styling and saving
plot_interaction_body <- ggplot(emm_df_body, aes(x = Habitat, y = emmean, color = Sex_label, group = Sex_label)) +
  geom_point(position = position_dodge(0.4), size = 3) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), 
                position = position_dodge(0.4), width = 0.2) +
  geom_line(position = position_dodge(0.4), linewidth = 1.2) +
  scale_color_manual(values = c("Female" = "#e377c2", "Male" = "#56B4E9")) +
  labs(
    x = NULL,
    y = "Estimated Marginal Mean",
    color = "Sex"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold", size = 20),
    axis.text = element_text(face = "bold", size = 16),
    legend.title = element_text(face = "bold", size = 18),
    legend.text = element_text(face = "bold", size = 16)
  )

# Print the plot
print(plot_interaction_body)

# Save the plot as PNG
ggsave("interaction_body_length.png", plot = plot_interaction_body, width = 14, height = 8, units = "in", dpi = 600, bg = "white")

# Estimated Marginal Means: Gnathopod 2 Length -------------------------------

emm_gnath <- emmeans(lmer_log_gna, ~ Sex_label * Habitat)
emm_df_gnath <- as.data.frame(emm_gnath)

emm_df_gnath$Habitat <- factor(emm_df_gnath$Habitat, 
                               levels = habitat_codes,
                               labels = habitat_labels)

# Interaction plot for Gnathopod 2 Length with consistent styling
plot_interaction_gnath <- ggplot(emm_df_gnath, aes(x = Habitat, y = emmean, color = Sex_label, group = Sex_label)) +
  geom_point(position = position_dodge(0.4), size = 3) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), 
                position = position_dodge(0.4), width = 0.2) +
  geom_line(position = position_dodge(0.4), linewidth = 1.2) +
  scale_color_manual(values = c("Female" = "#e377c2", "Male" = "#56B4E9")) +
  labs(
    x = NULL,
    y = "Estimated Marginal Mean (log-transformed)",
    color = "Sex"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold", size = 20),
    axis.text = element_text(face = "bold", size = 16),
    legend.title = element_text(face = "bold", size = 18),
    legend.text = element_text(face = "bold", size = 16)
  )

# Print the plot
print(plot_interaction_gnath)

# Save the plot
ggsave("interaction_gnathopod_length.png", plot = plot_interaction_gnath, width = 14, height = 8, units = "in", dpi = 600, bg = "white")

# Estimated Marginal Means: Gnathopod Ratio ----------------------------------

emm_ratio <- emmeans(model_ratio, ~ Sex_label * Habitat)
emm_df_ratio <- as.data.frame(emm_ratio)

emm_df_ratio$Habitat <- factor(emm_df_ratio$Habitat, 
                               levels = habitat_codes,
                               labels = habitat_labels)

# Interaction plot for Gnathopod 2 / Body Length Ratio
plot_interaction_ratio <- ggplot(emm_df_ratio, aes(x = Habitat, y = emmean, color = Sex_label, group = Sex_label)) +
  geom_point(position = position_dodge(0.4), size = 3) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), 
                position = position_dodge(0.4), width = 0.2) +
  geom_line(position = position_dodge(0.4), linewidth = 1.2) +
  scale_color_manual(values = c("Female" = "#e377c2", "Male" = "#56B4E9")) +
  labs(
    x = NULL,
    y = "Estimated Marginal Mean",
    color = "Sex"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold", size = 20),
    axis.text = element_text(face = "bold", size = 16),
    legend.title = element_text(face = "bold", size = 18),
    legend.text = element_text(face = "bold", size = 16)
  )

# Print the plot
print(plot_interaction_ratio)

# Save the plot
ggsave("interaction_gnathopod_ratio.png", plot = plot_interaction_ratio, width = 14, height = 8, units = "in", dpi = 600, bg = "white")

# ---- Save Estimated Mean Plots ----

ggsave("interaction_body_length.png", plot = plot_interaction_body, width = 18, height = 10, units = "in", dpi = 600, bg = "white")
ggsave("interaction_gnathopod_length.png", plot = plot_interaction_gnath, width = 18, height = 10, units = "in", dpi = 600, bg = "white")
ggsave("interaction_gnathopod_ratio.png", plot = plot_interaction_ratio, width = 18, height = 10, units = "in", dpi = 600, bg = "white")
ggsave("Gnathopod2_C_acanthifera.png", plot_gnath2_acanthifera, width = 18, height = 10, units = "in", dpi = 600, bg = "white")

