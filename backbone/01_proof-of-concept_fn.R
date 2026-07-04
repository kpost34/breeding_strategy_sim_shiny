# Proof of Concept: Functionalized


# Load Packages & Functions=========================================================================
## Packages
library(pacman)
pacman::p_load(
  here,
  AlphaSimR, 
  tidyverse, 
  # broom,
  emmeans
)


## Functions
source(here("utils.R"))



# Create Founders, Run Sim, and Extract Values======================================================
set.seed(22)

## Create first population
### 1. Create founding population
pop_sp <- create_founders()


### 2. Run simulation
df <- repeat_breed_cycles(pop=pop_sp$pop,
                          SP=pop_sp$SP,
                          selected_count=5,
                          n_gen=10,
                          n_cross=100)


### 3. Extract values
df_values <- extract_breeding_info(df, include_advanced=TRUE, keep_objects=TRUE)


## Create two populations
sel_counts <- c(5, 15)
scenarios <- c("High Intensity (5)", "Medium Intensity (15)")

df_values2 <- purrr::map2_df(
  .x=sel_counts,.y=scenarios, 
  .f=function(x, y){
    repeat_breed_cycles(pop=pop_sp$pop,
                        SP=pop_sp$SP,
                        selected_count=x,
                        n_gen=10,
                        n_cross=100) %>%
      extract_breeding_info(include_advanced=TRUE, keep_objects=TRUE) %>%
      add_column(scenario = y, .before="sel_count") %>%
      mutate(scenario_short = str_remove(scenario, "(?<= Int).*$"), .after=scenario)
  })

df_values3 <- run_mult_selections(
  population=pop_sp$pop, 
  start_params=pop_sp$SP, 
  sel_counts=c(5, 15),
  scenarios=c("High Intensity (5)", "Medium Intensity (15)")
)



# Generate Plots====================================================================================
scenario_color_fill_map <- c(
  "High Intensity (5)" = "red",
  "Medium Intensity (15)" = "blue"
)


## 1. Genetic Gain Curve
plot_gen_gain(df_values2, y="mean_gv", n_gen=10)


## 2. Breeder's Dilemma
#as genetic value (performance) increases (over generations) while genetic variance
#decreases
plot_breeders_dilemma(df_values2, n_gen=10)



## 3. Signal vs Noise (Vg vs Vp)
#area chart of Vp = Vg + Ve
#Pattern: Vp (top line) is noisy as it comprises large Ve
plot_signal_vs_noise(df_values2, n_gen=10)



## 4. Accuracy Decay
# Patterns: 1) selection accuracy declines over generations; 2) slope of relationship becomes
  #steeper as selection intensity becomes stronger
plot_acc_decay(df_values2, y="sel_acc", n_gen=10)




