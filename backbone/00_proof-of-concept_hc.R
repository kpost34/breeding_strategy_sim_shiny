# Proof of Concept: Hard-coded


# Load Packages=====================================================================================
## Packages
library(pacman)
pacman::p_load(
  AlphaSimR, 
  tidyverse, 
  broom
)



# Build Population==================================================================================
## 1. Create a founder population
#simulate 100 individuals with 10 chromosomes each
founders <- runMacs(nInd=100, nChr=10, species="WHEAT")


## 2. Set simulation parameters
SP <- SimParam$new(founders)

#add an additive trait (e.g., yield) controlled by 100 genes per chromosome
SP$addTraitA(nQtlPerChr=100)

#set the initial error variance so heritability (h2) is 0.5
SP$setVarE(h2=0.5)


## 3. Initialize population
#create the first generation of actual plants from the founder DNA
pop <- newPop(founders, simParam=SP)



# Run Simulation & Extract Information==============================================================
## One selected count: 10
pop_1 <- pop

#create and seed empty df
df <- tibble(
  gen = 1:10,
  mean_gv = numeric(10),
  var_g = numeric(10),
  var_p = numeric(10),
  sel_acc = numeric(10)
)

#populate df
for(gen in df$gen){
  #a. add environmental noise to the plants (phenotyping)
  pop_1 <- setPheno(pop_1, simParam = SP)
  
  #b. select the top 10% based on their phenotype
  parents <- selectInd(pop_1, nInd=10, use="pheno", simParam=SP)
  
  #c. cross the parents to make 100 new offspring for the next generation
  pop_1 <- randCross(parents, nCrosses=100, simParam=SP)
  
  #d. record the parameters of interest
  df[gen, "mean_gv"] <- meanG(pop_1)
  df[gen, "var_g"] <- varG(pop_1)
  df[gen, "var_p"] <- varP(pop_1)
  df[gen, "sel_acc"] <- cor(pop_1@gv, pop_1@pheno)
}


## Multiple selected counts (fractions)
#initialize the df
n_gen <- 10

df2_map <- tibble(
  sel_count = numeric(n_gen),
  gen = 1:n_gen,
  mean_gv = numeric(n_gen),
  var_g = numeric(n_gen),
  var_p = numeric(n_gen),
  sel_acc = numeric(n_gen)
)

sel_count <- c(5, 10, 15)

nm_sel_count <- paste0("sel_count_", sel_count)


#run the loop
df2 <- sel_count %>%
  purrr::map_df(function(x){
    pop_2 <- pop
    for(g in 1:n_gen){
      
      #a. add environmental noise to the plants (phenotyping)
      pop_2 <- setPheno(pop_2, simParam = SP)
      
      #b. select the top 10% based on their phenotype
      parents <- selectInd(pop_2, nInd=x, use="pheno", simParam=SP)
      
      #c. cross the parents to make 100 new offspring for the next generation
      pop_2 <- randCross(parents, nCrosses=100, simParam=SP)
      
      #d. record the parameters of interest
      df2_map[g, "sel_count"] <- x
      df2_map[g, "mean_gv"] <- meanG(pop_2)
      df2_map[g, "var_g"] <- varG(pop_2)
      df2_map[g, "var_p"] <- varP(pop_2)
      df2_map[g, "sel_acc"] <- cor(pop_2@gv, pop_2@pheno)
    }
    df2_map
  })



# View Results======================================================================================
## Genetic Gain Curve
### Single selected count
#plot
df %>%
  ggplot(aes(x=gen, y=mean_gv)) +
  geom_point() +
  geom_smooth(method="lm") +
  theme_bw()
#Pattern: genetic value increases with generation

#model
mod1 <- lm(mean_gv ~ gen, data=df)
m_mod1 <- coef(mod1)[2]
m_mod1


## Breeder's Dilemma
df %>%
  ggplot() +
  geom_line(aes(x=gen, y=mean_gv),
            color="darkred") +
  geom_area(aes(x=gen, y=var_g),
            fill="navy") +
  scale_x_continuous(breaks=1:10,
                     labels=1:10) +
  scale_y_continuous(
    "Mean Genetic Value",
    sec.axis = sec_axis(~ .x, name = "Genetic Variance")
  ) +
  theme_bw()
#Pattern: as genetic value (performance) increases (over generations) while genetic variance
#decreases


## Signal vs Noise (Vg vs Vp)
#area chart of Vp = Vg + Ve
df %>%
  select(!c(mean_gv, sel_acc)) %>%
  pivot_longer(
    cols=c(var_g, var_p),
    names_to="var_part",
    values_to="variance"
  ) %>%
  mutate(
    var_part=ifelse(
      var_part=="var_g", 
      "Vg",
      "Ve"
    )
  ) %>%
  ggplot() +
  geom_area(aes(x=gen, y=variance, color=var_part, fill=var_part),
            color="black") +
  scale_fill_manual("Variance \nComponent",
                    values=c("Vg"="green", "Ve"="red")) +
  scale_x_continuous(labels=1:n_gen,
                     breaks=1:n_gen) +
  xlab("Generation") +
  theme_bw() +
  theme(
    legend.position="bottom"
  )
#Pattern: Vp (top line) is noisy as it comprises large Ve


## Accuracy Decay
#selection accuracy = correlation between GV and phenotype over time

### Single selected count
df %>%
  ggplot(aes(x=gen, y=sel_acc)) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE) +
  theme_bw()
#Pattern: selection accuracy declines over time (generations)


### Multiple selected counts
#plot
df2 %>%
  ggplot(aes(x=gen, y=sel_acc, color=as.factor(sel_count))) +
  geom_point() + 
  geom_smooth(method="lm", se=FALSE) +
  labs(x="Generation",
       y="Selection accuracy",
       color="Selected count") +
  scale_color_viridis_d(end=0.6) +
  theme_bw() +
  theme(legend.position="bottom")
#Pattern: slope increases with smaller selected count

#slope
m_sel_strength <- df2 %>%
  group_by(sel_count) %>%
  nest() %>%
  mutate(slope = map_dbl(data, ~coef(lm(mean_gv ~ gen, data = .x))[2])) %>%
  select(sel_count, slope)

m_sel_strength


