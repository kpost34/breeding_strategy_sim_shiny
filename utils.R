# Utility Functions


# Load Packages=====================================================================================
library(pacman)
pacman::p_load(
  AlphaSimR, 
  tidyverse, 
  broom
)



# Backbone Functions================================================================================
## Create founding population
create_founders <- function(n_ind=100, n_chr=10, n_gpchr=100, sp="WHEAT", h2=0.5){
  # Simulate population
  founders <- runMacs(nInd=n_ind, nChr=n_chr, species=sp)
  
  # Set simulation parameters
  sim_params <- SimParam$new(founders)
  
  # Add an additive trait (e.g., yield) controlled by n genes per chromosome
  sim_params$addTraitA(nQtlPerChr=n_gpchr)
  
  # Set the iniital error variance using heritability param
  sim_params$setVarE(h2=h2)
  
  # Create the first generation 
  population <- newPop(founders, simParam=sim_params)
  
  return(
    list(pop=population,
         SP=sim_params
    )
  )
}


## Run one breeding cycle
run_breed_cycle <- function(pop, SP, intensity, n_cross=100){
  # Add environmental noise to the plants (phenotyping)
  pop <- setPheno(pop, simParam = SP)
  
  # Select the top x% based on their phenotype
  parents <- selectInd(pop, nInd=intensity, use="pheno", simParam=SP)
  
  # Cross the parents to make x new offspring for the next generation
  pop <- randCross(parents, nCrosses=n_cross, simParam=SP)
  
  return(pop)
}


## Repeat breeding cycles
repeat_breed_cycles <- function(pop, SP, intensity, n_gen, n_cross=100){
  df <- tibble(
    gen = 1:n_gen,
    pop_object = vector("list", n_gen)
  )
  
  for(x in 1:n_gen){
    pop <- run_breed_cycle(pop, SP, intensity, n_cross)
    df$pop_object[[x]] <- pop
  }
  
  return(df)
}



## Extract information
extract_breeding_info <- function(df, include_advanced=FALSE, keep_objects=FALSE) {
  df1 <- df %>%
    #extract key info
    mutate(
      mean_gv = map_dbl(pop_object, meanG),
      var_g = map_dbl(pop_object, varG),
      var_p = map_dbl(pop_object, varP),
      sel_acc = map2_dbl(map(pop_object, gv), map(pop_object, pheno), cor)
    ) %>%
    #include advanced stats
    {if(include_advanced) 
      mutate(., 
             # parent_mean_p = meanP(pop_sp$pop),
             # gen_mean_p = map_dbl(pop_object, meanP),
             # parent_mean_g = meanG(pop_sp$pop),
             sel_response = mean_gv - meanG(pop_sp$pop),
             sel_diff = meanP(pop_sp$pop) - map_dbl(pop_object, meanP),
             h2_realized = sel_response/sel_diff
      ) else .} %>%
    # retain population objects
    {if(!keep_objects) select(., !pop_object) else .}
  
  return(df1)
}














