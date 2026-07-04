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
create_founders <- function(determ=TRUE, n_ind=100, n_chr=10, n_gpchr=10, sp="WHEAT", h2=0.5){
  # Add determinism
  if(determ) {set.seed(923)}
  
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
run_breed_cycle <- function(pop, SP, selected_count, n_cross=100){
  # Add environmental noise to the plants (phenotyping)
  pop <- setPheno(pop, simParam = SP)
  
  # Select the top x based on their phenotype
  parents <- selectInd(pop, nInd=selected_count, use="pheno", simParam=SP)
  
  # Cross the parents to make x new offspring for the next generation
  pop <- randCross(parents, nCrosses=n_cross, simParam=SP)
  
  return(pop)
}


## Repeat breeding cycles
repeat_breed_cycles <- function(pop, SP, selected_count, n_gen, n_cross=100){
  df <- tibble(
    sel_count = selected_count,
    gen = 1:n_gen,
    pop_object = vector("list", n_gen)
  )
  
  for(x in 1:n_gen){
    pop <- run_breed_cycle(pop, SP, selected_count, n_cross)
    df$pop_object[[x]] <- pop
  }
  
  return(df)
}


## Extract information
extract_breeding_info <- function(df, include_advanced=FALSE, keep_objects=FALSE) {
  df1 <- df %>%
    #calculate selection info
    mutate(
      gen_size = map_int(pop_object, nInd),
      sel_frac = sel_count/gen_size,
      sel_int = selInt(sel_frac),
      .before = gen
    ) %>%
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
             sel_response = mean_gv - meanG(pop_sp$pop),
             sel_diff = meanP(pop_sp$pop) - map_dbl(pop_object, meanP),
             h2_realized = sel_response/sel_diff
      ) else .} %>%
    # retain population objects
    {if(!keep_objects) select(., !pop_object) else .}
  
  return(df1)
}


## Wrapper function for multiple selected counts/scenarios
run_mult_selections <- function(population, start_params, sel_counts, scenarios){
  df <- purrr::map2_df(
    .x=sel_counts,.y=scenarios, 
    .f=function(x, y){
      repeat_breed_cycles(pop=population,
                          SP=start_params,
                          selected_count=x,
                          n_gen=10,
                          n_cross=100) %>%
        extract_breeding_info(include_advanced=TRUE, keep_objects=TRUE) %>%
        add_column(scenario = y, .before="sel_count") %>%
        mutate(scenario_short = str_remove(scenario, "(?<= Int).*$"), .after=scenario)
    })
  
  return(df)
  
}



# Plotting Functions================================================================================
## Genetic gain curve
#extract slopes from model
get_slopes <- function(df, y, x="gen", cat="scenario_short", round=TRUE){
  #build formula
  form <- as.formula(
    paste(y,
          paste(x, cat, sep=" * "),
          sep=" ~ ")
  )
  
  mod <- lm(form, data=df)
  
  #extract vector of slopes
  vec_slopes <- emtrends(mod, specs="scenario_short", var="gen") %>%
    as_tibble() %>%
    mutate(gen.trend=round(gen.trend, 3)) %>%
    pull(gen.trend, name=scenario_short) %>%
    paste(names(.), ., sep=": ")
  
  return(vec_slopes)
}

#genetic value increases with generations
plot_gen_gain <- function(df, y, x, cat, n_gen){
  #extract slopes
  vec_m <- get_slopes(df, y)
  
  annotation_text <- paste0("Slopes:\n", 
                            paste(vec_m, collapse = "\n"))
  
  #create plot
  p <- df %>%
    ggplot(aes(x=gen, y=mean_gv, color=scenario)) +
    geom_point() +
    geom_smooth(aes(fill=scenario), alpha=0.2, method="lm") +
    scale_x_continuous(breaks=seq(0, n_gen, 2), 
                       labels=seq(0, n_gen, 2), 
                       limits=c(0, NA),
                       expand=expansion(mult=c(0, 0.05))) +
    scale_y_continuous(limits=c(0, NA),
                       expand=expansion(mult=c(0, 0.05))) +
    scale_color_manual(values=scenario_color_fill_map) +
    scale_fill_manual(values=scenario_color_fill_map) +
    labs(x="Generation",
         y="Mean genetic value") +
    annotate(geom = "label",
             x = 0.5, 
             y = 1.1*max(df[["mean_gv"]]),
             hjust=0,
             vjust=1,
             label=annotation_text,
             size=3.25,
             fill = "yellow",
             color = "black") +
    theme_bw() +
    theme(
      legend.position="bottom",
      legend.title=element_blank()
    )
  
  return(p)
}


## Breeder's dilemma
plot_breeders_dilemma <- function(df, n_gen){
  p <- df %>%
    ggplot() +
    geom_line(aes(x=gen, y=mean_gv),
              color="darkred") +
    geom_point(aes(x=gen, y=mean_gv),
               color="darkred") +
    geom_area(aes(x=gen, y=var_g),
              fill="navy") +
    facet_wrap(~scenario, ncol=1) +
    scale_x_continuous(breaks=seq(0, n_gen, 2),
                       labels=seq(0, n_gen, 2)) +
    scale_y_continuous(
      "Mean genetic value",
      sec.axis = sec_axis(~ .x, name = "Genetic variance")
    ) +
    xlab("Generation") +
    theme_bw()
  
  return(p)
}


## Signal vs Noise
pivot_for_svn <- function(df){
  df2 <- df %>%
    #pivot data into correct format
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
    )
  
  return(df2)
}


plot_signal_vs_noise <- function(df, n_gen){
  p <- df %>%
    pivot_for_svn() %>%
    #create plot
    ggplot() +
    geom_area(aes(x=gen, y=variance, color=var_part, fill=var_part),
              color="black") +
    facet_wrap(~scenario, ncol=1) +
    scale_fill_manual("Variance \nComponent",
                      values=c("Vg"="green", "Ve"="red")) +
    scale_x_continuous(breaks=seq(0, n_gen, 2),
                       labels=seq(0, n_gen, 2)) +
    labs(x="Generation",
         y="Variance") +
    theme_bw() +
    theme(
      legend.position="bottom"
    )
  
  return(p)
}


plot_signal_vs_noise <- function(df, n_gen){
  p <- df %>%
    pivot_for_svn() %>%
    #create plot
    ggplot() +
    geom_area(aes(x=gen, y=variance, color=var_part, fill=var_part),
              color="black") +
    scale_fill_manual("Variance \nComponent",
                      values=c("Vg"="green", "Ve"="red")) +
    scale_x_continuous(breaks=seq(0, n_gen, 2),
                       labels=seq(0, n_gen, 2)) +
    labs(x="Generation",
         y="Variance") +
    theme_bw() +
    theme(
      legend.position="bottom"
    )
  
  return(p)
}


## Accuracy Decay
plot_acc_decay <- function(df, y, x="gen", cat="scenario_short", n_gen) {
  #extract slopes
  vec_m <- get_slopes(df, y="sel_acc")
  
  annotation_text <- paste0("Slopes:\n", 
                            paste(vec_m, collapse = "\n"))
  
  p <- df %>%
    ggplot(aes(x=gen, y=sel_acc, color=scenario)) +
    geom_point() + 
    geom_smooth(aes(fill=scenario), alpha=0.2, method="lm") +
    scale_x_continuous(breaks=seq(0, n_gen, 2),
                       labels=seq(0, n_gen, 2)) +
    scale_color_manual(values=scenario_color_fill_map) +
    scale_fill_manual(values=scenario_color_fill_map) +
    labs(x="Generation",
         y="Selection accuracy") +
    annotate(geom = "label",
             x = .65*n_gen, 
             y = 1.1*max(df[["sel_acc"]]),
             hjust=0,
             vjust=1,
             label=annotation_text,
             size=3.25,
             fill = "yellow",
             color = "black") +
    theme_bw() +
    theme(
      legend.position="bottom",
      legend.title=element_blank()
    )
  
  return(p)
}














