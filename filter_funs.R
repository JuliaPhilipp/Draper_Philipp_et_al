# ASTC Filter
filter_ASTC <- function(event_df, dpsi_df, species){
  sample_name_filter <- paste0(str_sub(species,1,1),"ips")
  result <- event_df %>% 
    filter(str_detect(sample_name, sample_name_filter)) %>% 
    left_join(dpsi_df, by = "event_id") %>% 
    filter(constants == species) %>% 
    filter(sign(min_dpsi) != 0) %>% 
    filter(sign(min_dpsi) == sign(max_dpsi)) %>% 
    filter(abs_inner_dpsi > 0.1, min_jc_comparison >= 15)
  
  return(result)
}

# AS filter
filter_AS <- function(event_df, dpsi_df, species){
  sample_name_filter <- paste0(str_sub(species,1,1),"ips")
  result <- event_df %>%
    filter(str_detect(sample_name, sample_name_filter)) %>% 
    left_join(dpsi_df, by = "event_id") %>% 
    filter(constants == species) %>%
    group_by(event_id) %>%
    mutate(as_only = (all(inner_dpsi == 0 & abs(min_dpsi) <= 0.05) & all(abs(max_dpsi) <= 0.05) & any(min_jc_comparison >= 15)) & any(avg_psi < 1 & avg_psi > 0)) %>%
    filter(as_only == TRUE) 
  
  return(result)
}

filter_AS_strict <- function(event_df, dpsi_df, species){
  sample_name_filter <- paste0(str_sub(species,1,1),"ips")
  result <- event_df %>%
    filter(str_detect(sample_name, sample_name_filter)) %>% 
    left_join(dpsi_df, by = "event_id") %>% 
    filter(constants == species) %>%
    group_by(event_id) %>%
    mutate(as_only = (all(inner_dpsi == 0 & abs(min_dpsi) <= 0.05) & all(abs(max_dpsi) <= 0.05) & any(min_jc_comparison >= 15)) & all(avg_psi < 1 & avg_psi > 0)) %>%
    filter(as_only == TRUE)
  
  return(result)
}

# NMD filter and name conversion
filter_always_NMD <- function(event_df, name_conversion_df, species){
  name_conversion_species_df <- name_conversion_df %>% 
    filter(species_1 == species) %>% 
    select(final_event_id, original_event_id_species_1)

  result <- event_df %>% 
    filter(nmd_switch == "never") %>% 
    filter(nmd_form == "both") %>% 
    left_join(name_conversion_species_df, by = c("event_id" = "original_event_id_species_1")) %>% 
    mutate(event_id = final_event_id) %>% 
    select(-final_event_id)
  
  return(result)
}

filter_switch_NMD <- function(event_df, name_conversion_df, species){
  name_conversion_species_df <- name_conversion_df %>% 
    filter(species_1 == species) %>% 
    select(final_event_id, original_event_id_species_1)
  
  result <- event_df %>% 
    filter(nmd_switch %in% c("always","sometimes")) %>% 
    filter(nmd_form %in% c("included","excluded")) %>% 
    left_join(name_conversion_species_df, by = c("event_id" = "original_event_id_species_1")) %>% 
    mutate(event_id = final_event_id) %>% 
    select(-final_event_id)
  
  return(result)
}

filter_never_NMD <- function(event_df, name_conversion_df, species){
  name_conversion_species_df <- name_conversion_df %>% 
    filter(species_1 == species) %>% 
    select(final_event_id, original_event_id_species_1)
  
  result <- event_df %>% 
    filter(nmd_switch == "never") %>% 
    filter(nmd_form == "never") %>% 
    left_join(name_conversion_species_df, by = c("event_id" = "original_event_id_species_1")) %>% 
    mutate(event_id = final_event_id) %>% 
    select(-final_event_id)
  
  return(result)
}
