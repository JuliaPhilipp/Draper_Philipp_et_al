# RBP functions

# filter for events that are specific types of splicing in each species
# e.g. ASTC in human and AS in chimp, vice versa, or ASTC in both, or AS in both
RBP_by_event_type <- function(df, species1, species2, species1_ASTC, species1_AS, species1_event_type, species2_ASTC, species2_AS, species2_event_type, filter, event_type){
  comparison <- paste0(species1,"_",species2)
  result <- df %>% 
                mutate(species1_ET = case_when(event_id %in% species1_ASTC$V1 ~ "ASTC",
                                                   event_id %in% species1_AS$V1 ~ "AS",
                                                   TRUE ~ "neither"),
                      species2_ET = case_when(event_id %in% species2_ASTC$V1 ~ "ASTC",
                                                   event_id %in% species2_AS$V1 ~ "AS",
                                                   TRUE ~ "neither"),
                      source_new = case_when(source == species1 ~ "species1",
                                             source == species2 ~ "species2")) %>% 
                group_by(current_snp, current_motif, Gene.name, source_new, species1_ET, species2_ET) %>% 
                summarize(total_score = sum(score)) %>% 
                filter(species1_ET == species1_event_type & species2_ET == species2_event_type) %>% 
                spread(source_new, total_score) %>% 
                mutate(species2 = replace_na(species2, 0),
                       species1 = replace_na(species1, 0),
                       species = comparison,
                       event_type = event_type) %>% 
                select(c(current_snp, current_motif, Gene.name, species1_ET, species2_ET, species1, species2, species, event_type))
  
 
  if(filter == TRUE){
    result <- result %>% 
      mutate(diff = species1-species2) %>% 
      filter(diff > 2 | diff < -2)
    
    colnames(result) <- c("current_snp", "current_motif", "Gene.name", paste0(species1,"_","ET"), paste0(species2,"_","ET"), species1, species2, "species","event_type","diff")
  } else{
    colnames(result) <- c("current_snp", "current_motif", "Gene.name", paste0(species1,"_","ET"), paste0(species2,"_","ET"), species1, species2, "species","event_type")
  }
  
  return(result)
  
}

# RBP by conservation
RBP_by_conservation <- function(df, species1, species2, conservation_df, conservation_status, filter, event_type){
  comparison <- paste0(species1,"_",species2)
  result <- df %>% 
    filter(event_id %in% conservation_df$event_id) %>% 
    mutate(source_new = case_when(source == species1 ~ "species1",
                                  source == species2 ~ "species2")) %>% 
    group_by(current_snp, current_motif, Gene.name, source_new) %>% 
    summarize(total_score = sum(score)) %>% 
    spread(source_new, total_score) %>% 
    mutate(species2 = replace_na(species2, 0),
           species1 = replace_na(species1, 0),
           species = comparison,
           event_type = event_type,
           cat = conservation_status) %>% 
    select(c(current_snp, current_motif, Gene.name, species1, species2, species, event_type, cat))
  
  
  if(filter == TRUE){
    result <- result %>% 
      mutate(diff = species1-species2) %>% 
      filter(diff > 2 | diff < -2)
    
    colnames(result) <- c("current_snp", "current_motif", "Gene.name", species1, species2, "species","event_type","cat","diff")
  } else{
    colnames(result) <- c("current_snp", "current_motif", "Gene.name", species1, species2, "species","event_type","cat")
  }
  
  return(result)
  
}

summary_table <- function(rnacmpt_ids, data_a, data_b, data_c, data_d, data_cons, data_spec){
  table_0 <- rnacmpt_ids %>% 
    mutate(Var1 = ID) %>% 
    arrange(Var1) %>% 
    select(Gene.name, Var1)
  
  table_a <- as.data.frame(table(data_a$current_motif))
  
  table_b <- as.data.frame(table(data_b$current_motif))
  
  table_c <- as.data.frame(table(data_c$current_motif))
  
  table_d <- as.data.frame(table(data_d$current_motif))
  
  data_e <- rbind(data_a, data_b)
  table_e <- as.data.frame(table(data_e$current_motif))
  
  table_f <- as.data.frame(table(data_cons$current_motif))
  
  table_g <- as.data.frame(table(data_spec$current_motif))
  
  rbp_table <- table_0 %>% 
    left_join(table_c, by = "Var1", copy = TRUE) %>% 
    left_join(table_b, by = "Var1", copy = TRUE) %>% 
    left_join(table_a, by = "Var1", copy = TRUE) %>% 
    left_join(table_e, by = "Var1", copy = TRUE) %>% 
    left_join(table_d, by = "Var1", copy = TRUE) %>% 
    left_join(table_f, by = "Var1", copy = TRUE) %>% 
    left_join(table_g, by = "Var1", copy = TRUE)
  
  colnames(rbp_table) <- c("Gene.name","RBP","ASTC_ASTC","AS_ASTC","ASTC_AS","AS_ASTC_total","AS_AS","ASTC_cons","ASTC_spec")
  
  rbp_table_long <- rbp_table %>% 
    replace_na(list(ASTC_ASTC =0,
                    AS_ASTC = 0,
                    ASTC_AS = 0,
                    AS_ASTC_total = 0,
                    ASTC_cons = 0,
                    ASTC_spec = 0)) %>% 
    gather(type, count, ASTC_ASTC:ASTC_spec) %>% 
    mutate(total = case_when(type == "ASTC_ASTC" ~ nrow(data_c),
                             type == "AS_ASTC" ~ nrow(data_b),
                             type == "ASTC_AS" ~ nrow(data_a),
                             type == "AS_ASTC_total" ~ nrow(data_e),
                             type == "AS_AS" ~ nrow(data_d),
                             type == "ASTC_cons" ~ nrow(data_cons),
                             type == "ASTC_spec" ~ nrow(data_spec)
    ))  %>% 
    group_by(RBP) %>% 
    mutate(chisq = case_when(type == "ASTC_ASTC" ~ 
                               prop.test(x=c(count[1],count[5]),n=c(total[1],total[5]))$p.value, 
                             type == "AS_ASTC" ~
                               prop.test(x=c(count[2],count[5]),n=c(total[2],total[5]))$p.value,
                             type == "ASTC_AS" ~
                               prop.test(x=c(count[3],count[5]),n=c(total[3],total[5]))$p.value,
                             type == "AS_ASTC_total" ~
                               prop.test(x=c(count[4],count[5]),n=c(total[4],total[5]))$p.value,
                             type == "ASTC_cons" ~ prop.test(x=c(count[6],count[5]), n=c(total[6],total[5]))$p.value,
                             type == "ASTC_spec" ~ prop.test(x=c(count[7],count[5]), n=c(total[7],total[5]))$p.value)) %>% 
    mutate(sig = case_when(chisq < 0.01 ~ "***",
                           chisq < 0.05 ~ "**",
                           chisq < 0.1 ~ "*")) %>% 
    mutate(frequency = count/total)
  
  return(rbp_table_long)
  
}

add_sig_level <- function(df, rbp_table, mode){
  if(mode == "total"){
    result <- df %>% 
      group_by(current_motif) %>% 
      left_join((subset(rbp_table, type == "AS_ASTC_total")), by = c("current_motif" = "RBP")) %>% 
      select(-"Gene.name.y","type","count","total","chisq","frequency")
  } else if(mode == "AS_ASTC"){
    result <- df %>% 
      group_by(current_motif) %>% 
      left_join((subset(rbp_table, type == "AS_ASTC")), by = c("current_motif" = "RBP")) %>% 
      select(-"Gene.name.y","type","count","total","chisq","frequency")
  } else if(mode == "ASTC_AS"){
    result <- df %>% 
      group_by(current_motif) %>% 
      left_join((subset(rbp_table, type == "ASTC_AS")), by = c("current_motif" = "RBP")) %>% 
      select(-"Gene.name.y","type","count","total","chisq","frequency")
  }
  return(result)
  
}



