#map stats on remapped primate data

## sequencing statistics

map_stats <- read.delim("~/Documents/05_results/01_primates/andrew_pipeline/remap_01_2020/remap_map_stats.csv", 
                        sep = ",", 
                        stringsAsFactors = F) 

all_samples <- map_stats$X.Sample[1:36]

map_stats_new <- map_stats %>% 
  rowwise() %>% 
  separate(BothUCnt, into = c("bothucnt","m1"), sep ="M") %>% 
  separate(BothSCnt, into = c("junction_reads","m2"), sep ="M") %>% 
  mutate_at("junction_reads", as.numeric) %>% 
  mutate_at("bothucnt", as.numeric) %>% 
  mutate(uniquely_mapped = bothucnt - junction_reads) %>% 
  dplyr::select(X.Sample, junction_reads, uniquely_mapped) %>% 
  gather(key="group", value="count", junction_reads, uniquely_mapped) %>% 
  rowwise() %>% 
  mutate(species = case_when(str_detect(X.Sample,"hips") ~ "human",
                             str_detect(X.Sample,"cips") ~ "chimpanzee",
                             str_detect(X.Sample,"oips") ~ "orangutan")) %>% 
  separate(X.Sample, into = c("fake_species","fraction"), sep = "s_") %>% 
  filter(!str_detect(fraction, "untr"))

all_fractions <- map_stats_new$fraction[1:30]
#map_stats_new$species <- rep(c(rep("human",10),rep("chimpanzee",10),rep("orangutan",10)),2)
map_stats_new$species <- factor(map_stats_new$species, levels = c("human","chimpanzee","orangutan"))
map_stats_new$fraction <- factor(map_stats_new$fraction, levels = c("cyto_jd952","cyto_jd958","mono_jd954",
                                                                    "mono_jd960","poll_jd955","poll_jd961",
                                                                    "polm_jd956","polm_jd962","polh_jd957",
                                                                    "polh_jd963","cyto_jd940","cyto_jd946",
                                                                    "mono_jd942","mono_jd948","poll_jd943",
                                                                    "poll_jd949","polm_jd944","polm_jd950",
                                                                    "polh_jd945","polh_jd951","cyto_jd928",
                                                                    "cyto_jd934","mono_jd930","mono_jd936",
                                                                    "poll_jd931","poll_jd937","polm_jd932",
                                                                    "polm_jd938","polh_jd933","polh_jd939"))
pdf("~/Documents/05_results/01_primates/andrew_pipeline/remap_01_2020/remap_map_stats.pdf", 
    width = 8, height = 4)
ggplot(map_stats_new, aes(x=fraction, y=count, fill=group)) +
  geom_col(color = "black", width = 0.85) +
  scale_fill_manual(values = c("white","grey36"),
                    labels = c("junction reads","uniquely mapped reads")) +
  ylab("Read Count in M") +
  xlab("") +
  facet_grid(cols = vars(species), scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  #scale_x_discrete(labels=rep(c("rep1","rep2"),5)) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x.bottom = element_text(vjust = 0.5))
dev.off()
