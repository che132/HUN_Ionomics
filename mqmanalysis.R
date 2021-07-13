################################
# Multi-QTL analysis using MQM #
################################

MQM_ion_crossobject <-  read.cross(
  format = "csv",
  dir = "output/crossobject",
  file = "crossobject_ion_hun.csv",
  genotypes = c("AA", "AB", "BB"),
  alleles = c("CML", "W22"),
  estimate.map = F)

# creating a pheno.col dataframe
MQM_ion_pheno.col <- MQM_ion_crossobject %>%
  phenames %>%
  enframe(name = "pheno.col", value = "phenotype") %>%
  filter(row_number() > 2) 

# genotype imputation (augmentation)
MQM_ion_hun_add <- mqmaugment(MQM_ion_hun_add, minprob = 0.1)

# Set automatic cofactors across all genome (100)
MQM_ion_hun_add_cofactors <- mqmautocofactors(MQM_ion_hun, 100)


# Multiqtl scan using MQM
MQM_ion_hun_add_scan <- MQM_ion_pheno.col %>%
  mutate(
    mqmscan = 
      map(.x = pheno.col,
          .f = ~ mqmscan(
            cross = MQM_ion_hun_add, 
            cofactors = MQM_ion_hun_add_cofactors,
            pheno.col = .x, 
            step.size = 1,
            window.size = 25,
            n.clusters = 4)))

# Permutation test (1000 perms)
MQM_ion_hun_add_scan_perms <- MQM_RIL_scan_perms

MQM_ion_hun_add_scan_perms <- MQM_ion_hun_add_scan %>%
  mutate(perms = purrr::map(.x = pheno.col,
                            .f = ~ mqmpermutation(cross = MQM_ion_hun_add, 
                                                  pheno.col = .x, 
                                                  method = "permutation",
                                                  n.cluster = 4,
                                                  n.perm = 1000,
                                                  verbose = F)))

# save scanone
save(
  list = c("ion_crossobject_intcov.out.hk"), 
  file = "output/qtl_objects/scanone/MQM_ion_hun_add_scan_perms.RData"
)


MQM_RIL_scan_perms_2 <- MQM_RIL_scan_perms %>%
  mutate(perms_2 = purrr::map(.x = perms, .f = ~ mqmprocesspermutation(.x)))

MQM_RIL_QTL_summary <- MQM_RIL_scan_perms_2 %>%
  MQM_summary_to_tibble(., CROSS = MQM_ion_hun_add)

MQM_ion_hun_add_LOD_map <- MQM_RIL_scan_perms_2 %>%
  mutate(
    mapa = map(
      .x = mqmscan,
      .f = ~ mapa_base_MQM(
        SCANONE = .x, 
        TIBBLE = MQM_RIL_QTL_summary))) %>%
  .$mapa %>%
  Reduce(bind_rows, .)

MQM_lodplot_function <- function(data_df, ion) {
  
  title <- paste0("Ion: ", ion)
  
  plot <- data_df %>% 
    {
      ggplot(data =., 
             aes(x = cM, 
                 y = lod, 
                 group = interaction(chr), 
                 alpha = alpha)) +
        geom_line(size = 1) +
        geom_line(aes(x = cM, 
                      y = lod, 
                      group = (chr)), 
                  alpha = 0.15, 
                  size = 1.5) +
        geom_vline(xintercept = c(.$inicio, .$fin), 
                   linetype = "dotted") +
        geom_hline(aes(yintercept = treshold),
                   linetype = "dotted", 
                   size = 0.75, 
                   alpha = 0.75) +
        theme(panel.background = 
                element_rect(
                  fill = "white", 
                  colour = "grey50")) +
        scale_x_continuous(
          name = "chr", 
          breaks = unique(.$breaks), 
          labels = unique(as.character((.$chr)))) +
        scale_alpha_identity() + 
        theme(legend.position = "top",
              text = element_text(size = 12), 
              plot.title = element_text(hjust = 0.5)) +
        xlab(NULL) +
        ylab("LOD") +
        ggtitle(title) +
        facet_grid(trait ~ .) 
        
    }
  
  return(plot)
}

MQM_ion_hun_add_lodplots <- MQM_ion_hun_add_LOD_map %>%
  group_by(trait) %>%
  nest %>%
  separate(trait, into = c("ion", NA), remove = F) %>%
  unnest(data) %>%
  group_by(ion) %>%
  nest %>%
  mutate(plot = map2(.x = data,
                     .y = ion, 
                     .f = ~ MQM_lodplot_function(.x, .y))) 
  
p <- MQM_ion_hun_add_lodplots$plot

# create directory for plots
# dir.create(file.path(paste0(getwd(), "/output/plots"), "MQM_lodplots"))

# Saving plots into a single pdf
ggsave(
  filename = 
    "output/plots/MQM_lodplots/ion_MQM_lodplots.pdf", 
  plot = gridExtra::marrangeGrob(p, nrow=1, ncol=1), 
  width = 15, height = 9
)




