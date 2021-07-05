# load libraries
library(ABHgenotypeR)
library(ASMap)
library(tidyverse)
library(googledrive)
library(googlesheets4)


# Import ABH calls
hun_ABH_prev_url <- "https://docs.google.com/spreadsheets/d/1fU3ZhOLmytAM1mBq5HkFOBM6Dolcv-42mDlAxDHjRj0/edit?usp=drive_web&ouid=101129460041007378390"

hun_ABH_prev <- hun_ABH_prev_url %>%
  as_id %>%
  range_read(., col_types = "c")

# Physical position of the chip markers
chip_map <- read_delim(
  "data/CHIP_genotyping/RSDNAMGM010819_raw.map", 
  "\t", 
  escape_double = FALSE, 
  col_names = FALSE, 
  trim_ws = TRUE) %>%
  rename_with(~c("chr", "marker", "gen_pos", "phys_pos")) %>%
  arrange(chr, phys_pos) %>%
  select(-gen_pos)

chip_marker_pos <- hun_ABH_prev %>%
  names %>%
  enframe(value = "marker") %>% 
  left_join(chip_map)%>%
  filter(!is.na(chr)) %>%
  mutate(marker2 = paste0(chr, "_", phys_pos))

# new crossobject with ABH calls and physical position of markers
hun_ABH_crossobject <- rbind(c("",chip_marker_pos$chr),
                             c("", chip_marker_pos$phys_pos)) %>%
  as_tibble %>%
  rename_with(~ c("id", chip_marker_pos$marker)) %>%
  bind_rows(hun_ABH_prev %>% filter(id != "NA")) %>%
  rename_with(~ c("id", chip_marker_pos$marker2)) 
  
# Export crossobject
write_csv(hun_ABH_crossobject, 
          "output/crossobject/hun_ABH_crossobject.csv")

# Import as crossobject into r-qtl
FUNGUS_CROSS <-read.cross(format = 'csv',
                          dir = "output/crossobject/",
                          file = 'hun_ABH_crossobject.csv',
                          estimate.map = FALSE,
                          F.gen = 2,
                          BC.gen = 0)

# Identify markers with the same pattern of segregation
# (duplicate or non-informative markers) 
fungus_dupmarkers <- FUNGUS_CROSS %>%
  findDupMarkers() %>%
  unlist(use.names = F) 

# remove cuplicate markers
FUNGUS_CROSS2 <- FUNGUS_CROSS %>%
  drop.markers(fungus_dupmarkers)

# Calculate first iteration of the genetic map using mstmap
FUNGUS_MAP <- mstmap(
  FUNGUS_CROSS2, 
  id = 'id', 
  bychr = TRUE, 
  dist.fun = 'kosambi',
  anchor = FALSE,
  detectBadData = TRUE)

# eliminate markers in linkage groups outside of chromosomes (.2)
fungus_marker_.2 <- FUNGUS_MAP %>%
  pull.map %>%
  map_df(
    .x =., 
    .f = ~ enframe(.x, name = "marker", value = "gpos"), 
    .id = "chr") %>%
  filter(grepl("\\.2$", chr)) %>%
  .$marker

# Remove markers 
FUNGUS_CROSS3 <- FUNGUS_CROSS2 %>%
  drop.markers(fungus_marker_.2)

#second iteration of genetic map
FUNGUS_MAP <- mstmap(
  FUNGUS_CROSS3, 
  id = 'id', 
  bychr = TRUE, 
  dist.fun = 'kosambi',
  anchor = T,
  detectBadData = TRUE)

# Identify markers with segregation distortion
gt <- geno.table(FUNGUS_MAP)

fungus_marker_weird <- rownames(gt[gt$P.value < 0.05/totmar(FUNGUS_MAP),])
fungus_marker_weird <- 
  fungus_marker_weird[! grepl("5_", fungus_marker_weird)]
# weird markers identified manually
fungus_marker_manual <- c(
  "2_221501477", "5_17708059", "5_31015101", "5_42255994", 
  "5_195455637", "6_89070446", "8_8103592", "8_60822982", 
  "10_38943265", "10_76422787"
  )

# Drop weir markers and markers with segregation distortion 
FUNGUS_CROSS4 <- FUNGUS_CROSS3 %>%
  drop.markers(c(fungus_marker_weird, fungus_marker_manual))

# Third iteration of the genetic map
FUNGUS_MAP <- mstmap(
  FUNGUS_CROSS4, 
  id = 'id', 
  bychr = TRUE, 
  dist.fun = 'kosambi',
  anchor = T,
  detectBadData = TRUE)

# Note: Genotypes AA: 25.2; AB: 507; BB: 24.1 . This genotype 
# frequencies are within expected for an F2 cross

# Physical vs genetic map (Marey map)
fungus_marey_map <- FUNGUS_MAP %>%
  pull.map %>%
  map_df(
    .x =., 
    .f = ~ enframe(.x, name = "marker", value = "gpos"), 
    .id = "chr") %>%
  mutate(chr = as.integer(chr)) %>%
  mutate(gpos = as.numeric(gpos)) %>%
  separate(marker, into = c(NA, "ppos"), remove = F) %>%
  mutate(ppos = as.integer(ppos)) %>%
  arrange(chr, ppos) %>%
  ggplot(data =., aes(x = ppos, y = gpos)) +
  geom_point() +
  facet_wrap(. ~ chr, scales = "free") +
  xlab("Physical position (MB)") +
  ylab("Genetic position (cM)") +
  ggtitle("CML312 x W22 HUN Marey maps")

dir.create(file.path(paste0(getwd(), "/output/plots"), 
                     "genetic_map"))

# Save marey_plots
ggsave(
  plot = fungus_marey_map,
  filename = "output/plots/genetic_map/Marey_plot.pdf",
  device = "pdf",
  scale = 2
  )

# Save genetic map

pdf("output/plots/genetic_map/CML312xW22_HUN_genmap.pdf") 
plotMap(FUNGUS_MAP)
# Close the pdf file
dev.off() 

# Save summary of the genetic map

sink("output/CML312xW22_HUN_genmap_summary.txt")
print(summary(FUNGUS_MAP))
sink()  # returns output to the console
  
# save crossobject
write.cross(
  cross = FUNGUS_MAP,
  format = "csv",
  filestem = "output/crossobject/cml312xw22_hun_genmap"
  )

# save genetic map in google drive
drive_upload(
  media = "output/crossobject/cml312xw22_hun_genmap.csv",
  path = input_data_drive_id,
  name = "cml312xw22_hun_genmap.csv",
  type = "spreadsheet",
  overwrite = T
)

# Identify markers removed
markers_removed <- c(
  fungus_marker_.2, 
  fungus_dupmarkers,
  fungus_marker_weird, 
  fungus_marker_manual
  ) %>%
  enframe(value = "marker", name = NULL) %>%
  separate(marker, into = c("chr", "pos"), remove = F) %>%
  mutate(across(c(chr, pos), ~ as.integer(.))) %>%
  arrange(chr, pos)
  




  
