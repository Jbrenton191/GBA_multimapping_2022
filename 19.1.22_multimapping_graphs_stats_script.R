library(tidyverse)
library(optparse)
# library(Biostrings)
library(magrittr)
# library(Rsamtools)
library(GenomicRanges)
library(IRanges)
library(wesanderson)
library("RColorBrewer")
library(MetBrewer)

setwd('~/multi_res.19.1.22/GBA/')
# x<-read_csv(file = "Multimapping_results_all_mates.csv")

res_df<-read_csv(file = "Multimapping_results.csv")

res_df$name<-str_extract(string = res_df$name, pattern = "PD[C][0-9]*|PD[0-9]*")

res_df$"Total reads mapping to the GBA locus (unique + multimapping)"<-res_df$unique_gene_counts+res_df$multi_map_gene_counts
res_df$"Total reads mapping to the GBAP1 locus (unique + multimapping)"<-res_df$unique_pseudogene_counts+res_df$multi_map_pseudogene_counts

res_df %<>% dplyr::rename("Reads mapping to the GBA locus uniquely" = unique_gene_counts)
res_df %<>% dplyr::rename("Reads mapping to the GBAP1 locus uniquely" = unique_pseudogene_counts)
res_df %<>% dplyr::rename("Reads multi-mapping to the GBA locus" = multi_map_gene_counts)
res_df %<>% dplyr::rename("Reads multi-mapping to the GBAP1 locus" = multi_map_pseudogene_counts)

res_df$"% of all reads mapping to the GBA locus that map uniquely"<-(res_df$`Reads mapping to the GBA locus uniquely`/res_df$`Total reads mapping to the GBA locus (unique + multimapping)`)*100

res_df$"% of all reads mapping to the GBAP1 locus that map uniquely"<-(res_df$`Reads mapping to the GBAP1 locus uniquely`/res_df$`Total reads mapping to the GBAP1 locus (unique + multimapping)`)*100


res_df %<>% dplyr::rename("% Reads multi-mapping to the GBA locus who also multi-map to the GBAP1 locus (ie the read name is also present in multi-mapping GBAP1 reads)" = multi_gene_in_pseudo.vs.gene_ratio)
res_df %<>% dplyr::rename("% Reads multi-mapping to the GBAP1 locus who also multi-map to the GBA locus (ie the read name is also present in multi-mapping GBA reads)" = multi_pseudo_in_gene.vs.pseudo_ratio)

# Adding group variable

metadata_file<-"../20201229_MasterFile_SampleInfo.csv"
# metadata_key<-arguments$args[2]

meta_df<-read_csv(file = metadata_file, col_names = T, skip=1)
name_df<-res_df$name

# name_vec<-as.vector(unlist(name_df))

key<-which(meta_df$CaseNo %in% name_df)

selected_metadata<-meta_df[key,]

selected_metadata<-selected_metadata[,c(1:2)]

selected_metadata<-unique(selected_metadata)

names(selected_metadata)[1]<-"name"

x<-left_join(res_df, selected_metadata, by="name")

x %<>% select(name, Disease_Group, `Reads mapping to the GBA locus uniquely`, `Reads multi-mapping to the GBA locus`, `Total reads mapping to the GBA locus (unique + multimapping)`, `% of all reads mapping to the GBA locus that map uniquely`, `% Reads multi-mapping to the GBA locus who also multi-map to the GBAP1 locus (ie the read name is also present in multi-mapping GBAP1 reads)`, 
             `Reads mapping to the GBAP1 locus uniquely`, `Reads multi-mapping to the GBAP1 locus`, `Total reads mapping to the GBAP1 locus (unique + multimapping)`, `% of all reads mapping to the GBAP1 locus that map uniquely`,`% Reads multi-mapping to the GBAP1 locus who also multi-map to the GBA locus (ie the read name is also present in multi-mapping GBA reads)`)

write_csv(file = "Multimapping_results_groups_added.csv", x = x)

multimap_full_df<-x

library(car)

a<-car::Anova(glm(data=multimap_full_df, formula = `% of all reads mapping to the GBA locus that map uniquely` ~ Disease_Group), 
           type=2, test="F")

b<-car::Anova(glm(data=multimap_full_df, formula = `% Reads multi-mapping to the GBA locus who also multi-map to the GBAP1 locus (ie the read name is also present in multi-mapping GBAP1 reads)` ~ Disease_Group), 
           type=2, test="F")

c<-car::Anova(glm(data=multimap_full_df, formula = `% of all reads mapping to the GBAP1 locus that map uniquely` ~ Disease_Group), 
           type=2, test="F")

d<-car::Anova(glm(data=multimap_full_df, formula = `% Reads multi-mapping to the GBAP1 locus who also multi-map to the GBA locus (ie the read name is also present in multi-mapping GBA reads)` ~ Disease_Group), 
           type=2, test="F")

knitr::kable(a, digits = 3, caption = 'ANOVA (glm and F Test (type 2 SS)) of: % of all reads mapping to the GBA locus that map uniquely', format = 'html')
# Then copied output and made html file using vim
knitr::kable(b, digits = 3, caption = 'ANOVA (glm and F Test (type 2 SS)) of: % Reads multi-mapping to the GBA locus who also multi-map to the GBAP1 locus', format = 'html')

knitr::kable(c, digits = 3, caption = 'ANOVA (glm and F Test (type 2 SS)) of: % of all reads mapping to the GBAP1 locus that map uniquely', format = 'html')

knitr::kable(d, digits = 3, caption = 'ANOVA (glm and F Test (type 2 SS)) of: % Reads multi-mapping to the GBAP1 locus who also multi-map to the GBA locus', format = 'html')


# x$`% of all reads mapping to the GBA locus that map uniquely`
# x$`% Reads multi-mapping to the GBA locus who also multi-map to the GBAP1 locus (ie the read name is also present in multi-mapping GBAP1 reads)`
# x$`% of all reads mapping to the GBAP1 locus that map uniquely`
# x$`% Reads multi-mapping to the GBAP1 locus who also multi-map to the GBA locus (ie the read name is also present in multi-mapping GBA reads)`


png(filename = "percent_GBA_unique_mapping.png",
    type = "cairo", width = 4000, height = 2500, res = 300)

ggplot(multimap_full_df, aes(Disease_Group, `% of all reads mapping to the GBA locus that map uniquely`,
              colour=Disease_Group))+
  # geom_violin()+
  geom_boxplot(outlier.shape = NA)+
    geom_point()+
  scale_color_manual(values=met.brewer("Egypt", n=4))+
  coord_flip()+
  guides(color = guide_legend(reverse = T))+
  labs(x = "Disease Group", y="% Reads (unique/(unique+multi-mapped))")+
  ggtitle('Percentage of all reads mapping to the GBA locus uniquely')+
  theme(plot.title = element_text(size = 24,face = "bold.italic", hjust = 0.5, vjust = 0), 
        plot.subtitle = element_text(size = 14, face = "italic"),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(face = "bold", size = 14, margin = margin(3)),
        axis.title.y = element_text(face = "bold", size = 14, vjust=1.5, margin = margin(3)),
        legend.text = element_text(size = 12),
        legend.position = "right",
        legend.title = element_text(size=14),
        legend.spacing = unit(x = 0, units = "pt"), 
        legend.box.spacing = margin(0.0001))

dev.off()

png(filename = "percent_GBA_multimapping_to_GBAP1.png",
    type = "cairo", width = 4000, height = 2500, res = 300)

ggplot(multimap_full_df, aes(Disease_Group, 
                             `% Reads multi-mapping to the GBA locus who also multi-map to the GBAP1 locus (ie the read name is also present in multi-mapping GBAP1 reads)`*100,
                             colour=Disease_Group))+
  # geom_violin()+
  geom_boxplot(outlier.shape = NA)+
  geom_point()+
  scale_color_manual(values=met.brewer("Egypt", n=4))+
  coord_flip()+
  guides(color = guide_legend(reverse = T))+
  labs(x = "Disease Group", y="% Reads multi-mapping to GBAP1")+
  ggtitle('Percentage of all reads multi-mapping to the GBA locus who also multi-map to the GBAP1 locus')+
  theme(plot.title = element_text(size = 20,face = "bold.italic", hjust = 0.5, vjust = 0), 
        plot.subtitle = element_text(size = 14, face = "italic"),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(face = "bold", size = 14, margin = margin(3)),
        axis.title.y = element_text(face = "bold", size = 14, vjust=1.5, margin = margin(3)),
        legend.text = element_text(size = 12),
        legend.position = "right",
        legend.title = element_text(size=14),
        legend.spacing = unit(x = 0, units = "pt"), 
        legend.box.spacing = margin(0.0001))

dev.off()

png(filename = "percent_GBAP1_unique_mapping.png",
    type = "cairo", width = 4000, height = 2500, res = 300)

ggplot(multimap_full_df, aes(Disease_Group, `% of all reads mapping to the GBAP1 locus that map uniquely`,
                             colour=Disease_Group))+
  # geom_violin()+
  geom_boxplot(outlier.shape = NA)+
  geom_point()+
  scale_color_manual(values=met.brewer("Egypt", n=4))+
  coord_flip()+
  guides(color = guide_legend(reverse = T))+
  labs(x = "Disease Group", y="% Reads (unique/(unqiue+multi-mapped))")+
  ggtitle('Percentage of all reads mapping to the GBAP1 locus uniquely')+
  theme(plot.title = element_text(size = 24,face = "bold.italic", hjust = 0.5, vjust = 0), 
        plot.subtitle = element_text(size = 14, face = "italic"),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(face = "bold", size = 14, margin = margin(3)),
        axis.title.y = element_text(face = "bold", size = 14, vjust=1.5, margin = margin(3)),
        legend.text = element_text(size = 12),
        legend.position = "right",
        legend.title = element_text(size=14),
        legend.spacing = unit(x = 0, units = "pt"), 
        legend.box.spacing = margin(0.0001))

dev.off()

png(filename = "percent_GBAP1_multimapping_to_GBA.png",
    type = "cairo", width = 4000, height = 2500, res = 300)

ggplot(multimap_full_df, aes(Disease_Group, 
                             `% Reads multi-mapping to the GBAP1 locus who also multi-map to the GBA locus (ie the read name is also present in multi-mapping GBA reads)`*100,
                             colour=Disease_Group))+
  # geom_violin()+
  geom_boxplot(outlier.shape = NA)+
  geom_point()+
  scale_color_manual(values=met.brewer("Egypt", n=4))+
  coord_flip()+
  guides(color = guide_legend(reverse = T))+
  labs(x = "Disease Group", y="% Reads multi-mapping to GBA")+
  ggtitle('Percentage of all reads multi-mapping to the GBAP1 locus who also multi-map to the GBA locus')+
  theme(plot.title = element_text(size = 20,face = "bold.italic", hjust = 0.5, vjust = 0), 
        plot.subtitle = element_text(size = 14, face = "italic"),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(face = "bold", size = 14, margin = margin(3)),
        axis.title.y = element_text(face = "bold", size = 14, vjust=1.5, margin = margin(3)),
        legend.text = element_text(size = 12),
        legend.position = "right",
        legend.title = element_text(size=14),
        legend.spacing = unit(x = 0, units = "pt"), 
        legend.box.spacing = margin(0.0001))

dev.off()






