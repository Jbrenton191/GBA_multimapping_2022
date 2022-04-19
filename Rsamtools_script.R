library(tidyverse)
library(optparse)
# library(Biostrings)
library(magrittr)
library(Rsamtools)
library(GenomicRanges)
library(IRanges)

bam_dir<-"./multi_map_tests/"
Gene_start<-155234452
Gene_end<-155244699
Gene_chromosome<-"chr1"
Gene_strand<-"-"
Pseudogene_start<-155213821
Pseudogene_end<-155227574
Pseudogene_chromosome<-"chr1"
Pseudogene_strand<-"-"

# Load Arguments ----
arguments <- parse_args(OptionParser(), positional_arguments = 9)

# bam_dir<-arguments$args[1]
# 
# Gene_start<-arguments$args[2]
# Gene_end<-arguments$args[3]
# Gene_strand<-arguments$args[4]
# Gene_chromosome<-arguments$args[5]

Gene_start<-as.numeric(Gene_start)
Gene_end<-as.numeric(Gene_end)
Gene_strand<-as.character(Gene_strand)
Gene_chromosome<-as.character(Gene_chromosome)

# Pseudogene_start<-arguments$args[6]
# Pseudogene_end<-arguments$args[7]
# Pseudogene_strand<-arguments$args[8]
# Pseudogene_chromosome<-arguments$args[9]

Pseudogene_start<-as.numeric(Pseudogene_start)
Pseudogene_end<-as.numeric(Pseudogene_end)
Pseudogene_strand<-as.character(Pseudogene_strand)
Pseudogene_chromosome<-as.character(Pseudogene_chromosome)

# Create df/tibble of bam files and their indexes ----

bam_list<-list.files(path = bam_dir, pattern = "*.bam$")
bam_list<-file.path(bam_dir, bam_list)

index_list<-c()
count<-1
for (i in 1:length(bam_list)) {
  index_list[count]<-paste0(bam_list[i], ".bai")
  count<-count+1
} 

bam_df<-tibble(bams=bam_list, indexes=index_list)

# Make a GRange object of the Gene and Pseudogene ----

gene_grange<-GRanges(seqnames = Gene_chromosome,
                     ranges = IRanges(start = Gene_start, end = Gene_end),
                     strand = Gene_strand)
pseudogene_grange<-GRanges(seqnames = Pseudogene_chromosome,
                           ranges = IRanges(start = Pseudogene_start, end = Pseudogene_end),
                           strand = Pseudogene_strand)

# For loop to determine overlaps within gene and pseudogene loci ----

Results_df<-tibble(name="", # name of file/bam
                   unique_gene_counts=NA, # how many unique/non-multimapping counts (any not just pairs) found within gene locus
                   unique_pseudogene_counts=NA, # how many counts (any not just pairs) found within pseudogene locus
                   multi_map_gene_counts=NA, # how many multimapping counts (any not just pairs) found within gene locus
                   multi_map_pseudogene_counts=NA, # how many multimapping counts (any not just pairs) found within pseudogene locus
                   gene_unique.vs.multi=NA, # ratio of: unique gene counts / multimapping gene counts
                   pseudogene_unique.vs.multi=NA, # ratio of: unique pseudogene counts / multimapping pseudogene counts
                   unique_gene_vs_pseudo=NA, # ratio of unique counts in gene vs unique counts in pseudo
                   multi_gene_vs_pseudo=NA, # ratio of multimap counts in gene vs multi counts in pseudo
                   multi_gene_in_pseudo.vs.gene_ratio=NA, # ratio of: multimapping reads within gene locus whose read 
                                                          # names are also present in pseudogene reads / all multimapping 
                                                          # gene counts
                   multi_pseudo_in_gene.vs.pseudo_ratio=NA, # ratio of: multimapping reads within pseudogene locus whose read 
                                                            # names are also present in gene reads / all multimapping reads within locus
                                                            # pseudogene counts
                   # multimap_gene_vs_other_ratio=NA, # ratio of read names of multimapping read within gene locus that map to any other location
                   # multimap_pseudogene_vs_other_ratio=NA # ratio of read names of multimapping read within pseudogene locus that map to any other location
                   )

# c(name, #name of file/bam
#   length(unique_map_gene_within_df$names), # 
#   length(unique_map_pseudogene_within_df$names),
#   length(multi_map_gene_within_df$names),
#   length(multi_map_pseudogene_within_df$names),
#   unique_count_ratio,
#   multi_count_ratio,
#   multi_gene_ratio
#   multi_pseudogene_ratio,
#   multimap_gene_vs_other_ratio,
#   multimap_pseudogene_vs_other_ratio)


for (i in 1:dim(bam_df)[1]) {
  
# Gene Counts/search ----
  
bam_file<-BamFile(file = bam_df$bams[i], index = bam_df$indexes[i])
p_G <- ScanBamParam(what = scanBamWhat(),
                   tag = "NH", 
                   which=gene_grange)
# does an Irange for ref/seq, start and end only. doesn't seem to filter anything but, start, end, width

filt_bam_gene <- scanBam(bam_file, param = p_G)

filt_bam_gene_df<-as_tibble(filt_bam_gene[[1]][c("qname", "flag", "rname", "strand", "pos", "qwidth", "cigar")])
NH<-as_tibble(filt_bam_gene[[1]]$tag)

filt_bam_gene_df<-cbind(filt_bam_gene_df, NH)

filt_bam_gene_df<-filt_bam_gene_df %>% filter(strand=="-")

unique_map_gene_df<-filt_bam_gene_df %>% filter(NH==1)
multi_map_gene_df<-filt_bam_gene_df %>% filter(NH>1)

# Make Grange objects of single and multi-mapped hits - although these are 
# still those that don't lie completely within the gene locus and overlap

# this also doesn't use cigar string, it merely uses read size as the end!!!!
# can adjust to make cigar string 

unique_map_gene_within_df<-unique_map_gene_df %>% filter(pos>=Gene_start & pos+qwidth<=Gene_end)
unique_map_gene_within_df$Gene_or_Pseudogene<-"Gene"
unique_map_gene_within_df$Mapping<-"Unique"

multi_map_gene_within_df<-multi_map_gene_df %>% filter(pos>=Gene_start & pos+qwidth<=Gene_end)
multi_map_gene_within_df$Gene_or_Pseudogene<-"Gene"
multi_map_gene_within_df$Mapping<-"Multi"

gene_searched_tibble<-rbind(unique_map_gene_within_df, multi_map_gene_within_df)

# Gene search done

# Pseudogene search ----

bam_file<-BamFile(file = bam_df$bams[i], index = bam_df$indexes[i])
p_PG <- ScanBamParam(what = scanBamWhat(),
                   tag = "NH", 
                   which=pseudogene_grange)

filt_bam_pseudogene <- scanBam(bam_file, param = p_PG)

filt_bam_pseudogene_df<-as_tibble(filt_bam_pseudogene[[1]][c("qname", "flag", "rname", "strand", "pos", "qwidth", "cigar")])
NH<-as_tibble(filt_bam_pseudogene[[1]]$tag)

filt_bam_pseudogene_df<-cbind(filt_bam_pseudogene_df, NH)

filt_bam_pseudogene_df<-filt_bam_pseudogene_df %>% filter(strand=="-")

unique_map_pseudogene_df<-filt_bam_pseudogene_df %>% filter(NH==1)
multi_map_pseudogene_df<-filt_bam_pseudogene_df %>% filter(NH>1)

# Make Grange objects of single and multi-mapped hits - although these are 
# still those that don't lie completely within the gene locus and overlap

unique_map_pseudogene_within_df<-unique_map_pseudogene_df %>% filter(pos>=Pseudogene_start & pos+qwidth<=Pseudogene_end)
unique_map_pseudogene_within_df$Gene_or_Pseudogene<-"Pseudogene"
unique_map_pseudogene_within_df$Mapping<-"Unique"

multi_map_pseudogene_within_df<-multi_map_pseudogene_df %>% filter(pos>=Pseudogene_start & pos+qwidth<=Pseudogene_end)
multi_map_pseudogene_within_df$Gene_or_Pseudogene<-"Pseudogene"
multi_map_pseudogene_within_df$Mapping<-"Multi"

pseudogene_searched_tibble<-rbind(unique_map_pseudogene_within_df, multi_map_pseudogene_within_df)

# Creating Ratios ----

# Ratio of unique maps ----

uniq_gene_in_pseudo<-length(unique_map_gene_within_df$qname[unique_map_gene_within_df$qname %in% unique_map_pseudogene_within_df$qname])
uniq_pseudo_in_gene<-length(unique_map_pseudogene_within_df$qname[unique_map_pseudogene_within_df$qname %in% unique_map_gene_within_df$qname])

uniq_gene_ratio<-uniq_gene_in_pseudo/length(unique_map_gene_within_df$qname)
uniq_pseudogene_ratio<-uniq_pseudo_in_gene/length(unique_map_pseudogene_within_df$qname)

# Ratio of multi maps ----
multi_gene_in_pseudo<-length(multi_map_gene_within_df$qname[multi_map_gene_within_df$qname %in% multi_map_pseudogene_within_df$qname])
multi_pseudo_in_gene<-length(multi_map_pseudogene_within_df$qname[multi_map_pseudogene_within_df$qname %in% multi_map_gene_within_df$qname])

multi_gene_ratio<-multi_gene_in_pseudo/length(multi_map_gene_within_df$qname)
multi_pseudogene_ratio<-multi_pseudo_in_gene/length(multi_map_pseudogene_within_df$qname)

# multi_mapping elsewhere

# p_g_other <- ScanBamParam(what = "qname")
# 
# bam_other_gene<-scanBam(bam_file, param = p_g_other)
# 
# bam_other_gene_df<-as_tibble(bam_other_gene[[1]][c("qname")])
# 
# bam_other_gene_df<-bam_other_gene_df$qname[bam_other_gene_df$qname %in% multi_map_gene_within_df$qname]
#   
# gene_other_multimap_num_reads<-length(bam_other_gene_df)
# 
# multimap_gene_vs_other_ratio<-length(multi_map_gene_within_df$qname)/gene_other_multimap_num_reads
# 
# # multimap ratio mapping to pseudogene vs elsewhere ----
# 
# p_pg_other <- ScanBamParam(what = "qname")
# 
# bam_other_pseudogene<-scanBam(bam_file, param = p_pg_other)
# 
# bam_other_pseudogene_df<-as_tibble(bam_other_pseudogene[[1]][c("qname")])
# 
# bam_other_pseudogene_df<-bam_other_pseudogene_df$qname[bam_other_pseudogene_df$qname %in% multi_map_gene_within_df$qname]
# 
# pseudogene_other_multimap_num_reads<-dim(bam_other_pseudogene_df)[1]
# 
# multimap_pseudogene_vs_other_ratio<-length(multi_map_pseudogene_within_df$qname)/pseudogene_other_multimap_num_reads

# Amount of names of reads multi-mapping to gene and elsewhere and same with pseudogene

# x<-"NM3350_PD683_A1B5_GM-T_S24_trimmed_mapped.BAM_Aligned.sortedByCoord.out.bam"
# name<-sub(x = x, pattern = ".bam", replacement = "")
# sub(bam_df$bams[i]
name<-sub(x = bam_df$bams[i], pattern = ".bam", replacement = "")

Results_df[i,]  <- as.list(c(name,
                   length(unique_map_gene_within_df$qname),
                   length(unique_map_pseudogene_within_df$qname),
                   length(multi_map_gene_within_df$qname),
                   length(multi_map_pseudogene_within_df$qname),
                   length(unique_map_gene_within_df$qname)/length(multi_map_gene_within_df$qname),
                   length(unique_map_pseudogene_within_df$qname)/length(multi_map_pseudogene_within_df$qname),
                   length(unique_map_gene_within_df$qname)/length(unique_map_pseudogene_within_df$qname),
                   length(multi_map_gene_within_df$qname)/length(multi_map_pseudogene_within_df$qname),
                   multi_gene_ratio,
                   multi_pseudogene_ratio))
}
                   
                   # ,
                   # multimap_gene_vs_other_ratio,
                   # multimap_pseudogene_vs_other_ratio
                   # ))

# }

# Want Ratio of total gene to pseudogene unique maps, multimaps and total within each sample
# Then go and check the multi-mapping reads and how many times they map elsewhere? 
# This encoded in the NH but does take into account mates?

# do venn diagram
# plyranges for granges for these objects

# use ensembl id to get all the 


# make csv
# write_csv(Results_df, "Multimapping_results.csv")



