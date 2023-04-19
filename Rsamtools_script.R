library(tidyverse)
library(optparse)
# library(Biostrings)
library(magrittr)
library(Rsamtools)
library(GenomicRanges)
library(IRanges)

bam_dir<-"~/GBA_multimapping_2022/nextflow_pipeline/output/Samtools/"
Gene_start<-155234452
Gene_end<-155244699
Gene_chromosome<-"chr1"
Gene_strand<-"-"
Pseudogene_start<-155213821
Pseudogene_end<-155227574
Pseudogene_chromosome<-"chr1"
Pseudogene_strand<-"-"

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
                     ranges = IRanges(start = Gene_start, end = Gene_end))
pseudogene_grange<-GRanges(seqnames = Pseudogene_chromosome,
                           ranges = IRanges(start = Pseudogene_start, end = Pseudogene_end))

# does an Irange for ref/seq, start and end only. doesn't seem to filter anything but, start, end, width. Don't need strand - had in before - if breaks put back


Results_df<-tibble(name="", # name of file/bam
                   unique_gene_counts=NA, # how many unique/non-multimapping counts (any not just pairs) found within gene locus
                   unique_pseudogene_counts=NA, # how many counts (any not just pairs) found within pseudogene locus
                   multi_map_gene_counts=NA, # how many multimapping counts (any not just pairs) found within gene locus
                   multi_map_pseudogene_counts=NA, # how many multimapping counts (any not just pairs) found within pseudogene locus
                   gene_unique.vs.multi=NA, # ratio of: unique gene counts / multimapping gene counts
                   pseudogene_unique.vs.multi=NA, # ratio of: unique pseudogene counts / multimapping pseudogene counts
                   unique_gene_vs_pseudo=NA, # ratio of unique counts in gene vs unique counts in pseudo
                   multi_gene_vs_pseudo=NA, # ratio of multimap counts in gene vs unique counts in pseudo
                   multi_gene_in_pseudo.vs.gene_ratio=NA, # ratio of: multimapping reads within gene locus whose read names are also present in pseudogene reads / all multimapping gene counts
                   multi_pseudo_in_gene.vs.pseudo_ratio=NA # ratio of: multimapping reads within pseudogene locus whose read names are also present in gene reads / all multimapping pseudogene counts
                   #                  multimap_gene_vs_other_ratio=NA, # ratio of read names of multimapping read within gene locus that map to any other location
                   #                   multimap_pseudogene_vs_other_ratio=NA # ratio of read names of multimapping read within pseudogene locus that map to any other location
)

#---- Define Cigar extaction function

cig_extract_length<-function(x) {
  cig_vec<-as.vector(unlist(str_extract_all(string = x, pattern = '[0-9]*[A-Z|=]')))
  y<-cig_vec[which(grepl(x = cig_vec, pattern = "M|D|N|=|X"))]
  z<-str_extract_all(string = y, pattern = '[0-9]+')
  zz<-stringr::str_flatten(z, collapse = "+")
  x2<-parse(text = zz)
  eval(x2)
}

#---- For loop through data

for (i in 1:dim(bam_df)[1]) {
  
  # Gene Counts/search ----
  
  bam_file<-BamFile(file = bam_df$bams[i], index = bam_df$indexes[i])
  
  flag<-scanBamFlag(isPaired = T, isProperPair = T, isUnmappedQuery = F,  
                    hasUnmappedMate = F, isMinusStrand = T, isMateMinusStrand = F, 
                    isFirstMateRead = T, isSecondMateRead = F,
                    isSecondaryAlignment = NA, isNotPassingQualityControls = F, isDuplicate = F)
  p_G <- ScanBamParam(what = scanBamWhat(),
                      tag = "NH", 
                      which=gene_grange,
                      flag = flag)
  filt_bam_gene <- scanBam(bam_file, param = p_G)
  # does an Irange for ref/seq, start and end only. 
  # doesn't seem to filter anything but, start, end, width
  
  
  filt_bam_gene_df<-as_tibble(filt_bam_gene[[1]][c("qname", "flag", "rname", "strand", "pos", "qwidth", "cigar")])
  NH<-as_tibble(filt_bam_gene[[1]]$tag)
  
  filt_bam_gene_df<-cbind(filt_bam_gene_df, NH)
  
  filt_bam_gene_df<-filt_bam_gene_df %>% filter(strand=="-")
  
  
  
  cig<-""
  count<-0
  for (j in filt_bam_gene_df$cigar) {
    count<-count+1
    cig[count]<-cig_extract_length(j)
  }
  cig<-as.numeric(cig)
  
  filt_bam_gene_df$actual_width<-cig
  
  
  
  unique_map_gene_df<-filt_bam_gene_df %>% filter(NH==1)
  multi_map_gene_df<-filt_bam_gene_df %>% filter(NH>1)
  
  unique_map_gene_within_df<-unique_map_gene_df %>% filter(pos>=Gene_start & pos+actual_width<=Gene_end)
  unique_map_gene_within_df$Gene_or_Pseudogene<-"Gene"
  unique_map_gene_within_df$Mapping<-"Unique"
  
  multi_map_gene_within_df<-multi_map_gene_df %>% filter(pos>=Gene_start & pos+actual_width<=Gene_end)
  multi_map_gene_within_df$Gene_or_Pseudogene<-"Gene"
  multi_map_gene_within_df$Mapping<-"Multi"
  
  gene_searched_tibble<-rbind(unique_map_gene_within_df, multi_map_gene_within_df)
  
  
  # Pseudogene seach ----
  
  bam_file<-BamFile(file = bam_df$bams[i], index = bam_df$indexes[i])
  p_PG <- ScanBamParam(what = scanBamWhat(),
                       tag = "NH",
                       which=pseudogene_grange,
                       flag=flag)
  
  filt_bam_pseudogene <- scanBam(bam_file, param = p_PG)
  
  filt_bam_pseudogene_df<-as_tibble(filt_bam_pseudogene[[1]][c("qname", "flag", "rname", "strand", "pos", "qwidth", "cigar")])
  NH<-as_tibble(filt_bam_pseudogene[[1]]$tag)
  
  filt_bam_pseudogene_df<-cbind(filt_bam_pseudogene_df, NH)
  
  filt_bam_pseudogene_df<-filt_bam_pseudogene_df %>% filter(strand=="-")
  
  cig<-""
  count<-0
  for (j in filt_bam_pseudogene_df$cigar) {
    count<-count+1
    cig[count]<-cig_extract_length(j)
  }
  cig<-as.numeric(cig)
  filt_bam_pseudogene_df$actual_width<-cig
  
  unique_map_pseudogene_df<-filt_bam_pseudogene_df %>% filter(NH==1)
  multi_map_pseudogene_df<-filt_bam_pseudogene_df %>% filter(NH>1)
  
  unique_map_pseudogene_within_df<-unique_map_pseudogene_df %>% filter(pos>=Pseudogene_start & pos+actual_width<=Pseudogene_end)
  unique_map_pseudogene_within_df$Gene_or_Pseudogene<-"Pseudogene"
  unique_map_pseudogene_within_df$Mapping<-"Unique"
  
  multi_map_pseudogene_within_df<-multi_map_pseudogene_df %>% filter(pos>=Pseudogene_start & pos+actual_width<=Pseudogene_end)
  multi_map_pseudogene_within_df$Gene_or_Pseudogene<-"Pseudogene"
  multi_map_pseudogene_within_df$Mapping<-"Multi"
  
  pseudogene_searched_tibble<-rbind(unique_map_pseudogene_within_df, multi_map_pseudogene_within_df)
  
  uniq_gene_in_pseudo<-length(unique_map_gene_within_df$qname[unique_map_gene_within_df$qname %in% unique_map_pseudogene_within_df$qname])
  uniq_pseudo_in_gene<-length(unique_map_pseudogene_within_df$qname[unique_map_pseudogene_within_df$qname %in% unique_map_gene_within_df$qname])
  
  uniq_gene_ratio<-uniq_gene_in_pseudo/length(unique_map_gene_within_df$qname)
  uniq_pseudogene_ratio<-uniq_pseudo_in_gene/length(unique_map_pseudogene_within_df$qname)
  
  # Ratio of multi maps ----
  multi_gene_in_pseudo<-length(multi_map_gene_within_df$qname[multi_map_gene_within_df$qname %in% multi_map_pseudogene_within_df$qname])
  multi_pseudo_in_gene<-length(multi_map_pseudogene_within_df$qname[multi_map_pseudogene_within_df$qname %in% multi_map_gene_within_df$qname])
  
  multi_gene_ratio<-multi_gene_in_pseudo/length(multi_map_gene_within_df$qname)
  multi_pseudogene_ratio<-multi_pseudo_in_gene/length(multi_map_pseudogene_within_df$qname)
  
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


# Want Ratio of total gene to pseudogene unique maps, multimaps and total within each sample
# Then go and check the multi-mapping reads and how many times they map elsewhere? 
# This encoded in the NH but does take into account mates?

# do venn diagram
# plyranges for granges for these objects

# use ensembl id to get all the 


# make csv
write_csv(Results_df, "Multimapping_results.csv")



