# hm450hg38manifest
Manifest prepared for HM450 array with hg38 genome build

# Purpose
This repo simply stores a new HM450 manifest freshly compiled for hg38 using R/Bioc. 

# Code/Programming
This section outlines the steps used to create the manifest. In summary, we use R/Bioc. packages to load an existing hg19 HM450 manifest file, liftOver CpGs to updated positions for hg38, and then work extensively with Annotation DB objects and code to complete the annotation.

```{r}
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
man19 <- as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))
grman19 <- makeGRangesFromDataFrame(data.frame(chr=man19$chr,
                                               start=man19$pos,
                                               end=man19$pos,
                                               cg.id=man19$Name,
                                               stringsAsFactors = F),
                                    keep.extra.columns = T)
```

Use liftOver to change CpG locations for hg38.

```{r}
library(rtracklayer)
path = system.file(package="liftOver", "extdata", "hg19ToHg38.over.chain")
ch = import.chain(path)
grman38 <- liftOver(grman19, ch)
grman38 <- unlist(grman38)
length(grman38) # 485344, = 485344/485512 = 99.97%, 168 CGs not re-mapped/dropped by liftOver

# subset, sort, and assemble hg38 manifest with hg19 position info included
grman19.lo <- grman19[grman19$cg.id %in% grman38$cg.id]
identical(grman19.lo$cg.id,grman38$cg.id)
identical(seqnames(grman19.lo),seqnames(grman38))

man38 <- data.frame(chr=seqnames(grman38),
                    pos=start(grman38),
                    Name=grman38$cg.id,
                    chr.hg19=seqnames(grman19.lo),
                    pos.hg19=start(grman19.lo),
                    stringsAsFactors=F)
rownames(man38) <- man38$Name

# subset, sort, assemble the gene and islands annotation from hg19
man19.lo <- man19[man19$Name %in% man38$Name,]
identical(man19.lo$Name,man38$Name)
man38$strand.hg19 <- man19.lo$strand
man38$ucsc.refgene.group.hg19 <- man19.lo$UCSC_RefGene_Group
man38$ucsc.refgene.name.hg19 <- man19.lo$UCSC_RefGene_Name
man38$cpg.island.name.hg19 <- man19.lo$Islands_Name
man38$cpg.island.region.hg19 <- man19.lo$Relation_to_Island
```

Now load and subset the Annotation DBI objects from Bioc. We'll use TxDb hg38 and Ensembl v86/GRCh38 for this task.

```{r}
library(EnsDb.Hsapiens.v86);library(TxDb.Hsapiens.UCSC.hg38.knownGene)
tx38 <- TxDb.Hsapiens.UCSC.hg38.knownGene # hg38 build
ens38 <- EnsDb.Hsapiens.v86; seqinfo(ens38) # GRCh38 build

# get the gene features as GenomicRanges objects
genes.tx38 <- transcriptsBy(tx38,by="gene"); genes2.tx38 <- unlist(genes.tx38) # GenomicRangesList=>GenomicRanges object
pr.tx38 <- promoters(tx38)
tr.tx38 <- transcripts(tx38)
genes.ens38 <- genes(ens38);
pr.ens38 <- promoters(ens38)
tr.ens38 <- transcripts(ens38); 
# for ensembl db objects we need to fix the seqlevels labeling before querying
seqlevels(tx.ens38) <- paste0("chr",seqlevels(tx.ens38)) # # change seqlevels for subsetting "chr#"
seqlevels(pr.ens38) <- paste0("chr",seqlevels(pr.ens38))
seqlevels(genes.ens38) <- paste0("chr",seqlevels(genes.ens38))

```

Now, using findOverlaps, we quickly compile the annotation data into CpG-level data frames that we can then use to assemble the new manifest.

```{r}
# make a list of findOverlaps outputs referencing indices for gr objects
{
  fol.list <- list(findOverlaps(grman38,genes.ens38),
                   findOverlaps(grman38,genes2.tx38),
                   findOverlaps(grman38,tr.ens38),
                   findOverlaps(grman38,tr.tx38),
                   findOverlaps(grman38,pr.ens38),
                   findOverlaps(grman38,pr.tx38))
  names(fol.list) <- c("genes.ens",
                       "genes.tx",
                       "tr.ens",
                       "tr.tx",
                       "pr.ens",
                       "pr.tx")
  df.objlist <- list(grman38$cg.id,
                     genes.ens38$gene_id,
                     genes2.tx38$tx_id,
                     tr.ens38$tx_id,
                     tr.ens38$tx_biotype,
                     tr.ens38$gene_id,
                     tr.tx38$tx_id,
                     tr.tx38$tx_name,
                     pr.ens38$tx_id,
                     pr.ens38$tx_biotype,
                     pr.ens38$gene_id,
                     pr.tx38$tx_id,
                     pr.tx38$tx_name)
  names(df.objlist) <- c("cg.id",
                         "genes.ens",
                         "trgenes.tx",
                         "trid.ens",
                         "trbio.ens",
                         "trgene.ens",
                         "trid.tx",
                         "trname.tx",
                         "prtrid.ens",
                         "prbio.ens",
                         "prgene.ens",
                         "prtrid.tx",
                         "prtrname.tx")
}

df.list <- list() # list of dfs matching names using indices from findOverlaps query/subject
{
  # df genes ensembl
  df.list[[1]] <- data.frame(cg.id=df.objlist[[1]][queryHits(fol.list[[1]])],
                             geneid.ens=df.objlist[[2]][subjectHits(fol.list[[1]])],
                             stringsAsFactors = F)
  names(df.list)[1] <- "df.genes.ens"
  
  # df genes tx
  df.list[[2]] <- data.frame(cg.id=df.objlist[[1]][queryHits(fol.list[[2]])],
                             geneid.tx=df.objlist[[3]][subjectHits(fol.list[[2]])],
                             stringsAsFactors = F)
  names(df.list)[2] <- "df.genes.tx"
  
  # df tr ensembl
  df.list[[3]] <- data.frame(cg.id=df.objlist[[1]][queryHits(fol.list[[3]])],
                             trid.ens=df.objlist$trid.ens[subjectHits(fol.list[[3]])],
                             trbio.ens=df.objlist$trbio.ens[subjectHits(fol.list[[3]])],
                             trgene.ens=df.objlist$trgene.ens[subjectHits(fol.list[[3]])],
                             stringsAsFactors = F)
  names(df.list)[3] <- "df.tr.ens"
  
  # df tr tx
  df.list[[4]] <- data.frame(cg.id=df.objlist[[1]][queryHits(fol.list[[4]])],
                             trid.tx=df.objlist$trid.tx[subjectHits(fol.list[[4]])],
                             trname.tx=df.objlist$trname.tx[subjectHits(fol.list[[4]])],
                             stringsAsFactors = F)
  names(df.list)[4] <- "df.tr.tx"
  
  # df pr ens
  df.list[[5]] <- data.frame(cg.id=df.objlist[[1]][queryHits(fol.list$pr.ens)],
                             prid.ens=df.objlist$prtrid.ens[subjectHits(fol.list$pr.ens)],
                             prbio.ens=df.objlist$prbio.ens[subjectHits(fol.list$pr.ens)],
                             prgene.ens=df.objlist$prgene.ens[subjectHits(fol.list$pr.ens)],
                             stringsAsFactors = F)
  names(df.list)[5] <- "df.pr.tx"
  
  # df pr tx
  df.list[[6]] <- data.frame(cg.id=df.objlist[[1]][queryHits(fol.list$pr.tx)],
                             prtrid.tx=df.objlist$prtrid.tx[subjectHits(fol.list$pr.tx)],
                             prtrname.tx=df.objlist$prtrname.tx[subjectHits(fol.list$pr.tx)],
                             stringsAsFactors = F)
  names(df.list)[6] <- "df.pr.tx"
}


# new man file
man38new <- data.frame(cgid=grman38$cg.id,
                       chr=seqnames(grman38),
                       pos=start(grman38),
                       stringsAsFactors=F)
# manifest variables for genome annotation data
{
  man38new$geneid.ens <- "NA"
  man38new$geneid.tx <- "NA"
  man38new$trid.ens <- "NA"
  man38new$trbio.ens <- "NA"
  man38new$trgene.ens <- "NA"
  man38new$trid.tx <- "NA"
  man38new$trname.tx <- "NA"
  man38new$prid.ens <- "NA"
  man38new$prbio.ens <- "NA"
  man38new$prgene.ens <- "NA"
  man38new$prtrid.tx <- "NA"
  man38new$prtrname.tx <- "NA"
}

```

Now we're ready to populate the manifest using some data coercion with aggregate() function. We will also clean the data entries to remove extra characters and divide features using ';'.

```{r}
#===============================================================
# aggregate data frames by looping, then populate the manifest
#===============================================================
for(i in 1:length(df.list)){
  dfi <- df.list[[i]]
  coli <- intersect(colnames(dfi),colnames(man38new))
  for(x in 1:length(coli)){
    dfix <- aggregate(as.character(dfi[,coli[x]]) ~ cg.id, data = dfi, c)
    class(dfix[,2])
    dfix[,2] <- as.character(dfix[,2])
    #dfix[,2] <- unlist(dfix[,2]); class(dfix[,2])
    colnames(dfix)[2] <- coli[x]
    head(dfix)
    
    dfxnull <- data.frame(cg.id=man38new$cgid[!man38new$cgid %in% dfix$cg.id],
                          b=rep("NA",length(man38new$cgid[!man38new$cgid %in% dfix$cg.id])))
    colnames(dfxnull)[2] <- coli[x]
    
    dfix <- rbind(dfix,dfxnull)
    
    dim(dfix)
    dfix <- dfix[order(match(dfix$cg.id,man38new$cgid)),]
    identical(dfix$cg.id,man38new$cgid)
    
    man38new[,coli[x]] <- dfix[,2]
    
    message(i,":",x)
  }
  
  message(i)
}

man38new.final <- man38new

#=================================
# clear extra character and gsub
#=================================
col.oi <- colnames(man38new)[4:15]

for(i in 1:length(col.oi)){
  man38new[,col.oi[i]] <- gsub(",",";",gsub('c\\(|\\)|\\"|[[:space:]]',"",man38new[,col.oi[i]]))
  message(i)
}

# an example of what the above does
man38new[,col.oi[i]][994]
#[1] "c(\"ENSG00000123728\", \"ENSG00000232160\")"
gsub(",",";",gsub('c\\(|\\)|\\"|[[:space:]]',"",man38new[,col.oi[i]][994]))
#[1] "ENSG00000123728;ENSG00000232160"

# save the manifest!
save(man38new,file="hm450_hg38_manifest_new-from-findOL.rda")

```

# Session info for the above

- Session info -------------------------------------------------------------------------------------------------------
 setting  value                       
 version  R version 3.5.1 (2018-07-02)
 os       Windows >= 8 x64            
 system   x86_64, mingw32             
 ui       RStudio                     
 language (EN)                        
 collate  English_United States.1252  
 tz       America/Los_Angeles         
 date     2018-07-29                  

- Packages -----------------------------------------------------------------------------------------------------------
 package              * version    date       source                            
 AnnotationDbi          1.42.1     2018-05-08 Bioconductor                      
 AnnotationFilter       1.4.0      2018-05-01 Bioconductor                      
 assertthat             0.2.0      2017-04-11 CRAN (R 3.5.1)                    
 bindr                  0.1.1      2018-03-13 CRAN (R 3.5.1)                    
 bindrcpp               0.2.2      2018-03-29 CRAN (R 3.5.1)                    
 Biobase                2.40.0     2018-05-01 Bioconductor                      
 BiocGenerics         * 0.26.0     2018-05-01 Bioconductor                      
 BiocParallel           1.14.2     2018-07-09 Bioconductor                      
 biomaRt                2.36.1     2018-05-24 Bioconductor                      
 Biostrings             2.48.0     2018-05-01 Bioconductor                      
 bit                    1.1-14     2018-05-29 CRAN (R 3.5.0)                    
 bit64                  0.9-7      2017-05-08 CRAN (R 3.5.0)                    
 bitops                 1.0-6      2013-08-17 CRAN (R 3.5.0)                    
 blob                   1.1.1      2018-03-25 CRAN (R 3.5.1)                    
 clipr                * 0.4.1      2018-06-23 CRAN (R 3.5.1)                    
 clisymbols             1.2.0      2017-05-21 CRAN (R 3.5.1)                    
 crayon                 1.3.4      2017-09-16 CRAN (R 3.5.1)                    
 curl                   3.2        2018-03-28 CRAN (R 3.5.1)                    
 DBI                    1.0.0      2018-05-02 CRAN (R 3.5.1)                    
 DelayedArray           0.6.1      2018-06-15 Bioconductor                      
 devtools               1.13.6     2018-06-27 CRAN (R 3.5.1)                    
 digest                 0.6.15     2018-01-28 CRAN (R 3.5.1)                    
 dplyr                * 0.7.6      2018-06-29 CRAN (R 3.5.1)                    
 EnsDb.Hsapiens.v86     2.99.0     2018-07-25 Bioconductor                      
 ensembldb              2.4.1      2018-05-08 Bioconductor                      
 GenomeInfoDb         * 1.16.0     2018-05-01 Bioconductor                      
 GenomeInfoDbData       1.1.0      2018-07-22 Bioconductor                      
 GenomicAlignments      1.16.0     2018-05-01 Bioconductor                      
 GenomicFeatures        1.32.0     2018-05-01 Bioconductor                      
 GenomicRanges        * 1.32.6     2018-07-20 Bioconductor                      
 git2r                  0.23.0     2018-07-17 CRAN (R 3.5.1)                    
 glue                   1.3.0      2018-07-17 CRAN (R 3.5.1)                    
 hms                    0.4.2      2018-03-10 CRAN (R 3.5.1)                    
 httr                   1.3.1      2017-08-20 CRAN (R 3.5.1)                    
 IRanges              * 2.14.10    2018-05-16 Bioconductor                      
 lattice                0.20-35    2017-03-25 CRAN (R 3.5.1)                    
 lazyeval               0.2.1      2017-10-29 CRAN (R 3.5.1)                    
 magrittr               1.5        2014-11-22 CRAN (R 3.5.1)                    
 Matrix                 1.2-14     2018-04-13 CRAN (R 3.5.1)                    
 matrixStats            0.53.1     2018-02-11 CRAN (R 3.5.1)                    
 memoise                1.1.0      2017-04-21 CRAN (R 3.5.1)                    
 pillar                 1.3.0      2018-07-14 CRAN (R 3.5.1)                    
 pkgconfig              2.0.1      2017-03-21 CRAN (R 3.5.1)                    
 plyr                 * 1.8.4      2016-06-08 CRAN (R 3.5.1)                    
 prettyunits            1.0.2      2015-07-13 CRAN (R 3.5.1)                    
 progress               1.2.0      2018-06-14 CRAN (R 3.5.1)                    
 ProtGenerics           1.12.0     2018-05-01 Bioconductor                      
 purrr                  0.2.5      2018-05-29 CRAN (R 3.5.1)                    
 R6                     2.2.2      2017-06-17 CRAN (R 3.5.1)                    
 Rcpp                   0.12.17    2018-05-18 CRAN (R 3.5.1)                    
 RCurl                  1.95-4.11  2018-07-15 CRAN (R 3.5.1)                    
 rlang                  0.2.1      2018-05-30 CRAN (R 3.5.1)                    
 Rsamtools              1.32.2     2018-07-03 Bioconductor                      
 RSQLite                2.1.1      2018-05-06 CRAN (R 3.5.1)                    
 rstudioapi             0.7        2017-09-07 CRAN (R 3.5.1)                    
 rtracklayer            1.40.3     2018-06-02 Bioconductor                      
 S4Vectors            * 0.18.3     2018-06-08 Bioconductor                      
 sessioninfo          * 1.0.1.9000 2018-07-30 Github (r-lib/sessioninfo@c871d01)
 stringi                1.1.7      2018-03-12 CRAN (R 3.5.0)                    
 stringr                1.3.1      2018-05-10 CRAN (R 3.5.1)                    
 SummarizedExperiment   1.10.1     2018-05-11 Bioconductor                      
 tibble                 1.4.2      2018-01-22 CRAN (R 3.5.1)                    
 tidyselect             0.2.4      2018-02-26 CRAN (R 3.5.1)                    
 withr                  2.1.2      2018-03-15 CRAN (R 3.5.1)                    
 XML                    3.98-1.12  2018-07-15 CRAN (R 3.5.1)                    
 XVector                0.20.0     2018-05-01 Bioconductor                      
 zlibbioc               1.26.0     2018-05-01 Bioconductor                      

# Citations and Resources

https://www.bioconductor.org/help/course-materials/2014/SeattleOct2014/B02.4_Annotation.html
https://bioconductor.org/packages/devel/bioc/vignettes/GenomicFeatures/inst/doc/GenomicFeatures.pdf
http://genomicsclass.github.io/book/pages/bioc1_annoOverview.html
https://bioconductor.org/packages/release/data/annotation/
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5408848/
https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.html
