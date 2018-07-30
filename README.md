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

# Citations and Resources

https://www.bioconductor.org/help/course-materials/2014/SeattleOct2014/B02.4_Annotation.html
https://bioconductor.org/packages/devel/bioc/vignettes/GenomicFeatures/inst/doc/GenomicFeatures.pdf
http://genomicsclass.github.io/book/pages/bioc1_annoOverview.html
https://bioconductor.org/packages/release/data/annotation/
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5408848/
https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.html
