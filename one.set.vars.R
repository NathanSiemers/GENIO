## http and https proxies?
use.proxy = TRUE
if(use.proxy){
    Sys.setenv(http_proxy='http://proxy-server.hpw.pri.bms.com:8080')
    Sys.setenv(https_proxy='http://proxy-server.hpw.pri.bms.com:8080')
}

my.cores = 5
fast.run = TRUE

if (fast.run) {
    skip.early = TRUE
    use.cached = TRUE
    compute.anova = FALSE
    do.download = FALSE
    bigcompute = FALSE
    compute.eval = TRUE
} else {
    ## skip all early loadings?
    skip.early = FALSE
    ## used cached R objects (files) instead of re-loading
    use.cached = FALSE
    ## compute the anovas or use cached?  they are getting expensive...
    compute.anova = TRUE
    ## download fresh data from ucsc?
    do.download = TRUE
    ## topnumber - top N mutations and CNVs to study
    ## use parallel compute power?
    bigcompute=TRUE
    if(bigcompute) {
        library(doParallel)
        registerDoParallel(cores=my.cores)
    }
}

topmut.cutoff = 3200
## restricted set of cohorts for paper
paper.cohorts.only = TRUE


## force (re)creation of integrated data objects even if objects exist?
force.reload = FALSE
## create 4-way clin-rna-mut-cnv data set?
## not necessary for this analysis, but generally useful
create.allmerge = FALSE
## created cohort reports?  costs some time and memory
create.dreports = FALSE
## generate supplementary tables? (uses ddply, a bit slow)
make.suptable = FALSE
## dir where data and results are read and written
data.dir = "data/"
## knitit output root filename
filename = "manuscript"
## adjust P values
adjust.P = TRUE
## set to true unless everything is already loaded into global R env
##compute.eval = TRUE




## function  you will need to run 'knitit'
library(knitr)

knitit = function(filename = filename, withgenes=FALSE, nerd.report=FALSE, myeval=TRUE, format='landscape', centralreport=withgenes) {
    system('mkdir Reports')
    if (format == 'landscape') {
        figwidth=10; figheight=5; docx = 'manuscript.docx'
    }
    if (format == 'portrait') {
        figwidth=7.5; figheight=5; docx = 'manuscript.docx'
    }

    if(nerd.report) {myresults = 'show';myecho = TRUE;mywarning = TRUE;mymessage = TRUE}
    if(!nerd.report) {myresults = 'hide'; myecho=FALSE; mywarning = FALSE; mymessage=FALSE}
    if(withgenes == FALSE) { genestring = "" }
    if(! withgenes == FALSE) {genestring = paste(withgenes, sep=".", collapse=".") }
    if(!nerd.report) nerdstring = ""
    if(nerd.report) nerdstring = "nerd."
    bigcompute=0
    rmdfile = paste(filename, ".Rmd", sep="")
    print(rmdfile)
    mdfile = paste(genestring , filename, ".md", sep="")
    print(mdfile)
    docxfile = paste("Reports/", genestring,  filename, ".", nerdstring,  Sys.Date(),  ".docx", sep="")
    centralfile = paste("~/Reports/", paste(genestring), "/", genestring,  filename, ".", nerdstring,  Sys.Date(),  ".docx", sep="")
    print(centralfile)
    centraldir = paste("~/Reports/", genestring, sep="")
    print(centraldir)
    print(docxfile)
    opts_chunk$set(comment=NA, fig.cap="", fig.keep='all', fig.show='asis', fig.width=figwidth, fig.height=figheight, dpi=300, eval=myeval, results=myresults, echo=myecho, warning=mywarning, message=mymessage)
    knit(rmdfile, output=mdfile)
    mystring = paste('pandoc',    " --reference-docx=", docx, " ", mdfile, ' -f markdown -t docx -o ', docxfile, sep="")
    system(mystring)
    if(! centralreport == FALSE) {
        mystring = paste('pandoc',    " --reference-docx=", docx, " ", mdfile, ' -f markdown -t docx -o ', centralfile, sep="")
        dirstring = paste("mkdir", centraldir)
        print(dirstring)
        system(dirstring)
        system(mystring)
        }
}


knitr::opts_chunk$set(fig.process = function(x) {
  x2 = sub('-\\d+([.][a-z]+)$', '\\1', x)
  if (file.rename(x, x2)) x2 else x
})


opts_knit$set(unnamed.chunk.label = "Fig")

## at this point (but see some more variables below) you can run
##
##      knitit(filename)
##

## make sure pandoc is installed on the server
## and in your path

## code will be executed, but see flags below
## markdown and docx files will be created
## word document will be in ./Reports directory and date stamped

## only used when my drive needs to be remounted on my mac
##  setwd('/stf/biosrc/siemersn/StorEDGE/TCGA/PANCAN/11.15')

## use server network methods? our servers use different methods than  my mac
##server.network=TRUE


paper.cohorts.not = c('Large_B-cell_Lymphoma', 'thymoma', 'Glioblastoma', 'Lower_Grade_Glioma', 'Acute_Myeloid_Leukemia')

## list of cohorts we don't want included (substring match)
nonpaper.cohorts = c("ilot", "eukemia", "ymphoma", "hymoma", "glioma", "Glioma")


## end of flags you need to worry about when/before running...
################################################################

## remove intermediate data frames to save memory?
remove.intermediate = FALSE

## leave heme in for now...
bad.cohorts = c("ilot")

## make gistic data frame for manhattan plots?
make.gistic = FALSE

## CNA gain/loss cutoff
## gain calculation will be x > -1 * cna.cutoff
## loss calculation will be x < cna.cutoff
cna.cutoff = 0.1


## minimum correlation of signatures to use for study
minimum.sig.coherence = 0.45

## minimum cohort size to include in calculations
min.samples = 20

## minimum number of occurences of specific mutations
fmut.min = 2

## P cutoff to show effect size in graphs
effect.cutoff.P = 3

## default P value threshold for viewing on manhattan plots
v.p.view.threshold = 1


## this is now the multiple test correction for CNVs
## clustering of cnv data, single linkage
## cutoff 0.05 (> 0.95 correlation is same cluster)
## number of clusters = 1023
topcnv.cutoff = 1023

## sigs to include in manuscript
## if I've coded this right, all sigs will be used in calculation....
## only the ones below will show up in the manuscript
good.sigs = c("TCD8.sig", "Treg.sig", "Tcell.sig",  "Bcell.sig",
    "Mono.sig", "MGran.sig", "MFm2.sig", "NK.sig", "TregCD8.sig", "NKCD8.sig")

## remove mgran for paper?? YES
paper.sigs = c("TCD8.sig", "Treg.sig", "Tcell.sig",  "Bcell.sig",
    "Mono.sig",  "MFm2.sig", "NK.sig", "TregCD8.sig", "NKCD8.sig")


figure.number = table.number = table.number.printed = figure.number.printed = stable.number = stable.number.printed = sfigure.number = sfigure.number.printed = 0


## my standard includes
library(ggplot2)
library(gplots)
library(grid)
library(gridExtra)
library(GGally)
library(ggthemes)
library(RColorBrewer)
library(plyr)
library(reshape2)
library(data.table)
library(RJSONIO)
library(limma)
library(RCurl)
library(org.Hs.eg.db)
library(GO.db)


################################################################
################################################################
################################################################
## functions

plus = function ( var.as.string ) {
    x = get(var.as.string, envir = .GlobalEnv) + 1
    assign(var.as.string, x, envir = .GlobalEnv)
    x
}

decon_signature_list = list(
    TCD8.sig = c("CD8A", "CD8B"),
    Treg.sig = c("FOXP3", "CCR8"),
    Tcell.sig = c("CD3D", "CD3E", "CD2"),
    Bcell.sig = c("CD19", "CD79A", "MS4A1"),
    NK.sig = c( "KIR2DL1" , "KIR2DL3" , "KIR2DL4" , "KIR3DL1" , "KIR3DL2" , "KIR3DL3" , "KIR2DS4" ),
    MGran.sig = c("CLEC4D", "CLEC4E", "CLEC6A"),
    Mono.sig = c("CD86", "CSF1R", "C3AR1"),
    MFm2.sig = c("CD163", "VSIG4", "MS4A4A"),
    Fib.sig = c("COL1A1", "COL1A2", "COL3A1"),
    IFNG.sig = c("IFNG")
    )
dsl = decon_signature_list

decon_genelist = unlist(decon_signature_list)

deconfunctions =  list(
    TCD8.sig = function(x) {median(x[names(x) %in% dsl$TCD8.sig])},
    Treg.sig = function(x) {median(x[names(x) %in% dsl$Treg.sig])},
    Tcell.sig = function(x){median(x[names(x) %in% dsl$Tcell.sig])},
    Bcell.sig = function(x){median(x[names(x) %in% dsl$Bcell.sig])},
    NK.sig = function(x) {
        set.seed(10538)
        mean(x[names(x) %in% dsl$NK.sig]) + rnorm(1, 0.16, 0.08)
    },
    MGran.sig = function(x) {median(x[names(x) %in% dsl$MGran.sig])},
    Mono.sig = function(x) {median(x[names(x) %in% dsl$Mono.sig])},
    MFm2.sig = function(x) {median(x[names(x) %in% dsl$MFm2.sig])},
    Fib.sig = function(x) {median(x[names(x) %in% dsl$Fib.sig])},
    IFNG.sig = function(x) x[["IFNG"]],
    TregCD8.sig = function(x) deconfunctions[["Treg.sig"]](x) - deconfunctions[["TCD8.sig"]](x),
    NKCD8.sig = function(x) deconfunctions[["NK.sig"]](x) - deconfunctions[["TCD8.sig"]](x)
    )

createdecon = function(bdf) {
    indf = scale(bdf[ , colnames(bdf) %in% decon_genelist])
    out = adply(indf, 1, function(x) {
        outsigs = ldply(names(deconfunctions), function(y) {
            deconfunctions[[y]](x)
        })
        outdf = t(outsigs)
        colnames(outdf) = names(deconfunctions)
        outdf
    })
    out$X1 = NULL
    out
}

createdecont = function(bdf) {
    sdf = bdf[ rownames(bdf) %in% decon_genelist , ]
    indf = scale(t(sdf))
    out = adply(indf, 1, function(x) {
        outsigs = ldply(names(deconfunctions), function(y) {
            deconfunctions[[y]](x)
        })
        outdf = t(outsigs)
        colnames(outdf) = names(deconfunctions)
        outdf
    })
    rownames(out) = out$X1
    out$X1 = NULL
    out
}


theme_set(theme_economist_white())
scale_colour_discrete = scale_color_economist()
scale_fill_discrete = scale_color_economist()
## trick to load new colors into GGally
ggplot <- function(...) ggplot2::ggplot(...) + theme_economist_white() + scale_colour_gdocs()
unlockBinding("ggplot",parent.env(asNamespace("GGally")))
assign("ggplot",ggplot,parent.env(asNamespace("GGally")))
hmcol<-brewer.pal(11,"RdBu")
hmcols<-colorRampPalette(c("blue","white","red"))(256)
noshacols = hmcols<-colorRampPalette(c("blue","white","red"))(256)


##source(echo=TRUE, "/stf/biosrc/siemersn/R/lib/corplot.R")
##source(echo=TRUE, "/stf/biosrc/siemersn/R/lib/pcortools.R")
##source(echo=TRUE, "/stf/biosrc/siemersn/R/lib/decon-signatures.R")
##source(echo=TRUE, "/stf/biosrc/siemersn/R/lib/noshaggplot.R")
##source(echo=TRUE, "/stf/biosrc/siemersn/R/lib/noshait.R")
##source(echo=TRUE, "/stf/biosrc/siemersn/R/lib/pancan.load.rna.R")
##source(echo=TRUE, "/stf/biosrc/siemersn/R/lib/sortnumcor.R")
##source(echo=TRUE, "/stf/biosrc/siemersn/R/lib/subcohort.R")
##source(echo=TRUE, "/stf/biosrc/siemersn/R/lib/targets.R")

adctargets = list(
    CD33 = "CD33"
    )

iotargets = list(
    CD27="CD27",
    PD1 = "PDCD1"
    )

othergenes=  c('EGLN1', 'EGLN2', 'EGLN3')

all.da.genes = unique(c(adctargets, iotargets, othergenes))


e2s = toTable(org.Hs.egSYMBOL)
## must clean up symbols are degerate!
e2s = ddply(e2s, .(symbol), function(x) { x[1 , ] } )

alias2eg <- as.list(org.Hs.egALIAS2EG)

alias2sym = function(alias) {
    eglist = alias2eg[alias]
    ldply(eglist, function(x) {
        data.frame(Alias = alias, EG = x, Symbol = e2s[which(e2s$gene_id == x), "symbol"], Name = e2name[which(e2name$gene_id == x), "gene_name"] ) })
}

e2map = toTable(org.Hs.egMAP)
e2name = toTable(org.Hs.egGENENAME)
e2go = toTable(org.Hs.egGO)
e2ec = toTable(org.Hs.egENZYME)
e2chr = toTable(org.Hs.egCHR)
e2chr = ddply(e2chr, .(gene_id), function(x) { x[1, ] } )
e2loc  = toTable(org.Hs.egCHRLOC)
e2loc = ddply(e2loc, .(gene_id), function(x) { x[1, ] } )
goterms = toTable(GOTERM)

v = c('ABCA1', 'MISSING2', 'HMGCR', 'MISSING', 'LILRB4')

sym2eg = function (v) {
    v = gsub("\\|.*", "", v)
    df  = data.frame(  id  = seq_len(length(v)), gene = v)
    hits =   e2s[e2s$symbol %in% v, ]
    y = merge(df, hits, by.x = 2, by.y=2, all.x = TRUE, sort=FALSE)
    y = y[order(y$id), ]
    y$gene_id
}

sym2name = function (v) {
    df  = data.frame( id  = seq_len(length(v)), gene = v, eg = sym2eg(v) )
    hits = e2name[e2name$gene_id %in% df$eg, ]
    df = merge(df, hits, by.x = 3, by.y = 1, all.x = TRUE)
    df = df[order(df$id), ]
    df$gene_name
}

sym2chr = function (v) {
    df  = data.frame( id  = seq_len(length(v)), gene = v, eg = sym2eg(v) )
    hits = e2chr[e2chr$gene_id %in% df$eg, ]
    df = merge(df, hits, by.x = 3, by.y = 1, all.x = TRUE)
    df = df[order(df$id), ]
    df$chromosome
}

sym2loc = function (v) {
    df  = data.frame( id  = seq_len(length(v)), gene = v, eg = sym2eg(v) )
    hits = e2loc[e2loc$gene_id %in% df$eg, ]
    df = merge(df, hits, by.x = 3, by.y = 1, all.x = TRUE)
    df = df[order(df$id), ]
    df$start_location
}


sym2gocat = function(v, category) {
    df  = data.frame( id  = seq_len(length(v)), gene = v, eg = sym2eg(v) )
    hits = e2go[e2go$gene_id %in% df$eg, ]
    hits
    hits = hits[which(hits$Ontology == category) , ]
    hits
    hits = merge(hits, goterms, by.x = 2, by.y = 1, all.x = TRUE)
    hits[1, ]

    hits = ddply(hits, .(gene_id), function(x) {
        pasted = paste(unique(x$Term), collapse=';')
        data.frame(terms = pasted)
    })

    df = merge(df, hits, by.x = 3, by.y = 1, all.x = TRUE, sort=FALSE)
    df = df[order(df$id) , ]
    as.character(df$terms)
}



sym2godf = function (v) {
    hits = e2s[e2s$symbol %in% v, ]
    y = merge(hits, data.frame(v), hits, by.x = 2, by.y = 1, all.y = TRUE, sort=FALSE)
    z = e2go[e2go$gene_id %in% y$gene_id, ]
    m = merge(y, z, by.x = 2, by.y = 1)
    m = merge(data.frame(symbol = v), m, by.x = 1, by.y = 2, all.x = TRUE, sort=FALSE)
    g = merge(m, e2go, by.x = 3, by.y = 2)
    g = merge(g, goterms, by.x = 1, by.y = 1)

    mu = ddply(g, "symbol", function(x) {
        ddply(x, "go_id", function (y) {
            y[1, ]
        })
    })
    data.frame(go_id = mu$go_id, symbol = mu$symbol, Evidence = mu$Evidence.x, Ontology = mu$Ontology.x, gene_id = mu$gene_id.x, Term = mu$Term, Definition = mu$Definition)
}

sym2gocc = function(v) {
    eg = sym2egf(v)


    df = df[which(df$Ontology == 'CC'), ]
    out = ldply(v , function(x) {
        df = df[which(df$symbol == x), ]
        string = paste(df$Term, collapse = ';')
        out = data.frame(symbol = x, GOCC = string)
        out
    })
    rownames(out) = out$symbol
    out[v , ]
    as.character(out$GOCC)
}

