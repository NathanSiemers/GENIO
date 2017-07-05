
################################################################
## all major computations happen in this block
## the other areas generally make figures., etc
################################################################



## list of meta-mutations (like activating p53 mutations). Signatures will be created for them.
mut.expressions = c(
    KRAS.G12.mmut = 'KRAS\\.p\\.G12',
    TP53.act.mmut = 'TP53\\.p\\.R175|TP53\\.p\\.R248|TP53\\.p\\.R249|TP53\\.p\\.R273'
    )


## convert and expression (above) to a 1/0 signature
meta_mutation = function(df, expression ) {
    ind = grepl(expression, colnames(df) )
    out = ifelse(Reduce('|', df[, ind, drop=FALSE]), 1, 0)
    print(length(out))
    out
}

paper.extra.genes = c(
    'CD274',
    'PTPRC',
    'LAG3',
    'PDCD1'     )


## handle all mutational signatures
meta_mutations = function (df, expressions = mut.expressions) {
    mm = ldply(expressions, function(x) { print(x); meta_mutation(df, x) } )
    print(dim(mm))
    rownames(mm) = mm$.id
    mm$.id = NULL
    t(mm)
}

## this is not working at the moment
position_mutations = function ( df ) {  ## df is final genes in columns format of fmut...
    odf = df[ , grepl("\\.p\\.", colnames(df) ) ]
    odf = t(odf)
    pmut =  as.character(gsub("[a-zA-Z]+$", "", rownames(odf), perl=TRUE ))
    odf$pmut = as.character(odf$pmut)
    fdf = ddply(odf, .(pmut),  function(x){ x$pmut = NULL ; ifelse(Reduce('|', x, 1, 0)) })
    rownames(odf) = paste(odf$pmut, ".pmut") ; odf$pmut = NULL
    t(odf)
}


## plotting routines
myplot =  function(x, y, cohort = FALSE, mydata=rnacnv, facet = FALSE, o = o.g.sig) {
    cnv.loss = gsub("\\.cnv", ".l", x)
    cnv.gain = gsub("\\.cnv", ".g", x)
    if(! cohort == FALSE) { mydata = mydata[mydata$Cohort == cohort, ] }
    if(cohort == FALSE) { cohort = 'All TCGA' }
    mydata = mydata[, colnames(mydata) %in% c(x, y, "Cohort")]
    mydata.loss = mydata[ which(mydata[ , x] < cna.cutoff), ]
    mydata.gain = mydata[ which(mydata[ , x] > -1 * cna.cutoff), ]
    ## just labels
    l.x = gsub("\\.sig", " Estimate", x)
    l.x = gsub("\\.mut", " Mutation", l.x)
    l.x = gsub("\\.cnv", " Copy Number (log2)", l.x)
    l.y = gsub("\\.sig", " Estimate (S.D. units)", y)
    l.y = gsub("\\.mut", " Mutation", l.y)
    l.y = gsub("\\.cnv", " Copy Number (log2)", l.y)
    my.P.loss = o[ o$Cohort == cohort & o$Variant == cnv.loss & o$Sig == y, "P"]
    my.P.gain = o[ o$Cohort == cohort & o$Variant == cnv.gain & o$Sig == y, "P"]
    mytitle = gsub("_", " ", cohort)
    mytitle = paste(mytitle, ifelse(! is.na(my.P.loss),
        paste("-logP =", my.P.gain, "(gain),", my.P.loss, "(loss)"),
        "") )
    my.p = ggplot(data=mydata, aes_string(x=x, y=y) ) +
        ggtitle( mytitle   ) +
            theme_economist_white() +
                theme( plot.background = element_rect(fill = 'white'),
                      panel.grid.minor=element_blank(),
                      panel.grid.major=element_blank(),
                      text = element_text(size=14)
                      ) +
                          xlab(l.x) + ylab(l.y) +
                              theme(plot.title=element_text(size=12, face='italic')) +
                                  geom_point() +
                                      geom_smooth(method='lm', data= mydata.loss, color='darkslateblue') +
                                          geom_smooth(method='lm', data= mydata.gain, color="darkslateblue")
    my.p
}



mybox =  function( x, y, cohort, mydata=rnamut, o = o.g.sig ) {
    mydata = mydata[which(mydata$Cohort == cohort), ]
    mydata = mydata[ which(! mydata[ , x] == 0.5), ]
    ##mydata = mydata[ which(! mydata[ , x] == 'NA'), ]
    ##mydata = mydata[ which(! is.na(mydata[ , x])), ]
    mydata[ , x] = gsub(0, "WT", mydata[,x])
    mydata[ , x] = gsub(1, "Mut", mydata[,x])
    mydata[ , x ] = factor(mydata[ , x], levels = c('WT', 'Mut') )
    mydata = droplevels(mydata)
    ##mydata = cbind(mydata, createdecon(mydata))
    l.x = gsub("\\.sig", " Estimate (S.D. units)", x)
    l.x = gsub("\\.mut", " Mutation", l.x)
    l.x = gsub("\\.cnv", " Copy Number (log2)", l.x)
    l.y = gsub("\\.sig", " Estimate (S.D. units)", y)
    l.y = gsub("\\.mut", " Mutation", l.y)
    l.y = gsub("\\.cnv", " Copy Number (log2)", l.y)
    ## my.p = melt.o.p[
    ##     melt.o.p$Variant == x &
    ##         melt.o.p$variable == cohort &
    ##             melt.o.p$Sig == y, "value"]
     my.p = o[
         o$Variant == x &
             o$Cohort == cohort &
                 o$Sig == y, "P"]

    mytitle = gsub("_", " ", cohort)
    mytitle = paste(mytitle, "    -logP =", my.p)
    my.p = ggplot(data=mydata, aes_string(x=x, y=y) ) +
        ggtitle( mytitle) +
            theme_economist_white() +
                theme( plot.background = element_rect(fill = 'white'),
                      panel.grid.minor=element_blank(),
                      panel.grid.major=element_blank(),
                      text = element_text(size=12)
                      ) +
                          xlab(l.x) + ylab(l.y) +
                              theme(plot.title=element_text(size=12, face='italic')) +
                                  geom_boxplot(outlier.shape = NA) + geom_jitter(position = position_jitter(width = .1))
    my.p
}


## some table making functions

kablelist = function(lis, digits=5) {
    df = data.frame(Value = lis)
    df
}
kablelol = function(lol, digits=5) {
    ldply(lol, function(lis) {
        df = data.frame(Value = lis)
        df
    })
}


## remove pipes and dashes from gene names - might want to reverse this for reports...
cleangenes = function ( x ) {
    x = gsub('\\|', '.pi.', x)
    x = gsub('\\-', '.da.', x)
    x = gsub('\\*', '.st.', x)
    x = gsub('>', '.gt.', x)
    x = gsub('<', '.lt.', x)
    x
}
## add back pipes and dashes...
uncleangenes = function ( x ) {
    x = gsub('\\.pi\\.' , '|',  x)
    x = gsub('\\.da\\.', '-',  x)
    x = gsub('\\.st\\.', '*', x)
    x = gsub('\\.gt\\.', '>', x)
    x = gsub('\\.lt\\.', '<', x)
    x
}

## get rid of everything after and including a pipe foo|ENSTxxxxxxx...  becomes foo
nopipe = function ( x ) {
    x = gsub("\\|.*", "", x)
}
nodot = function ( x ) {
    x = gsub("\\..*", "", x)
}

## blck: simple head-like function for data frames
blck = function ( x , n=5) {
    x[1:n, 1:n]
}

## filter cohorts by minimum cohort size and also remove some heme cancers, pilots
filter.cohorts = function ( df ) {
    cstats = ddply(df, .(Cohort), .parallel = TRUE, function(x){data.frame(N=length(x$Sample_Type))})
    good.cohorts = cstats[which(cstats$N >= min.samples) , "Cohort", drop=FALSE]
    good.cohorts = good.cohorts[which(! grepl(paste(bad.cohorts, collapse='|'), good.cohorts$Cohort) ), "Cohort"]
    df = df[which(df$Cohort %in% good.cohorts) , ]
    df = droplevels(df)
    df
}

## Report:  details on data files and imported data frames.
d.report = function(  df, file = NA ) {
    out.prefix = paste(filename, file, Sys.Date(), sep=".")
    if ( ! is.na(file) ) {
        my.ll = (system( paste( "ls ", paste0(data.dir,file)), intern=TRUE))
        write(my.ll, file=paste(out.prefix, 'll', 'txt', sep="."))
        my.md5 = system( paste( "md5", file), intern=TRUE)
        write(my.md5, file=paste(out.prefix, 'md5', 'txt', sep="."))
    }
    out = ddply(df, .(Cohort, Sample_Type), .parallel = TRUE, function(x) {
        length(x[, 1])
    })
    print(out)
    write.csv(out, file=paste0( data.dir, paste(out.prefix, "cohortstats", "txt", sep=".")))
    out
}

## URLs and downloaded file name of data files from UCSC
##  mutations
muturl = 'https://genome-cancer.ucsc.edu/download/public/xena/TCGA/TCGA.PANCAN.sampleMap/PANCAN_mutation_gene'
mutfile = paste0(data.dir, "PANCAN_mutation_gene.txt")
mutmeta = "https://genome-cancer.ucsc.edu/download/public/xena/TCGA/TCGA.PANCAN.sampleMap/PANCAN_mutation_gene.json"

## clinical
clinurl = 'https://genome-cancer.ucsc.edu/download/public/xena/TCGA/TCGA.PANCAN.sampleMap/PANCAN_clinicalMatrix'
clinfile = paste0(data.dir, "PANCAN_clinicalMatrix.txt")
clinmeta = "https://genome-cancer.ucsc.edu/download/public/xena/TCGA/TCGA.PANCAN.sampleMap/PANCAN_clinicalMatrix.json"

## cnv
cnvurl = 'https://genome-cancer.ucsc.edu/download/public/xena/TCGA/TCGA.PANCAN.sampleMap/Gistic2_CopyNumber_Gistic2_all_data_by_genes'
cnvfile = paste0(data.dir, "Gistic2_CopyNumber_Gistic2_all_data_by_genes.txt")
cnvmeta = "https://genome-cancer.ucsc.edu/download/public/xena/TCGA/TCGA.PANCAN.sampleMap/Gistic2_CopyNumber_Gistic2_all_data_by_genes.json"

## rna
rnaurl = 'https://genome-cancer.ucsc.edu/download/public/xena/TCGA/TCGA.PANCAN.sampleMap/HiSeqV2'
rnafile = paste0(data.dir, "HiSeqV2.txt" )
rnameta = "https://genome-cancer.ucsc.edu/download/public/xena/TCGA/TCGA.PANCAN.sampleMap/HiSeqV2.json"

## mutations at specific positions
fmuturl = "https://genome-cancer.ucsc.edu/download/public/xena/TCGA/TCGA.PANCAN.sampleMap/PANCAN_mutation"
fmutfile = paste0(data.dir, "PANCAN_mutation")
fmutmeta = "https://genome-cancer.ucsc.edu/download/public/xena/TCGA/TCGA.PANCAN.sampleMap/PANCAN_mutation.json"


colonclin = 'https://genome-cancer.ucsc.edu/download/public/xena/TCGA/TCGA.COADREAD.sampleMap/COADREAD_clinicalMatrix'
colonclinmeta = 'https://genome-cancer.ucsc.edu/download/public/xena/TCGA/TCGA.COADREAD.sampleMap/COADREAD_clinicalMatrix.json'
colon.scid = 'CDE_ID_3226963'

breastclin = 'https://genome-cancer.ucsc.edu/download/public/xena/TCGA/TCGA.BRCA.sampleMap/BRCA_clinicalMatrix'
breastclinmeta = 'https://genome-cancer.ucsc.edu/download/public/xena/TCGA/TCGA.BRCA.sampleMap/BRCA_clinicalMatrix.json'
breast.scid = 'PAM50Call_RNAseq'

stomachclin = 'STADMasterPatientTable20140207.csv'
stomach.scid = 'Molecular.Subtype'

## download data from UCSC
## this is not sufficient to stop errors if there's a netowrk problem because of str() of objects.
rna.json = cnv.json = clin.json = mut.json = fmut.json = ""

if (do.download) {
    download.file(muturl, mutfile, cacheOK = FALSE, method='curl')
    download.file(method='curl', cacheOK = FALSE, cnvurl,  cnvfile)
    download.file(method='curl', cacheOK = FALSE, rnaurl,  rnafile)
    download.file(method='curl', cacheOK = FALSE, clinurl,  clinfile)
    download.file(method='curl', cacheOK = FALSE, fmuturl,  fmutfile)
    download.file(method='curl', cacheOK = FALSE, mutmeta,  paste0(data.dir, "mut.meta"))
    download.file(method='curl', cacheOK = FALSE, cnvmeta,  paste0(data.dir,"cnv.meta"))
    download.file(method='curl', cacheOK = FALSE, rnameta, paste0(data.dir, "rna.meta"))
    download.file(method='curl', cacheOK = FALSE, clinmeta,  paste0(data.dir,"clin.meta"))
    download.file(method='curl', cacheOK = FALSE, fmutmeta,  paste0(data.dir,"fmut.meta"))
    download.file(method='curl', cacheOK = FALSE, colonclin,  paste0(data.dir,"colon.clin"))
    download.file(method='curl', cacheOK = FALSE, colonclinmeta,  paste0(data.dir,"colon.clin.meta"))
    download.file(method='curl', cacheOK = FALSE, breastclin,  paste0(data.dir,"breast.clin"))
    download.file(method='curl', cacheOK = FALSE, breastclinmeta,  paste0(data.dir,"breast.clin.meta"))
}


## load TCGA meta data
colon.json = fromJSON(paste0(data.dir,"colon.clin.meta"))
breast.json = fromJSON(paste0(data.dir,"breast.clin.meta"))
rna.json = fromJSON( paste0(data.dir,"rna.meta"))
cnv.json = fromJSON( paste0(data.dir,"cnv.meta"))
clin.json = fromJSON( paste0(data.dir,"clin.meta"))
mut.json = fromJSON( paste0(data.dir,"mut.meta"))
fmut.json = fromJSON( paste0(data.dir,"fmut.meta"))

################################################################
## compile a subcohort classification
cc = read.table(paste0(data.dir,"colon.clin"), sep="\t", header=TRUE, comment.char = '')
cc$sampleID = gsub('-', '.', cc$sampleID)
ccdf = data.frame(sampleID = cc$sampleID, subcohort = cc[ , colon.scid] )
ccdf$subcohort = paste0('COAD_', ccdf$subcohort)


bc = read.table(paste0(data.dir,"breast.clin"), sep="\t", header=TRUE, comment.char = '')
bc$sampleID = gsub('-', '.', bc$sampleID)
bcdf = data.frame(sampleID = bc$sampleID, subcohort = bc[ , breast.scid] )
bcdf$subcohort = paste0('BRCA_', bcdf$subcohort)


sc = read.table(paste0(data.dir,stomachclin), sep=",", header=TRUE, comment.char = '')
sc$sampleID = gsub('-', '.', sc[ , 'TCGA.barcode'] )
sc$sampleID = paste0(sc$sampleID, '.01' )
scdf = data.frame(sampleID = sc$sampleID, subcohort = sc[ , stomach.scid] )
scdf$subcohort = paste0('STAD_', scdf$subcohort)

subcohorts = rbind(ccdf, bcdf, scdf)
dim(subcohorts)

################################################################
## read iogenes io genes list
iogenes = read.csv("iogenes.csv", header = FALSE, stringsAsFactors = FALSE)
iogenes = iogenes$V1
iogenes = iogenes[! grepl('\\.sig', iogenes)]
iogenes = iogenes[! grepl('\\.dif', iogenes)]


################################################################
## load USCS TCGA data files into data frames


## Mutations:  specific position mutations, all 1.5 million of them.
## must work the data into a similar format as gene-level mutations
## without causing regional power outages

if(! skip.early) {
    if( ! use.cached ) {
        if ( force.reload |  ! exists("myfmut") ) {

            myfmut.o = read.table(fmutfile, sep="\t", quote="", header=TRUE, as.is=TRUE, stringsAsFactors=FALSE, comment.char = '')
                                        # ugliness, clean it up now, not later!
            myfmut.o$Amino_Acid_Change = gsub('^\\.$', 'NULL', myfmut.o$Amino_Acid_Change)
            myfmut.o$Amino_Acid_Change = gsub('^$', 'NULL', myfmut.o$Amino_Acid_Change)
            myfmut.o$sample = gsub('\\-', '.', myfmut.o$sample)
            save(myfmut.o, file=paste0(data.dir,"myfmut.o.Robject"))
            ## must filter! 1.5M different mutations!  also try to label AA change consistently (they aren't)

            myfmut.oo = ddply(myfmut.o, .(gene), .parallel = TRUE, function (y) {
                ##print(y)
                out =ddply( y,  .( start, alt), .parallel = FALSE, failwith(NULL, function(x) {
                    ##print(x)
                    aa = x[ , "Amino_Acid_Change", drop=FALSE]
                    my.t = table(aa)
                    my.t = my.t[names(my.t) != 'NULL']
                    if (length(my.t) > 0 ) {
                        my.best = names(sort(my.t[names(my.t) != 'NULL'], decreasing = TRUE)[1])
                        x$Amino_Acid_Change = paste0(x$gene, '.', my.best)
                    } else {
                        x$Amino_Acid_Change = paste0(x$gene, '.', x$start, x$reference, x$alt)
                    }
                    x
                }) )
                out
            })

            head(myfmut.oo, 200)

            totalfmut = ddply(myfmut.oo, .(sample), .parallel = FALSE, summarize, totalmut = log2(length(sample)) )
            rownames(totalfmut) = totalfmut$sample; totalfmut$sample = NULL
            ## don't add totalfmut to df until after cast...

                                        # TRIM
            each.mut = table(myfmut.oo$Amino_Acid_Change)

            ##ddply(myfmut.oo, .(Amino_Acid_Change), .parallel = FALSE, summarize, totalmut = length(Amino_Acid_Change) )
            each.name = names(each.mut[ each.mut   > fmut.min  ])

            myfmut.ooo = myfmut.oo[myfmut.oo$Amino_Acid_Change %in% each.name , ]

            library(data.table)
            setDT(myfmut.ooo)
            myfmut = dcast.data.table(myfmut.ooo, sample ~ Amino_Acid_Change, fun=length)
            myfmut = as.data.frame(myfmut)

            ## the big cast
            ##myfmut = dcast(myfmut.ooo, sample ~ Amino_Acid_Change, fun.aggregate = length )

            rownames(myfmut) = myfmut$sample; myfmut$sample = NULL
            myfmut[myfmut > 1] = 1;
            colnames(myfmut) = cleangenes(colnames(myfmut))

            ## add totalmut column
            myfmut = merge(myfmut, totalfmut, by=0); rownames(myfmut) = myfmut$Row.names; myfmut$Row.names = NULL
            ## add composite mutations
            myfmut = cbind(myfmut, meta_mutations(myfmut) )

            save(myfmut, file = paste0(data.dir, paste(filename, "myfmut", "Robject", sep=".")))
            ## send to global namespace so knitit() will not have to generate again...
            assign('myfmut', myfmut, envir = .GlobalEnv)

            ## mypmut.oo = ddply( myfmut.o,  .(gene), .parallel = TRUE, function(y) {
            ##   out = ddply(y, .(start), .parallel = FALSE, failwith(NULL, function(x) {
            ##     ##if( length( x$gene )  < fmut.min ) {return(NULL)}
            ##     ## do the above later!!!!
            ##     my.skip = FALSE
            ##     aa = x$Amino_Acid_Change
            ##     aavec = aa[ ! ( is.null(aa)  | aa == 'NA' | aa == '.'  | aa == '' | aa == 'NULL' )  ]
            ##     if(length(aavec) == 0) {aavec = c('NULL' = 'NULL') }
            ##     aavec = as.vector(aa)
            ##     aachange = names(sort(table(aavec) ,decreasing=TRUE, na.last = TRUE)[1])
            ##     if(length(aachange) == 0){print(x); my.skip = TRUE}
            ##     ##cat(aachange); cat(' ')
            ##     mut = ifelse(
            ##         (my.skip | aachange == NA | is.null(aachange) |
            ##              aachange == 'NULL' | aachange == '' | aa == '.' ) ,
            ##         paste0(x$gene[[1]], '.', x$start[[1]], x$reference[[1]], x$alt[[1]] ),
            ##         paste0(x$gene[[1]], '.', aachange) )
            ##     ##mut = gsub('\\-', '.d.', mut); mut = gsub('\\*', '.s.', mut)
            ##     if(my.skip) { print ( mut ) }
            ##     ## different!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ##     ##cat(mut); cat ( ' ' )
            ##     ##print(x)
            ##     x
            ## }) )
            ##   out  })

            ## mypmut.oo$Amino_Acid_Change = gsub('[A-Z]+$', '', mypmut.oo$Amino_Acid_Change)
            ## mypmut = dcast(mypmut.oo, sample ~ Amino_Acid_Change, fun.aggregate = length )
            ## ## remove some dups, move sample column to row name
            ## rownames(mypmut) = mypmut$sample; mypmut$sample = NULL
            ## mypmut[mypmut > 1] = 1;
            ## colnames(mypmut) = cleangenes(colnames(mypmut))  ## perhaps redundant, see above
            ## ## add totalmut column
            ## mypmut = merge(mypmut, totalfmut, by=0); rownames(mypmut) = mypmut$Row.names; mypmut$Row.names = NULL
            ## save(mypmut, file = paste(filename, "mypmut", "Robject", sep="."))
            ## ## send to global namespace so knitit() will not have to generate again...
            ## assign('mypmut', mypmut, envir = .GlobalEnv)



        }
    }

    if ( use.cached | ! exists('myfmut') ) {
        load(file = paste0(data.dir,paste(filename, "myfmut", "Robject", sep=".")))
        ##load(file = paste(filename, "mypmut", "Robject", sep="."))
        load(file = paste0(data.dir,'myfmut.o.Robject'))
    }
    dim(myfmut)

    allfmut = colnames(myfmut)


    ## Mutations: 0, 1, and in a few cases 0.5 (transposed)
    ## ucsc xena mutations, should be about 7000 samples


    if( ! use.cached ) {
        if ( force.reload |  ! exists("mymut") ) {
            mymut = read.table(mutfile, sep="\t", quote="", header=TRUE, as.is=TRUE, stringsAsFactors=FALSE, comment.char = "")
            ## sample column is gene
            mymut$sample = cleangenes(mymut$sample)
            mymut$sample = paste(mymut$sample, ".mut", sep="")
            rownames(mymut) = mymut$sample
            mymut$sample = NULL
            ## you can't do this here, matrix is transposed!!!
            ##mymut = merge(mymut, totalfmut, by=0); rownames(mymut) = mymut$Row.names; mymut$Row.names = NULL
            save(mymut, file = paste0(data.dir, paste(filename, "mymut", "Robject", sep=".")))
            ## send to global namespace so knitit() will not have to generate again...
            assign('mymut', mymut, envir = .GlobalEnv)
        }
    }


    if ( use.cached | ! exists('mymut') ) {
        load(file = paste0(data.dir,paste(filename, "mymut", "Robject", sep=".")))
    }
    dim(mymut)





    ## not used...
    ## older june pancan release for mutation
    ##load(file='new.pancan.load.rnamut.Robject')

    ## load GISTIC CNV calls, continuous (transposed)

    if( ! use.cached ) {
        if ( force.reload |  ! exists("mycnv") ) {
            mycnv = read.table(cnvfile, sep="\t", quote="", header=TRUE, as.is=TRUE, stringsAsFactors=FALSE, comment.char = '')
            mycnv$Sample = cleangenes(mycnv$Sample)
            mycnv$Sample = paste(mycnv$Sample, ".cnv", sep="")
            rownames(mycnv) = mycnv$Sample
            mycnv$Sample=NULL
            dim(mycnv)
            blck(mycnv)
            save(mycnv, file = paste0(data.dir,paste(filename, "mycnv", "Robject", sep=".")))
            assign('mycnv', mycnv, envir = .GlobalEnv)
        }
    }

    if ( use.cached  | ! exists('mycnv')) {
        load(file = paste0(data.dir,paste(filename, "mycnv", "Robject", sep=".")))
    }
    dim(mycnv)


    ## load RNA-seq, log2 (transposed)


    if( ! use.cached ) {
        if ( force.reload |  ! exists("myrna") ) {
            myrna = read.table(rnafile, sep="\t", quote="", header=TRUE, as.is=TRUE, stringsAsFactors=FALSE, comment.char = '')
            myrna$Sample = cleangenes(myrna$Sample)
            rownames(myrna) = myrna$Sample
            myrna$Sample = NULL
            myrna = as.data.frame(t(myrna))
            dim(myrna)
            blck(myrna)
            ## load Butte data set.
            butte = read.csv(paste0(data.dir,"butte-estimates-ncomms9971-s2.csv"), stringsAsFactors = FALSE, header = TRUE)
            butte$Sample = gsub('\\-', '.', butte$Sample.ID)
            butte$Sample = gsub("[A-Z]$", '', butte$Sample)
            butte = butte[ ! is.nan(butte$CPE), ]
            butte$Estimate = log2(100 * butte$CPE)
            myrnab = merge(myrna, butte[ , 9:10], by.x = 0, by.y='Sample')
            ## should we just calibrate all CNVs vs butte predictor?
            ##mycnvb = merge(as.data.frame(t(mycnv)), butte[ , 9:10], by.x = 0, by.y = 'Sample')
            mycnvb = merge(butte[ , 9:10], as.data.frame(t(mycnv)), by.y = 0, by.x = 'Sample')
            mycnvb$Sample = NULL; dim(mycnvb)

            bmelt = melt(mycnvb, id.vars = 'Estimate')
            b.stdev = sd(bmelt$value)
            b.mean = mean(bmelt$value)
            ## trim to reasonable purity levels > 5
            bmelt = bmelt[bmelt$Estimate > 5, ]

            bsample = bmelt[sample(nrow(bmelt), 8000000), ];
            ## now for weirdness, look at sd of bins

            n = 25
            bsample$bin = cut(bsample$Estimate, n)    ##, labels = 1:n)
            bins = ddply(bsample, .(bin), function(x) { sd(x$value) } )
            ddply(bsample, .(bin), count)
            ## below is a plot of the effect of purity (binned) on variance of CNV estimates
            ## qplot(bin, V1, data = bins)
            bin.range = max(bsample$Estimate) - min(bsample$Estimate)
            bin.dim = (1:n * (bin.range/n) ) + min(bsample$Estimate)
            ## qplot(bin.dim, bins$V1) + geom_smooth(method = 'lm', formula = y ~  x )
            lmf  = lm (bins$V1 ~ bin.dim )
            summary(lmf)


            pmtrim = 2 * b.stdev
            bmeltminus = bmelt[bmelt$value < -1 * pmtrim, ]
            bmeltplus = bmelt[bmelt$value > pmtrim, ]


            ## qplot(Estimate, value, data = bsample) + geom_smooth()
            plus.sampled = bmeltplus[sample(nrow(bmeltplus), 20000), ]
            lmf = lm ( value ~ 0 + Estimate, data = plus.sampled )
            plus.coeffs = lmf$coefficients
            summary(lmf)
            ## qplot(Estimate, value, data = plus.sampled) + geom_smooth(method = 'lm', formula = y ~ 0 + x ) + geom_smooth()
            minus.sampled = bmeltminus[sample(nrow(bmeltminus), 20000), ]
            lmf = lm ( value ~ Estimate, data = minus.sampled )
            summary(lmf)
            minus.coeffs = lmf$coefficients
            ## qplot(Estimate, value, data = minus.sampled) + geom_smooth(method = 'lm')


### Trials to find an alternate purity predictor.
            ind = sapply(myrnab, is.numeric)
            ## find positively correlated genes to Butte estimate
            zb = na.omit(as.data.frame(as.table(cor(myrnab$Estimate, myrnab[, ind]))))
            ## no sigs allowed
            zb = zb[! grepl("\\.sig", zb$Var2), ]
            ## no perfect matches
            zb = zb[zb$Freq < 0.99, ]
            zb = zb[order(zb$Freq, decreasing = TRUE), ]
            ## create lm formula
            zb.formula = paste(zb[1:50, "Var2"], collapse = ' + ' )
            ##zb.formula = paste( " Estimate ~  0 + Cohort + " , zb.formula )
            zb.formula = paste( " Estimate ~ " , zb.formula )
            ## run lm
            lmfit = lm(zb.formula, myrnab)
            ## apply predictor to complete myrna data set
            myrna$PEst = predict(lmfit, myrna)
            myrna = cbind(myrna, createdecon(myrna))
            ## clean up
            rm(myrnab); gc()
            save(myrna, file = paste0(data.dir,paste(filename, "myrna", "Robject", sep=".")))
            assign('myrna', myrna, envir = .GlobalEnv)
        }
    }




    if ( use.cached  | ! exists('myrna') ) {
        load(file = paste0(data.dir,paste(filename, "myrna", "Robject", sep=".")))
    }
    dim(myrna)



    ## load Clinical information

    if( ! use.cached ) {
        if ( force.reload |  ! exists("myclin") ) {
            c = read.table(clinfile, sep="\t", header=TRUE, comment.char = '')
            c$sampleID = gsub('-', '.', c$sampleID)
            ## clean up cohort names
            c$X_cohort = gsub("^TCGA ", "", c$X_cohort)
            c$X_cohort = gsub(" \\& ", "_", c$X_cohort)
            c$X_cohort = gsub(" ", "_", c$X_cohort)
            ## remove bad cohorts
            c = c[! c$X_cohort == "", ]  ## blank cohort name
            c = c[! grepl("ilot", c$X_cohort), ]  ## paraffin pilot cohort
            myclin = data.frame(
                Sample = c$sampleID,
                Sample_Type = factor(c$sample_type),
                Event = factor(c$X_EVENT, levels = c(0, 1)),
                OS = as.numeric(c$X_OS),
                Anatomical_Origin = factor(c$X_anatomical_origin),
                Cohort = factor(c$X_cohort),
                Gender = factor(c$gender),
                Patient = c$X_PATIENT )
            rownames(myclin) = myclin$Sample
            myclin$sampleID = NULL
            ## remove samples with no sample_type annotation
            print(dim(myclin))
            myclin = myclin[ ! myclin$Sample_Type == "" , ]
            ## new - add subcohort classifications
            myclin = merge(myclin, subcohorts, by.x = 'Sample', by.y = 'sampleID', all.x = TRUE)
            print(dim(myclin))
            dim(myclin)
            blck(myclin)
            colnames(myclin)
            rownames(myclin) = myclin$Sample
            save(myclin, file = paste0(data.dir,paste(filename, "myclin", "Robject", sep=".")))
            assign('myclin', myclin, envir = .GlobalEnv)
        }
    }



    if ( use.cached  | ! exists('myclin')) {
        load(file = paste0(data.dir,paste(filename, "myclin", "Robject", sep=".")))
    }
    dim(myclin)
    ## create tables of stats for the cohorts?
    if(create.dreports){
        d.report(myclin)
        ## report statistics on cohorts with RNA, Mutation, CNV
        d.report(merge(myclin, myrna, by="row.names"), rnafile)
        gc()
        d.report(merge(myclin, t(mymut), by="row.names"), mutfile)
        gc()
        d.report(merge(myclin, t(mycnv), by="row.names"), cnvfile)
        gc()
    }



################################################################
    ## create merged data sets
################################################################

    ## rna and clinical
    if( ! use.cached ) {
        if ( force.reload |  ! exists("rna") ) {
            ## rna
            rna = merge(myclin, myrna, by="row.names")
            rownames(rna) = rna$Row.names; rna$Row.names = NULL
            dim(rna)
################################################################
            save(rna, file = paste0(data.dir,paste(filename, "rna", "Robject", sep=".")))
            assign('rna', rna, envir = .GlobalEnv)
        }
    }

if ( use.cached  | ! exists('rna')) {
    load(file = paste0(data.dir,paste(filename, "rna", "Robject", sep=".")))
}

dim(rna)
blck(rna)

}

## rna-clin-cnv

if( ! use.cached ) {
    if ( force.reload |  ! exists("rnacnv") ) {
        rnacnv = merge(rna, t(mycnv), by="row.names")
        rownames(rnacnv) = rnacnv$Row.names; rnacnv$Row.names = NULL
        rnacnv = filter.cohorts(rnacnv)
        rnacnv = droplevels(rnacnv)
        save(rnacnv, file = paste0(data.dir,paste(filename, "rnacnv", "Robject", sep=".")))
        assign('rnacnv', rnacnv, envir = .GlobalEnv)
        gc()
    }
}


if ( use.cached  | ! exists('rnacnv')) {
load(file = paste0(data.dir,paste(filename, "rnacnv", "Robject", sep=".")))
}

if(remove.intermediate){
gc()
}

## rna-clin-mutation

if( ! use.cached ) {
    if ( force.reload |  ! exists("rnamut") ) {
        rnamut = merge(rna, t(mymut), by="row.names")
        rownames(rnamut) = rnamut$Row.names; rnamut$Row.names = NULL
        rnamut = filter.cohorts(rnamut)
        rnamut = droplevels(rnamut)
        dim(rnamut)
        rnamut = merge(rnamut, totalfmut, by=0); rownames(rnamut) = rnamut$Row.names; rnamut$Row.names = NULL
        dim(rnamut)
        save(rnamut, file = paste0(data.dir,paste(filename, "rnamut", "Robject", sep=".")))
        assign('rnamut', rnamut, envir = .GlobalEnv)
        gc()
    }
}

if ( use.cached  | ! exists('rnamut')) {
load(file = paste0(data.dir,paste(filename, "rnamut", "Robject", sep=".")))
}
if(remove.intermediate){
rm(mymut)
gc()
}

## rna-clin-fine mutation
if( ! use.cached ) {
    if ( force.reload |  ! exists("rnafmut") ) {
        rnafmut = merge(rna, myfmut, by="row.names")
        rownames(rnafmut) = rnafmut$Row.names; rnafmut$Row.names = NULL
        rnafmut = filter.cohorts(rnafmut)
        rnafmut = droplevels(rnafmut)
        save(rnafmut, file = paste0(data.dir,paste(filename, "rnafmut", "Robject", sep=".")))
        assign('rnafmut', rnafmut, envir = .GlobalEnv)
        gc()
    }
}

if ( use.cached  | ! exists('rnfamut')) {
load(file = paste0(data.dir,paste(filename, "rnafmut", "Robject", sep=".")))
}
if(remove.intermediate){
rm(myfmut)
gc()
}



##  rna-clin-mutation-cnv merged - expensive
if( ! use.cached ) {
    if(create.allmerge) {
        mymerge = merge(rnamut, t(mycnv), by="row.names")
        rownames(mymerge) = mymerge$Row.names; mymerge$Row.names = NULL
        ##mymerge = filter.cohorts(mymerge)
        dim(mymerge)
        gc()
        save(mymerge, file = paste0(data.dir,paste(filename, "mymerge", "Robject", sep=".")))
        save(mymerge, file = paste0(data.dir,paste(filename, "mymerge", Sys.Date(), "Robject", sep=".")))
        ## Cohort summaries for different two-way genomic data sets
        cstats = ddply(mymerge, .(Cohort, Sample_Type), .parallel = TRUE, function(x){data.frame(N=length(x$Sample_Type))})
        print(kable(cstats))
        cstats = ddply(mymerge, .(Cohort, Sample_Type), .parallel = TRUE, function(x){data.frame(N=length(x$Sample_Type))})
        print(kable(cstats))
        gc()
    }
}


################################################################
## assess correlation of signature markers across cohorts
################################################################
if ( ! use.cached | ! exists('sig.cut') | force.reload ) {

    ##sig.df = rnacnv[ , c('Cohort', 'subcohort', decon_genelist) ]
    sig.df = rnacnv[ , c('Cohort', decon_genelist) ]

    sig.cor =
        ddply( sig.df, .(Cohort), function (x) {
        my.cohort = as.character(x[ 1, 1 ])
        x$Cohort = NULL
        print(head(x, 1))
        my.cor = as.data.frame(cor ( x, x ) )
        my.cor$Cohort = my.cohort
        my.cor
    })


    sig.val = ddply(sig.cor, .(Cohort), function(x) {
        ldply(decon_signature_list, function(y) {
            tmp.df = x[ , y]
            print(y)
            all.cors = unlist(cor(tmp.df, tmp.df))
            median ( all.cors[which(all.cors < 0.99999)] )
        })
    })


    colnames(sig.val) = c('Cohort', 'Sig', 'Cor' )
    sig.cut = sig.val[sig.val$Cor > minimum.sig.coherence , ]
    sig.cut$Cor = NULL
    sig.cut = sig.cut[complete.cases(sig.cut), ]
    sig.cut = droplevels(sig.cut)



################################################################
    ## ugly code to see if ratios should be included


    sig.nkcd8 = ddply(sig.cut, .(Cohort), function(y) {
        if ( length ( grep ( 'NK.sig', fixed = TRUE, y$Sig ) ) > 0 &
            length ( grep ( 'TCD8.sig', fixed = TRUE, y$Sig ) ) > 0  ) {
            out = data.frame ( Cohort = y$Cohort[[1]], Sig = 'NKCD8.sig' )
        } else {
            out = NULL
        }
        out
    } )

    sig.tregcd8 = ddply(sig.cut, ~ Cohort, function(x) {
        if ( length ( grep ( 'Treg.sig', fixed = TRUE, x$Sig ) ) > 0 &
            length ( grep ( 'TCD8.sig', fixed = TRUE, x$Sig ) ) > 0  ) {
            data.frame ( Cohort = x$Cohort[[1]], Sig = 'TregCD8.sig' )
        }
    } )

    sig.cut = rbind(sig.cut, sig.nkcd8, sig.tregcd8)

    sub.cut = ddply(sig.cut, .(Cohort, Sig), function (x) {
        df = NULL
        if(x$Cohort == 'Colon_Cancer') {
            df = data.frame(
                Cohort = c(
                    'COAD_Indeterminate',
                    'COAD_MSI.L',
                    'COAD_MSS',
                    'COAD_MSI.H'),
                Sig = x$Sig
                )
        }
        if(x$Cohort == 'Breast_Cancer') {
            df = data.frame(
                Cohort = c(
                    'BRCA_Normal',
                    'BRCA_Basal',
                    'BRCA_LumA',
                    'BRCA_LumB'),
                Sig = x$Sig
                )
        }
        if(x$Cohort == 'Stomach_Cancer') {
            df = data.frame(
                Cohort = c(
                    'STAD_MSI',
                    'STAD_CIN',
                    'STAD_EBV',
                    'STAD_GS'
                    ),
                Sig = x$Sig
                )
        }
        df
    })

    sig.cut = rbind(sig.cut, sub.cut)


}

################################################################




################################################################
## identify most frequent mutant genes and most variant CNVs
################################################################

## mutations and copy number alterations to be studied
topmutated = head( sort ( apply(rnamut[ , grepl("\\.mut", colnames(rnamut))], 2, sum), decreasing=TRUE ), topmut.cutoff)

## topcnv cutoff has been removed.

##topcnv = head( sort (apply(rnacnv[, grepl("\\.cnv", colnames(rnacnv))], 2, var), decreasing=TRUE ), topcnv.cutoff)

## topcnv.g =  head( sort (apply(   rnacnv[, grepl("\\.cnv", colnames(rnacnv))],    2, function(x){  var(x[x > -1 * cna.cutoff])   }), decreasing=TRUE ), topcnv.cutoff)


##     head( sort (apply(   rnacnv[, grepl("\\.cnv", colnames(rnacnv))],    2, function(x){  var(x[x < cna.cutoff])   }), decreasing=TRUE ), topcnv.cutoff)
##     ) )


## names(topcnv.l) = gsub("\\.cnv", ".l", names(topcnv.l))
## names(topcnv.g) = gsub("\\.cnv", ".g", names(topcnv.g))

allcnv = apply(rnacnv[, grepl("\\.cnv", colnames(rnacnv))], 2, var)
##filtered.variants = c(topmutated, topcnv.l, topcnv.g)

## clean up memory on my mac
if(remove.intermediate){
rm(myrna); gc()
}


################################################################
## Statistical Analysis
################################################################
## from Han and Nathan
####################################################################################


################################################################
## contrast making
################################################################

## make multiple contrasts, wt vs mut in each levels of type
## based on design matrix specified above in "mut_io function"
## levels in mut2 factor are in order of "WT", "MUT"
contrastType<-function(type){
n<-length(levels(type))
## base level/type
m<-c(0,1, rep(0,n-1), rep(0,n))
### addtional levels/types
for(i in 2:n){
    c1<-c(0, 1, rep(0,n-1), rep(0,i-1), 1, rep(0,n-i))
    m<-cbind(m,c1)
}
colnames(m)<-paste(rep("V", n), levels(type), sep=".")
m
}

## Nathan
## contrasts when mutation load not used as covariate (for CNVs)
contrastTypeCNV<-function(type){
n<-length(levels(type))
m<-c(0,1, rep(0,n-1), rep(0,n-1))
for(i in 2:n){
    c1<-c(0, 1, rep(0,n-1), rep(0,i-2), 1, rep(0,n-i))
    m<-cbind(m,c1)
}
colnames(m)<-paste(rep("V", n), levels(type), sep=".")
m
}

## Han
### R functions to compute associations of cancer genetics vs. IO gene/signature expression
####### mutation vs. IO gene/signature expressions
#################### inputs  ##################
####  io --- matrix: io gene/signature expression, log2 transformed
                                        # rows are samples, columns are gene/signatures.  colnames are gene/signature names
####  mut --- factor: mutation status of samples (matching rows in "io").
                                        # levels: WT, Mut (or n/no, y/yes).  make sure WT/n/no is the first level
#### type --- factor: cohort/type labels of samples (matching rows in "io")
#### ml --- co-variable,  mutation load of samples (matching rows in "io")
                                        # log2 transformed
mut_io<-function(io, mut, type, ml){
### make design matrix
design<-model.matrix(~mut*type+ml)
### make contrasts matrix
cm<-contrastType(type)
### fit lm
fit<-lmFit(t(io),design)
### contrasts fit
fit2<-contrasts.fit(fit, cm)
fit2<-eBayes(fit2)
### extract overall P value (F test of all contrasts == 0)
mytable<-topTable(fit2, number=500, sort.by="none")
mytable = cbind(Sig = rownames(mytable), mytable)
### extract contrast for Mut vs. WT in each cohort/type
n<-length(colnames(cm))
for(i in 1:n){
    ## get t, p-value, adjust.p.value for each contrast
    mt<-topTable(fit2, coef=i, number=500, sort.by="none")[,c(3:5)]
    colnames(mt)<-paste(rep(colnames(cm)[i], 3), colnames(mt), sep=".")
    ##mt<-topTable(fit2, coef=i, number=500, sort.by="none")[,c(4), drop=FALSE]
    ##colnames(mt)<-paste(rep(colnames(cm)[i], 1), colnames(mt), sep=".")
    ## combine matrix
    mytable<-cbind(mytable, mt)
}
mytable
}

## Nathan
## cnv_io - use more samples, don't include totalmut as covariate
## slightly altered contrast tables because of this
cnv_io<-function(io, cnv, type, estimate ){
### make design matrix
design<-model.matrix(~cnv*type+estimate)
### make contrasts matrix
##    cm<-contrastTypeCNV(type)
cm<-contrastType(type)
### fit lm
fit<-lmFit(t(io),design)
##print(topTable(fit))
### contrasts fit
fit2<-contrasts.fit(fit, cm)
fit2<-eBayes(fit2)
### extract overall P value (F test of all contrasts == 0)
mytable<-topTable(fit2, number=500, sort.by="none")
mytable = cbind(Sig = rownames(mytable), mytable)
##print(head(mytable))
### extract contrast for Cnv vs. WT in each cohort/type
n<-length(colnames(cm))
for(i in 1:n){
    ## get t, p-value, adjust.p.value for each contrast
    mt<-topTable(fit2, coef=i, number=500, sort.by="none")[,c(3:5)]
    colnames(mt)<-paste(rep(colnames(cm)[i], 3), colnames(mt), sep=".")
    ##mt<-topTable(fit2, coef=i, number=500, sort.by="none")[,c(4), drop=FALSE]
    ##colnames(mt)<-paste(rep(colnames(cm)[i], 1), colnames(mt), sep=".")
    ## combine matrix
    mytable<-cbind(mytable, mt)
}
mytable
}


####################################################################################
## functions to clean up cohorts
## nathan
## these routine remove cohorts whose stats are not estimable (need more than 1 level!)
################################################################

estimable_mut = function ( marker, df  ) {
##cat(paste(marker, ":"))
sdf = cbind( data.frame(
    mut = df[ , marker],
    Cohort = df$Cohort,
    totalmut = df$totalmut) ,
    df[ , grepl(".sig$", colnames(df) )],
    df[ , colnames(df) %in% iogenes ]
            )
out = ddply(sdf, .(Cohort), function(x) {
    ## print(unique(x[, "mut"]))
    if ( length(unique(x[, "mut"])) < 2 ) {return(NULL)}
    x
})
## print(dim(out))
if(dim(out)[[1]] < 2) {return(NULL)}
out
}

## you don't *really* need this for CNVs, but I do it anyway
estimable_cnv = function ( marker, df  ) {
sdf = cbind( data.frame(
    cnv = df[ , marker], Cohort = df$Cohort, Estimate = df$PEst )  ,
    df[ , grepl(".sig$", colnames(df) )] ,
    df[ , colnames(df) %in% iogenes ]
            )

out = ddply(sdf, .(Cohort), function(x) {
    if ( length(unique(x[, "cnv"])) < 2 ) {return(NULL)}
    x
})
if(dim(out)[[1]] < 2) {return(NULL)}
out
}

################################################################
## poly functions clean up data and run on one Variant
################################################################

## polymut and polycnv functions: remove non-estimable cohorts, drop levels, add decon, run stats
polymut = function(mut, df = rnamut) {
sdf = estimable_mut(mut, df = df)
if(is.null(sdf)){return(NULL)}
sdf = droplevels(sdf[complete.cases(sdf), ])
if(length(unique(sdf$Cohort)) < 2) {return(NULL)} # warning, this will skip mutations in a single cohort!
decon = cbind ( sdf[, grepl(".sig$", colnames(sdf))] , sdf[, colnames(sdf) %in% iogenes] )
mut_io(decon, sdf$mut, sdf$Cohort, sdf$totalmut)
}

polycnv.orig = function(cnv, df = rnacnv) {
sdf = estimable_cnv(cnv, df = df)
sdf = droplevels(sdf[complete.cases(sdf), ])
decon = c(sdf[, grepl(".sig$", colnames(sdf))], iogenes)
cnv_io(decon, sdf$cnv, sdf$Cohort)
}

polycnv = function(cnv, df = rnacnv, cna ) {
if ( cna == 'l' ) { ## l = loss
    sdf = df[which(df[ , cnv] < cna.cutoff), ]
}
if ( cna == 'g' ) { ## g = gain
    sdf = df[which(df[ , cnv] > -cna.cutoff), ]
}
sdf = estimable_cnv(cnv, df = sdf)
if(is.null(sdf)){print("WE GOT A NULL");return(NULL)}
sdf = droplevels(sdf[complete.cases(sdf), ])
decon = cbind ( sdf[, grepl(".sig$", colnames(sdf))] , sdf[, colnames(sdf) %in% iogenes] )
##print(dim(decon))
out = cnv_io(decon, sdf$cnv, sdf$Cohort, sdf$Estimate)
out
}

################################################################
## calculate all ANOVA results, using the above functions
##

if ( compute.anova ) {
    ## calculate everything!
    if ( paper.cohorts.only ) {
        iogenes = c()
    }

p.1 = ldply(names(topmutated), .parallel=TRUE, failwith(NULL, function(x) {
    out = polymut(x, df = rnamut)
    if(is.null(out)){return(NULL)}
    cbind(data.frame(Variant = x), out)
}))

p.1$type = 'mut'

p.1f = ldply(allfmut, .parallel=TRUE, failwith(NULL, function(x) {
    out = polymut(x, df = rnafmut)
    if(is.null(out)){return(NULL)}
    cbind(data.frame(Variant = x), out)
}))
p.1f$type = 'fmut'

p.2 = ldply(names(allcnv),  .parallel=TRUE, failwith(NULL, function(x) {
    cbind(data.frame(Variant = x), polycnv(x, df = rnacnv, cna='l'))
}))
p.2$type = 'cnv'
p.2$Variant = gsub("\\.cnv$", ".l", p.2$Variant)

p.3 = ldply(names(allcnv),  .parallel=TRUE,  failwith(NULL, function(x) {
    cbind(data.frame(Variant = x), polycnv(x, df = rnacnv, cna='g'))
}))
p.3$type = 'cnv'
p.3$Variant = gsub("\\.cnv$", ".g", p.3$Variant)

df = merge(subcohorts, rnacnv, by.x = 'sampleID', by.y = 'Sample')
df$Cohort = df$subcohort.x
print(dim(df))

p.4 = ldply(names(allcnv),  .parallel=TRUE,  failwith(NULL, function(x) {
    cbind(data.frame(Variant = x), polycnv(x, df = df, cna='g'))
}))
dim(p.4)

p.4$type = 'cnv'
p.4$Variant = gsub("\\.cnv$", ".g", p.4$Variant)

p.5 = ldply(names(allcnv),  .parallel=TRUE,  failwith(NULL, function(x) {
    cbind(data.frame(Variant = x), polycnv(x, df = df, cna='l'))
}))
p.5$type = 'cnv'
p.5$Variant = gsub("\\.cnv$", ".l", p.5$Variant)

df = merge(subcohorts, rnamut, by.x = 'sampleID', by.y = 'Sample')
df$Cohort = df$subcohort.x


p.6 = ldply(names(topmutated), .parallel=TRUE, failwith(NULL, function(x) {
    out = polymut(x, df = df)
    if(is.null(out)){return(NULL)}
    cbind(data.frame(Variant = x), out)
}))

p.6$type = 'mut'

df = merge(subcohorts, rnafmut, by.x = 'sampleID', by.y = 'Sample')
df$Cohort = df$subcohort.x
print(dim(df))

p.6f = ldply(allfmut, .parallel=TRUE, failwith(NULL, function(x) {
    out = polymut(x, df = df)
    if(is.null(out)){return(NULL)}
    cbind(data.frame(Variant = x), out)
}))

p.6f$type = 'fmut'

df = merge(subcohorts, rnacnv, by.x = 'sampleID', by.y = 'Sample')
df$Cohort = df$subcohort.x

p.7 = ldply(names(allcnv),  .parallel=TRUE,  failwith(NULL, function(x) {
    out = polycnv(x, df = df, cna='l')
    if(is.null(out)){return(NULL)}
    cbind(data.frame(Variant = x), out)
}))

p.7$type = 'cnv'
p.7$Variant = gsub("\\.cnv$", ".l", p.7$Variant)


p.8 = ldply(names(allcnv),  .parallel=TRUE,  failwith(NULL, function(x) {
    out = polycnv(x, df = df, cna='g')
    if(is.null(out)){return(NULL)}
    cbind(data.frame(Variant = x), out)
}))
p.8$type = 'cnv'
p.8$Variant = gsub("\\.cnv$", ".l", p.8$Variant)





###############################################################
## that's all the core compute
## everything else just massages things...
################################################################
o = rbind.fill(p.1, p.1f, p.2, p.3, p.4, p.5, p.6, p.6f, p.7, p.8)
## write analysis output
write.csv(o, file="full.final.unifiedmodel.csv", row.names = FALSE)
## also make a date-stamped file
write.csv(o, file=paste("full.final.unifiedmodel", Sys.Date(), "csv", sep="."), row.names = FALSE)
dim(o)
o.P = o[ , grepl("\\.P\\.Value$", colnames(o))]
o.P = o.P[ , order(colnames(o.P))]
o.P = -log10(o.P)
o.P = round(o.P, digits=1)
dim(o.P)
o = o[ , ! grepl("\\.P\\.Value$", colnames(o))]
o = o[ , ! grepl("\\.adj\\.P\\.Val$", colnames(o))]
o = o[ , ! grepl("\\.t$", colnames(o))]
dim(o)
o = o[ , ! grepl("\\.P\\.Value$", colnames(o))]
o.V = o[ , grepl("^V\\.", colnames(o))]
o.V = o.V[ , order(colnames(o.V))]
o.V = round(o.V, digits=2)
colnames(o.V)
o = o[ , ! grepl("^V\\.", colnames(o))]
colnames(o)
minP = apply( o.P, 1, max, na.rm=TRUE)
##o$Variant = gsub("\\.p\\.", "\\|", o$Variant)
variant = gsub("\\.d\\.", "\\-", o$Variant)
variant = o$Variant
variant = gsub("\\.p\\.", "\\|", o$Variant)
variant = gsub("\\.mut", "", variant)
variant = gsub("\\.cnv", "", variant)
variant = gsub("\\.l$", "", variant)
variant = gsub("\\.g$", "", variant)
chr = sym2chr(variant)
pos = abs(sym2loc(variant))
## tweak the output
o = cbind(
    data.frame(Variant = o$Variant, Type = o$type,Chr = chr, Pos = pos, Sig = o$Sig, F = round(o$F, digits=1), modP = round(-log10(o$P.Value), digits=1), minP = minP),
    o.P,
    o.V)
## how did a pipe get in fine mutant names?????!!!!
o$Variant = gsub('\\|', '.', o$Variant)
## write analysis output
save(o, file = paste0(data.dir,paste(filename, 'o', "Robject", sep=".")))
write.csv(o, file=paste0(data.dir,"final.unifiedmodel.csv"), row.names = FALSE)
## also make a date-stamped file
##write.csv(o, file=paste("final.unifiedmodel", Sys.Date(), "csv", sep="."))
## make a smaller.file
##write.csv(o[o$minP > 3, ], file = paste0(data.dir,"final.filtered.unifiedmodel.csv"))
}
if ( ! compute.anova )
    {
load(file = paste0(data.dir,paste(filename, "o", "Robject", sep=".")))
assign('o', o, envir = .GlobalEnv)
}






################################################################
## deprecated data structures - don't use
################################################################
## make melted p value table
## this is now only used in scatter and box plots to get the P value

if(FALSE) { ## no

if ( ! use.cached | ! exists('my.o') | force.reload ) {
    my.o = cbind(
        o[ , 1:7],
        o[ , grepl("\\.P\\.Value$", colnames(o) ) ]
        ##data.frame(O = o$minP)
        )
    melted = melt(my.o, id.vars=1:7)
    melted$variable = gsub("^V\\.", "", melted$variable)
    melted$variable = gsub("\\.P\\.Value$", "", melted$variable)
    melt.o.p = melted
    ## make melted effect size table
    my.o = o[ , ! grepl("\\.P\\.Value$", colnames(o) ) ]
    my.o = cbind(
        my.o[ , 1:7],
        my.o[ , grepl("^V\\.", colnames(my.o) ) ]
        )
    melted = melt(my.o, id.vars=1:7)
    melted$variable = gsub("^V\\.", "", melted$variable)
    melted$variable = gsub("\\.P\\.Value$", "", melted$variable)
    melt.o.e = melted
    best.p = function(variants, df = melt.o.p) {
        my.o = df[df$Variant %in% variants , ]
        my.o = my.o[order(my.o$minP, decreasing=TRUE), ]
        my.o = my.o[1, ]
        my.o
    }
}

}
################################################################
## make a good data structure for analytics (melted)
################################################################
gc()


if (! use.cached ) {

colnames(o) = gsub('\\-', '.', colnames(o))

o.g  = ldply(  colnames(o)[grepl("\\.P\\.Value$", colnames(o))] ,  function(h) {
    i = gsub("\\.P\\.Value$", "", h)
    print(i)
    cohort = gsub("^V\\." , "", i)
    print(cohort)
    cbind(
        data.frame(
            Type = o$Type,
            Variant = o$Variant,
            Chr = as.factor(o$Chr),
            Pos = as.integer(o$Pos),
            Sig = o$Sig,
            minP = as.numeric(o$minP)  ),
        data.frame( Cohort = as.factor(cohort), P = as.numeric(o[ , h]) ),
        data.frame( Effect = as.numeric(o[ , i]) )
        )
})

if(remove.intermediate){
    ## damn, I have an old table that uses 'o'
    ##rm(o);
    gc()
}


## make Effect size  of loss inversed
o.g[grepl( '\\.l$', o.g$Variant), "Effect"] = -1 * o.g[grepl( '\\.l$', o.g$Variant), "Effect"]
####
##o.g$Variant = factor(o.g$Variant)
o.g$Dir = factor(ifelse( o.g$Effect > 0, "Pos", "Neg"), levels = c("Pos", "Neg") )
o.g = o.g[ which(! is.na(o.g$P) ), ]

o.g.o = o.g  ## keep copy of unfiltered

## now that we've made copies of the full data set, filter.
## sorry about the next line
##o.g = o.g[gsub("\\.l$", ".cnv", gsub("\\.g$", ".cnv", o.g$Variant)) %in%
##    c(names(topcnv), names(topmutated)), ]
## remove effect size for insignificant results - very helpful for graphing
## don't need to do this - manhat2 already does it
## o.g$Effect[o.g$P < effect.cutoff.P ]  = 0

## nathan recent comment
##   o.g = o.g[ which(o.g$Sig %in% good.sigs) , ]
##  o.g$Sig = factor(o.g$Sig, levels = good.sigs)

##o.g = o.g[ which(o.g$Variant %in% names(filtered.variants) ) , ]
## make a protein coding data set.
o.g = o.g[! grepl("^SNOR", o.g$Variant), ]
o.g = o.g[! grepl("^MIR", o.g$Variant), ]
o.g = droplevels(o.g)



if ( adjust.P ) {
    out = o.g
    out$P.orig = out$P
    print(summary(out$P))
    out$P = 10 ** ( -1 * out$P )
    blah = out[out$P >= 0.01, ]
    muts = out[out$P < 0.01 & out$Type != 'cnv', ]
    cnvs = out[out$P < 0.01 & out$Type == 'cnv', ]
    blah$P = 1
    muts$P = sapply(muts$P, p.adjust, n=topmut.cutoff, method = 'holm')
    cnvs$P = sapply(cnvs$P, p.adjust, n=topcnv.cutoff, method = 'holm')
    out = rbind(muts, cnvs, blah)
}
out$P = round(-1 * log10(out$P), digits = 1)
o.g = out



assign('o.g', o.g, envir = .GlobalEnv)
assign('o.g.o', o.g.o, envir = .GlobalEnv)

write.csv(o.g, file = paste0(data.dir,"o.g.csv"))
write.csv(o.g.o, file = paste0(data.dir,"o.g.o.csv"))
save(o.g, file = paste0(data.dir,"o.g.Robject"))
save(o.g.o, file = paste0(data.dir,"o.g.o.Robject"))

dim(o.g)
o.g.p13 = o.g[o.g$P > 1.3, ]
dim(o.g.p13)
save(o.g.p13, file = paste0(data.dir,"o.g.p13.Robject"))
assign('o.g.p13', o.g.p13, envir = .GlobalEnv)


} else {
load(file = paste0(data.dir,"o.g.Robject"))
load(file = paste0(data.dir,"o.g.o.Robject"))
load(file = paste0(data.dir,"o.g.p13.Robject"))
assign('o.g', o.g, envir = .GlobalEnv)
assign('o.g.o', o.g.o, envir = .GlobalEnv)
assign('o.g.p13', o.g.p13, envir = .GlobalEnv)

}

## filter results with poor signature correlations
if ( ! use.cached |  force.reload ) {

    o.g.sig = merge(sig.cut, o.g)

    if ( paper.cohorts.only ) {

        o.g = o.g[ which (! o.g$Cohort %in% paper.cohorts.not ) , ]
        o.g.sig = o.g.sig[ which (! o.g.sig$Cohort %in% paper.cohorts.not ) , ]
        o.g.sig = o.g.sig [ which (o.g.sig$Sig %in% paper.sigs) , ]
    }

    o.g.sig = o.g.sig[ ! duplicated( o.g.sig ) , ]
    save(o.g.sig, file = paste0(data.dir,"o.g.sig.Robject") )
    write.csv(o.g.sig, file = "Sdata1.csv", row.names = FALSE)


    o.g.cnv = o.g.sig [ which(o.g.sig$Type == 'cnv'),  ]
    o.g.cnv = droplevels(o.g.cnv)
    o.g.cnv$cnv = gsub("\\..$", ".cnv", o.g.cnv$Variant)
    o.g.cnv$gene = gsub("\\..$", "", o.g.cnv$Variant)

    o.g.cor.tmp = ddply(.parallel = TRUE, o.g.cnv, .(Cohort), failwith( NULL, function(y) {
        my.colnames = colnames(rnacnv)
        my.cohort = y[1, "Cohort"]
        print(my.cohort)
        sub.cohort = rnacnv[rnacnv$Cohort == my.cohort , ]
        gc()
        ddply(y, .(cnv), function(x) {
            my.cohort = x[1, "Cohort"]
            my.sig = x[1, "Sig"]
            my.variant = x[1, "Variant"]
            my.gene = x[1, "gene"]
            my.cnv = x[1, 'cnv']
            if( length( which(  c(my.gene, my.cnv) %in% my.colnames ) )  < 2  )  {
                return(NULL)
            } else {
                rnacnvcor = round(digits = 2,
                    cor( use = 'pairwise.complete.obs',
                        sub.cohort[, my.gene], sub.cohort[, my.cnv] ) )
                return(
                    data.frame(
                        Cohort = my.cohort,
                        cnv =   my.cnv,
                        rnacnvcor = rnacnvcor ) )
            }
        })
    }))


    o.g.cor.next = merge(o.g.cnv, o.g.cor.tmp, by = c('Cohort', 'cnv')  )

    o.g.cor.rnasigcor = ddply(.parallel = TRUE, o.g.cor.next, .(Cohort), failwith( NULL, function(y) {
        my.colnames = colnames(rnacnv)
        my.cohort = y[1, "Cohort"]
        print(my.cohort)
        sub.cohort = rnacnv[rnacnv$Cohort == my.cohort , ]
        gc()
        ddply(y, .(Sig, cnv), function(x) {
            my.cohort = x[1, "Cohort"]
            my.sig = x[1, "Sig"]
            my.variant = x[1, "Variant"]
            my.gene = x[1, "gene"]
            my.cnv =  x[1, "cnv"]
            if( length( which(  c(my.gene, my.cnv) %in% my.colnames ) )  < 2  )  {
                return(NULL)
            } else {
                rnasigcor = round(digits = 2,
                    cor( use = 'pairwise.complete.obs',
                        sub.cohort[, my.gene], sub.cohort[, my.sig] ) )
                return(
                    data.frame(
                        Cohort = my.cohort,
                        cnv = my.cnv,
                        rnasigcor = rnasigcor ) )
            }
        })
    }))



    o.g.rich = merge(o.g.cor.next, o.g.cor.rnasigcor, by = c('Sig', 'Cohort', 'cnv')  )
    o.g.rich = o.g.rich [ ! duplicated(o.g.rich) , ]
    save(o.g.rich, file = paste0(data.dir, 'o.g.rich.Robject') )
    ##rm(o.g.cnv, o.g.cor.tmp, o.g.cor.next, o.g.cor.rnasigcor)
} else {
    load(file = paste0(data.dir, 'o.g.rich.Robject') )
    load(file = paste0(data.dir, 'o.g.sig.Robject') )
}





#v.p.view.threshold = 0; p.effect.threshold = 0; effect.cutoff.P = 0 ; good.sigs = c('TCD8.sig')


## bottom.label = function(df, labels = list(), l.pos = FALSE) {
##     ## need to place labels at bottom of graph
##     max.P = max(abs(df$P), na.rm=TRUE)
##     df$P = -1 * max.P
##     if(! l.pos == FALSE) { df$P = l.pos }
##     labels = gsub("\\..*", "", labels)
##     df$Variant = gsub("\\..*", "", df$Variant)
##     df = df[df$Variant %in% labels, ]
##     ## now we want unique??? a little voodoo...
##     df[match(unique(df$Variant), df$Variant), ]
## }


manhat3 = function(dat, my.y = "Effect", my.color = 'rnasigcor', my.shape = 'Type', my.size = 'P', variant = FALSE, hjust=0, vjust=NA, lSize = 3, rotate = 90 ) {
    dat = dat[complete.cases(dat), ]
    if(typeof( dat[ , my.color]) == 'double' | typeof( dat[ , my.color]) == 'integer' ) {
        my.scale = scale_color_gradient()
    } else {
        my.scale = scale_color_manual(values=c("#3366CC", "#DC3912", "#FF9900", "#109618",
                                               "#990099", "#0099C6",  "#000000",  "#AAAA11",
                                               "#6633CC", "#E67300", "#8B0707", "#651067",
                                               "#329262", "#5574A6", "#3B3EAC", "#DD4477",
                                               "#66AA00"))
    }
    dat$Type = gsub("^.*\\.", "", dat$Variant)
    dat[dat$Type == 'mut', "Type"]  = 'Mutant'
    dat[dat$Type == 'l', "Type"]  = 'Loss'
    dat[dat$Type == 'g', "Type"]  = 'Gain'
    dat[dat$Type == 'cnv', "Type"]  = 'CNA'
    dat[dat$Type == 'sig', "Type"]  = 'Sig'
    dat$Type = factor(dat$Type, levels = c("Mutant", "Loss", "Gain", "CNA", "Sig"))
    chromo = dat[, 'Chr'][[1]]
    cohort  = paste(unique(dat[, 'Cohort']), collapse = ',')
    dat$Pos = round(dat$Pos / 1000000, digits = 2)
                                        #dat$P = ifelse(dat$Effect > 0, dat$P, dat$P * -1)
    print(head(dat), 1)
    p = ggplot(aes_string(x="Pos", y=my.y, text = "Variant", size=my.size,  color = my.color, shape=my.shape), data=dat) + my.scale + geom_point(fill = NA) +  scale_shape_manual(values=c(8,1,2)) + theme_bw() +
        theme( legend.position = 'left',
              plot.background = element_rect(fill = 'white'),
              legend.text = element_text(size = 10),
              legend.background = element_rect(fill = 'white'),
              legend.key=element_rect(fill='NA'),
              axis.text = element_text(size = 8, vjust=1)) +
                  guides(colour = guide_legend(override.aes = list(size=5))) +
                      guides(shape=guide_legend(override.aes=list(size=5))) +
                          scale_shape_discrete(solid=F)
    if ( ! variant == FALSE ) {
        my.min = min(dat[ , my.y], na.rm = TRUE)
        variant.dat = dat[dat$Variant %in% variant, ]
        if(! is.na(vjust) ) {
            my.min = vjust
        }
        variant.dat[ , my.y ] = my.min
        variant.dat$Variant = gsub('\\.cnv$', '', (gsub('\\.[lg]$', '', gsub('\\.mut$', '', variant.dat$Variant))))
        variant.dat = variant.dat [ ! duplicated(variant.dat$Variant) , ]
        print ( head ( variant.dat ) )
        p = p + geom_text(
            aes_string(x='Pos', y=my.y, label="Variant"),
            angle=rotate,
            size= lSize, color='black',
            data = variant.dat )
    }
    ##               ifelse(variant == FALSE, NULL,
    ##                      geom_point(data = dat[dat$Variant == variant & dat$Cohort == cohort, ],
    ##                                 shape = 0, size=7, color='red'   ))  +
    ##                                     theme(legend.position = 'bottom') +
    ##                                         ggtitle(paste('Interactions', 'on Chr', chromo, 'in', cohort))
    p
}

## filter out non paper cohorts and non paper sigs too!





## make a gistic data set

if ( make.gistic ) {
mydf = cbind(data.frame(Cohort = rnacnv$Cohort), rnacnv[ , grepl("\\.cnv$", colnames(rnacnv)) ] )
mygistic = ddply(mydf, .(Cohort), .parallel = TRUE, function(x) {
    ##print(x[1, 1:3])
    melted = melt(apply(x[, 2:length(mydf)], 2, quantile, na.rm=TRUE, c(0.005, 0.01, 0.25, 0.5,  0.75, 0.99, 0.995) ))
    colnames(melted) = c("Sig", "Variant", "P")
    melted
} )
mygistic$Chr = sym2chr(nodot(nopipe(mygistic$Variant)))
mygistic$Pos = abs(sym2loc(nodot(nopipe(mygistic$Variant))))
mygistic$minP = mygistic$P
mygistic$Effect = mygistic$P
gc(); gc()
}

## make sig quantile data frame
## no, this is flawed and does not work

## mydf = cbind(data.frame(Cohort = rna$Cohort), createdecon(rna) )
## myqsig = ddply(mydf, .(Cohort), function(x) {
##     print(x[1, 1:3])
##     melted = melt(apply(x[, 2:length(mydf)], 2, quantile, na.rm=TRUE, c(0.02, 0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95, 0.98) ))
##     colnames(melted) = c("Sig", "Variant", "P")
##     melted
##     } )
## myqsig$Chr = sym2chr(nodot(nopipe(myqsig$Variant)))
## myqsig$Pos = abs(sym2loc(nodot(nopipe(myqsig$Variant))))
## myqsig$minP = myqsig$P
## myqsig$Effect = myqsig$P

##  subset data frames to give to ggplot to define
##  which data points are to be labeled

best.label = function(variants, df = o) {
df = df[ df$Variant %in% variants , ]
df = ddply(df, .(Variant), function (x) {
    out = x[order(x$Variant, decreasing = TRUE), ]
    out[1, ]
})
df
}

best.label.e = function(variants, df = melt.o.e) {
df = df[ df$Variant %in% variants , ]
df = ddply(df, .(Variant), function (x) {
    out = x[order(x$Variant, decreasing = TRUE), ]
    out[1, ]
})
df
}

## similar to above but takes df as input, not just list of variants
## subset data frame to plot labels in manhat2
get.b = function(df, labels = list() ) {
odf = df
if(length(labels) > 0 )  {
    df = df[df$Variant %in% labels, ]
}
df = ddply(df, .(Variant), function(x) {
    x = x[order(abs(x$P), decreasing = TRUE, na.last = NA) , ]
    x[ 1, , drop=FALSE]
})
if(  length(labels)  == 0  )  {df = odf[1, , drop=FALSE]; df = df[-1, , drop=FALSE] }
df$Variant = gsub("\\.cnv", ".c", df$Variant)
df$Variant = gsub("\\.mut", ".m", df$Variant)
df
}

## try to add gene labels at bottom of graph
bottom.label = function(df, labels = list(), l.pos = FALSE) {
## need to place labels at bottom of graph
max.P = max(abs(df$P), na.rm=TRUE)
df$P = -1 * max.P
if(! l.pos == FALSE) { df$P = l.pos }
labels = gsub("\\..*", "", labels)
df$Variant = gsub("\\..*", "", df$Variant)
df = df[df$Variant %in% labels, ]
## now we want unique??? a little voodoo...
df[match(unique(df$Variant), df$Variant), ]
}



## we have serious memory problems with ggplot now
## that we have included all CNVs...
## this function takes a standard output df like o.g
## and reduces it to only the data point for each gene and signature that
## has the best P value.
## early tests suggest this is just as bad as ggplot...

## best.gene.sig = return only the top n hits of any Variant/Signature combination
## good for multi-cohort plots
best.gene.sig = function(df, n=1) {
ddply(df, .(Variant, Sig), function(x) {
    x = droplevels(x)
    x = x[order(abs(x$P), decreasing = TRUE, na.last=NA), ]
    x[1:n , , drop=FALSE]
})
}




################################################################
## manhat2 - genomic manhattan plotting of results

## any o.g-form subset
################################################################

manhat2 = function(ndf, best.only = FALSE, p.view.threshold = v.p.view.threshold, p.effect.threshold = effect.cutoff.P, sigs = good.sigs, cohort=NA, gistic=FALSE, labels = list(), rotate=90, title = '' , hjust=0, vjust=0.5, l.pos = -30, size=NA, lSize = 3, posneg = TRUE, min = 1, max = FALSE, minSize = 0.5, maxSize = 8) {
##ndf = droplevels(ndf[! is.na(ndf$Chr), ])
ndf = ndf[which(! is.na(ndf$P) ), ]
ndf = ndf[which(! is.na(ndf$Chr) ), ]
ndf = ndf[which(ndf$P > p.view.threshold), ]
if (best.only) {   ndf = best.gene.sig(ndf)   }
ndf$Type = gsub("^.*\\.", "", ndf$Variant)
ndf[ndf$Type == 'mut', "Type"]  = 'Mutant'
ndf[ndf$Type == 'l', "Type"]  = 'Loss'
ndf[ndf$Type == 'g', "Type"]  = 'Gain'
ndf[ndf$Type == 'cnv', "Type"]  = 'CNA'
ndf[ndf$Type == 'sig', "Type"]  = 'Sig'
ndf$Type = factor(ndf$Type, levels = c("Mutant", "Loss", "Gain", "CNA", "Sig"))
gc()
if (! sigs[[1]] == 'all') {
    ndf = ndf[which(ndf$Sig %in% sigs), ]
    ndf$Sig = factor(gsub(".sig", "", ndf$Sig), levels = gsub(".sig", "", good.sigs) )
}
if ( sigs[[1]] == 'all') {
    ndf$Sig = factor(ndf$Sig)
}
ndf$Sig = droplevels(ndf$Sig)
ndf$Chr = factor(ndf$Chr, levels = c(
                              "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                              "11", "12", "13", "14", "15", "16",
                              "17", "18", "19", "20", "21", "22", "X", "Y") )
ndf$Pos = ndf$Pos / 1000000
if(posneg) {
    ndf[grepl( '^Neg$', ndf$Dir), "P"] = -1 * ndf[grepl( '^Neg$', ndf$Dir), "P"]
}
ndf[which(abs(ndf$P) < p.effect.threshold), "Effect"] = 0
ndf = droplevels(ndf)
plot.y.min = min(ndf$P, na.rm=TRUE)
plot.y.max = max(ndf$P, na.rm=TRUE)
if(plot.y.min > l.pos) {plot.y.min = l.pos}
gc()
p =
    ggplot(aes(x=Pos, y=P, color=Sig, size=abs(Effect), shape=Type), data=ndf) +
         scale_size_continuous(range = c(minSize, maxSize)) +
        geom_point() +
            ggtitle(title) +
                scale_shape_manual(values=c(8,1,2)) +
                    facet_wrap(~ Chr, scales='free_x') +
                        scale_color_manual(values=c("#3366CC", "#DC3912", "#FF9900", "#109618",
                                               "#990099", "#0099C6",  "#000000",  "#AAAA11",
                                               "#6633CC", "#E67300", "#8B0707", "#651067",
                                               "#329262", "#5574A6", "#3B3EAC", "#DD4477",
                                               "#66AA00")) +
theme( legend.position = 'left',
      plot.background = element_rect(fill = 'white'),
      legend.text = element_text(size = 8),
      legend.background = element_rect(fill = 'white'),
      legend.key=element_rect(fill='NA'),
      axis.text = element_text(size = 6, vjust=1)) +
          guides(colour = guide_legend(override.aes = list(size=5))) +
              guides(shape=guide_legend(override.aes=list(size=5))) +
                  geom_text(aes(x=Pos, y=P, label=Variant),
                            hjust=hjust, vjust=vjust, angle=rotate,
                            size= lSize, color='black',
                            data = bottom.label(ndf, labels=labels, l.pos=l.pos) ) +
                                ylab("-logP (sign: direction of effect)") +
                                    xlab("Chromosomal Position (megabases)")
## minor breaks not rendering right now...
if(! max == FALSE) {
    p = p  + scale_x_continuous(limits=c(min/ 1000000 , max / 1000000 ), minor_breaks = seq(0, 1000, 1) )
}
p = p + scale_y_continuous(limits=c(plot.y.min, plot.y.max))
## if(max == FALSE) {
##     p = p  + scale_x_continuous(minor_breaks = seq(0, 1000, 1) )
## }
rm(ndf)
gc()
p
}


## simplified manhattan plot

manhat4 = function ( data, p.thresh = 1, effect.thresh = 0.3, max.P = 30, mut.p.thresh = 1,
    noxy = TRUE, facet.legend = TRUE, make.sd = TRUE, show.mut = TRUE,
    scale.legend = FALSE, yscale = FALSE, small.theme = TRUE) {
    data$Chr = as.character(data$Chr)
    data$Chr = paste("Chromosome", data$Chr)
    print(str(data$Chr))
    data$Chr = factor(data$Chr, levels = paste("Chromosome", c(
                              "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                              "11", "12", "13", "14", "15", "16",
                                    "17", "18", "19", "20", "21", "22", "X", "Y") ) )
    droplevels(data)
    print(unique(levels(data$Chr)))
    print(str(data$Chr))
    if( noxy ) {
        data = data[data$Chr %in% paste("Chromosome", 1:22) , ]
    } else {
        data = data[data$Chr %in% paste("Chromosome", c(1:22, X, Y)) , ]
    }
    maxmin.data = ldply( unique(data$Chr) , function ( chrom ) {
        mymax = max( data[which(data$Chr == chrom), "Pos"], na.rm = TRUE ) / 1000000
        data.frame(Type = 'maxmin', Variant = 'none', Chr = chrom, Pos = c(0,mymax), Sig = 'none', P = 0, Effect = 0)
    })
    print(maxmin.data)
    data = data[data$P > p.thresh, ]
    data = data[abs(data$Effect) > effect.thresh, ]
    print(summary(data$P))
    data$P = ifelse( data$P > max.P, max.P, data$P )
    print(summary(data$P))
    data$Pos = data$Pos / 1000000
    data$P = ifelse( grepl( '\\.l$', data$Variant) , -1 * data$P, data$P )
    cnvdat = data[data$Type == 'cnv', ]
    if ( make.sd ) cnvdat$P = cnvdat$P / sd(cnvdat$P)
    mutdat = data[data$Type == 'mut' & data$P > mut.p.thresh, ]
    if ( make.sd ) mutdat$P = mutdat$P / sd(mutdat$P)
    print(summary(data$P))
    p = ggplot(aes(x = Pos, y = P, color = Effect), dat = data)
    p = p + geom_point(dat = maxmin.data, size = 0)
    p = p + geom_bar( aes(fill = Effect), stat = 'identity',  position = 'dodge', width = 1, dat = cnvdat)
    if (show.mut) {
        p = p + geom_point( shape = 8, dat = mutdat )
        labeldat = mutdat
        labeldat$P = -6
        labeldat$Variant = gsub('\\.mut', '', labeldat$Variant)
        p = p +  geom_text(aes( label = Variant ), hjust = -0.1, size = 1.5, color = 'black', angle = 90, dat = labeldat )
    }
    p = p + theme_gdocs()
    if (! yscale ) p = p + scale_y_continuous(breaks = c(0))
    p = p +
                xlab('Chromosome position (MBases)') +
                    ylab('Association ( sign: CNA increase vs decrease )' )
    if(facet.legend) {
        l.cnv = data.frame(
            Chr = "Legend",
            Variant = c(
                'CNA increase',
                'CNA decrease',
                'Increased infiltrate',
                "",
                'Decreased infiltrate',
                "",
                "" ),
            Pos = c(10, 30, 50, 50, 75, 75, 175),
            P = c(5, -5,  3, -3, 2, -2, 0),
            Effect = c( 0, 0, 1, 1, -1, -1, 0)
            )
        p = p + geom_bar( aes(fill = Effect, color = Effect), stat = 'identity', width = 5, dat = l.cnv)
        p = p + geom_bar( color = 'black' , fill = 'white', stat = 'identity', width = 5, dat = l.cnv[1:2 , ] )
            #scale_fill_manual(values=c("#000000", "#bb0000", "#0000bb")) +
                                        #scale_fill_discreet(low = '#0000bb', high = '#bb0000', mid = '#000000', guide = FALSE) +
            ##scale_fill_discrete(low = '#0000bb', high = '#bb0000', mid = '#000000', guide = FALSE) +
        p = p + geom_text ( aes( label = Variant ), size = 1.5, color = 'black', angle = 0, dat = l.cnv, hjust = -0.1 )
        l.mut = data.frame(
            Chr = "Legend",
            Variant = 'Mutation',
            Pos = 100,
            P = 4.5,
            Effect = 1
            )
        p = p + geom_point( shape = 8, color = 'black', size = 2, dat = l.mut) +
            geom_text ( aes( label = Variant ), size = 1.5, color = 'black', angle = 0, dat = l.mut, hjust = -0.35 )
    }
    p = p + facet_wrap(~ Chr, scales = 'free_x') +
        scale_colour_gradient2(name = "Immune\nEffect", low = '#0000bb', high = '#bb0000') +
    scale_fill_gradient2(name = "Immune\nEffect", low = '#0000bb', high = '#bb0000') 
    if (small.theme ) {
        p = p + theme(
            strip.text.x = element_text(size = 7,  margin = margin(.1, 0, .1, 0, "cm")),
            axis.text = element_text(size = 5),
            legend.title = element_text(size = 10), legend.key.width=unit(0.2,"cm"),
            legend.key.height=unit(0.4,"cm") 
            )
    }
    if (! scale.legend ) p = p +theme(legend.position = 'none')
    p
}



## plotting cnvs by different parameters
manhat3 = function(dat, my.y = "Effect", my.color = 'rnasigcor', my.shape = 'Type', my.size = 'P', variant = FALSE, hjust=0, vjust=NA, lSize = 3, rotate = 90 ) {
    dat = dat[complete.cases(dat), ]
    if(typeof( dat[ , my.color]) == 'double' | typeof( dat[ , my.color]) == 'integer' ) {
        my.scale = scale_color_gradient()
    } else {
        my.scale = scale_color_manual(values=c("#3366CC", "#DC3912", "#FF9900", "#109618",
                                                "#990099", "#0099C6",  "#000000",  "#AAAA11",
                                                "#6633CC", "#E67300", "#8B0707", "#651067",
                                                "#329262", "#5574A6", "#3B3EAC", "#DD4477",
                                                "#66AA00"))
    }
    dat$Type = gsub("^.*\\.", "", dat$Variant)
    dat[dat$Type == 'mut', "Type"]  = 'Mutant'
    dat[dat$Type == 'l', "Type"]  = 'Loss'
    dat[dat$Type == 'g', "Type"]  = 'Gain'
    dat[dat$Type == 'cnv', "Type"]  = 'CNA'
    dat[dat$Type == 'sig', "Type"]  = 'Sig'
    dat$Type = factor(dat$Type, levels = c("Mutant", "Loss", "Gain", "CNA", "Sig"))
    droplevels(dat)
    chromo = dat[, 'Chr'][[1]]
    cohort  = paste(unique(dat[, 'Cohort']), collapse = ',')
    dat$Pos = round(dat$Pos / 1000000, digits = 2)
                                        #dat$P = ifelse(dat$Effect > 0, dat$P, dat$P * -1)
    print(head(dat, 1))
    p = ggplot(aes_string(x="Pos", y=my.y, text = "Variant", size=my.size,  color = my.color, shape=my.shape), data=dat) +
        my.scale + geom_point(fill = NA) +
            scale_shape_manual(values=c(8,1,2)) + theme_bw() +
        theme( legend.position = 'left',
              plot.background = element_rect(fill = 'white'),
              legend.text = element_text(size = 10),
              legend.background = element_rect(fill = 'white'),
              legend.key=element_rect(fill='NA'),
              axis.text = element_text(size = 8, vjust=1)) +
                  guides(colour = guide_legend(override.aes = list(size=5))) +
                      guides(shape=guide_legend(override.aes=list(size=5))) +
                          scale_shape_discrete(solid=F) +
                              xlab("Chromosomal Position (megabases)") +
                                  ylab("Effect of CNA on Infiltrate (S.D. units)")
    if ( ! variant == FALSE ) {
        my.min = min(dat[ , my.y], na.rm = TRUE)
        variant.dat = dat[dat$Variant %in% variant, ]
        if(! is.na(vjust) ) {
            my.min = vjust
        }
        variant.dat[ , my.y ] = my.min
        variant.dat$Variant = gsub('\\.cnv$', '', (gsub('\\.[lg]$', '', gsub('\\.mut$', '', variant.dat$Variant))))
        variant.dat = variant.dat [ ! duplicated(variant.dat$Variant) , ]
        print ( head ( variant.dat ) )
        p = p + geom_text(
            aes_string(x='Pos', y=my.y, label="Variant"),
            angle=rotate,
            size= lSize, color='black',
            data = variant.dat )
    }
    ##               ifelse(variant == FALSE, NULL,
    ##                      geom_point(data = dat[dat$Variant == variant & dat$Cohort == cohort, ],
    ##                                 shape = 0, size=7, color='red'   ))  +
    ##                                     theme(legend.position = 'bottom') +
    ##                                         ggtitle(paste('Interactions', 'on Chr', chromo, 'in', cohort))
    p
}









################################################################
## manhattan plot around a gene

submanhat = function(variant, region = 10000000, cohort = NA,
effect.threshold = effect.cutoff.P, labels=c(variant),
p.view.threshold = v.p.view.threshold
                     ) {
my.pos = o.g.o[o.g.o$Variant == variant , "Pos"][[1]]
my.chr = o.g.o[o.g.o$Variant == variant , "Chr"][[1]]
my.upper = my.pos + region
my.lower = my.pos - region
df = o.g.o[o.g.o$Chr == my.chr &
               o.g.o$Pos < my.upper &
                   o.g.o$Pos > my.lower , ]
if(! is.na(cohort) ){
    df = df[df$Cohort == cohort, ]
}
p = manhat2(df, effect.threshold = effect.cutoff.P)
rm(df)
gc()
p
}

################################################################
## annotated manhattan plots (gene labels on lower edge of graph
################################################################


special.genes = c("HLA-A", "CDKN2A", "JAK2", "TP53", "MYC", "ERBB2", 'APC.mut',  'MIR587.cnv',
'ITIH5.mut', 'FLT1.mut', 'NOX5.mut',
'TUBD1.cnv', 'ZFN726.cnv',
'SNORA51', 'ERBB2.cnv', 'ATRX.mut')

annot.manhat = function(df, mut.no=80, l.pos = -50, lSize = 2,
p.view.threshold = v.p.view.threshold,
labels = special.genes, title = ""  ) {
tdf = df[order(abs(df$P), decreasing=TRUE, na.last=NA) , ]
muts = tdf[grepl("\\.mut$", tdf$Variant), ]
muts = muts[which( abs(muts$P) > 4 ), ]
good.muts = head(as.character(muts$Variant), mut.no)
genes = c(special.genes, good.muts)
p = manhat2(df, p.view.threshold = p.view.threshold, title = title, vjust = 0.5, hjust=-0.03, l.pos = l.pos, lSize= lSize, labels=genes)
p
}


################################################################
## Tables
## function to make tables of hits from list of variants
################################################################

hit.kable2 = function(l, P = 1.3, sort.by = 'Effect', decreasing = TRUE, sigs = good.sigs, coords = TRUE, data = o.g.sig, df.only = FALSE, fmut = TRUE) {
    if ( fmut ) {
        my.fmut.grep = gsub('\\.mut$', '\\\\.p\\\\.', l)
        ##print(l)
        ##print(my.fmut.grep)
        my.regexp =  paste(my.fmut.grep, collapse="|")
        ##print(my.regexp)
        matches = unique (   grep(
            my.regexp,
            data$Variant, value=TRUE ) )
        ##print(matches)
        l = unique( c(l, matches) )
    }
    x = data[which(data$Variant %in% l) ,  ]
    if( sigs != 'all' ) { x = x[which(x$Sig %in% sigs), , drop=FALSE] }
    x = x[ which(x$P > P), , drop = FALSE]
    if (coords) {
        x = data.frame(
            Variant = x$Variant,
            Chr. = x$Chr,
            Pos. = round(x$Pos / 1000000, digits=1),
            Signature = gsub("\\.sig$", "", x$Sig),
            ## Cohort = paste(x$Cohort, ' (Effect (&#8209;logP)', sep=""),
            Cohort = x$Cohort,
            value = paste(x$Effect, ' (', x$P, ')', sep="")
            )
        x$Cohort = gsub('BRCA.Basal', 'Breast (Basal)', x$Cohort)
        x$Cohort = gsub('BRCA.LumB', 'Breast (LumB)', x$Cohort)
        x$Cohort = gsub('COAD.MSS', 'Colon (MSS)', x$Cohort)
        x = dcast(x, Variant + Chr. + Pos. + Signature ~ Cohort)
    } else {
        x = data.frame(
            Variant = x$Variant,
            Signature = gsub("\\.sig$", "", x$Sig),
            ##Cohort = paste(x$Cohort, ' Effect (&#8209;logP)', sep=""),
            Cohort = x$Cohort,
            value = paste(x$Effect, ' (', x$P, ')', sep="")
            )
        x$Cohort = gsub('BRCA.Basal', 'Breast (Basal)', x$Cohort)
        x$Cohort = gsub('BRCA.LumB', 'Breast (LumB)', x$Cohort)
        x$Cohort = gsub('COAD.MSS', 'Colon (MSS)', x$Cohort)
        x = dcast(x, Variant + Signature ~ Cohort)
    }
    colnames(x)  = gsub("_Carcinoma", "", gsub("_Cancer", "", colnames(x) ) )
    x$Variant = gsub("\\.l", " loss", gsub("\\.g", " gain", gsub("\\.mut", " mutant", x$Variant)))
    out = x
    out[is.na(out)] = ""
    colnames(out) = gsub('_', " ", colnames(out))
    out = ddply(out, .(Variant), function(y) { y = y[order(y$Signature, decreasing = FALSE), ]; y })
    out$Variant = gsub('.p.', ' ', out$Variant, fixed = TRUE)
    out$Variant = gsub(".st.", ' stop', out$Variant, fixed = TRUE)
    out$Variant = gsub("_splice", ' splice', out$Variant, fixed = TRUE)
    if(df.only) {
        out
    } else {kable(out) }
}

## function to make table from filtered data frame
hitdf.kable2 = function(df = o.g.sig, P = 1.3, sort.by = 'Effect', decreasing = TRUE, df.only = FALSE) {
x = df
x = x[which(x$P > P) , ]
x = data.frame(
    Variant = x$Variant,
    Chr = x$Chr,
    Pos = round(x$Pos / 1000000, digits=1),
    Sig = gsub("\\.sig$", "", x$Sig),
    Cohort = x$Cohort,
    ##Cohort = paste(x$Cohort, ' <br/> effect(&#8209;logP)', sep=""),
    value = paste(x$Effect, ' (', x$P, ')', sep="")
    )
x$Cohort = gsub("_Carcinoma", "", gsub("_Cancer", "", x$Cohort) )
x$Variant = gsub("\\.l", " loss", gsub("\\.g", " gain", gsub("\\.mut", " mutant", x$Variant)))
out = dcast(x, Variant + Chr + Pos + Sig ~ Cohort)
out[is.na(out)] = ""
colnames(out) = gsub('_', " ", colnames(out))
out = ddply(out, .(Variant), function(y) { y = y[order(y$Sig, decreasing = TRUE), ]; y })
out$Variant = gsub('.p.', ' ', out$Variant, fixed = TRUE)
out$Variant = gsub(".st.", ' stop', out$Variant, fixed = TRUE)
out$Variant = gsub("_splice", ' splice', out$Variant, fixed = TRUE)
out$Cohort = gsub('BRCA_Basal', 'Breast (Basal)', out$Cohort)
out$Cohort = gsub('BRCA_LumB', 'Breast (LumB)', out$Cohort)
out$Cohort = gsub('COAD_MSS', 'Colon (MSS)', out$Cohort)
rm(x)
gc()
if(df.only) {
    out
} else {kable(out) }
}

gc()
