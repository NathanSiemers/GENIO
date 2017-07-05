source('one.set.vars.R')

pheno.url = 'https://icgc.xenahubs.net/download/donor.all_projects.xena'
pheno.meta = 'https://icgc.xenahubs.net/download/donor.all_projects.xena.json'
mut.url = 'https://icgc.xenahubs.net/download/SNV.donor.codingMutation-allProjects'
mut.meta = 'https://icgc.xenahubs.net/download/SNV.donor.codingMutation-allProjects.json'
mutfile = paste0(data.dir, 'icgc.mut')
mutmetafile = paste0(data.dir, 'icgc.mut.json')
phenofile = paste0(data.dir, 'icgc.pheno')
phenometafile = paste0(data.dir, 'icgc.pheno.json')

do.download = TRUE
system('mkdir data')

if (do.download) {
    download.file(pheno.url, phenofile, cacheOK= FALSE, method = 'curl')
    download.file(pheno.meta, phenometafile, cacheOK= FALSE, method = 'curl')
    download.file(mut.url, mutfile, cacheOK= FALSE, method = 'curl')
    download.file(mut.meta, mutmetafile, cacheOK= FALSE, method = 'curl')
}

## patient/sample data (pheno)

pheno = read.csv(sep = '\t', quote = '', header = TRUE, stringsAsFactors = FALSE, file = phenofile )
dim(pheno)
head(pheno)
## Non-US only
pheno = pheno[! grepl('US$', pheno$project_code) , ]
dim(pheno)
head(pheno)

## mutation data

mut = read.csv(sep = '\t', quote = '', header = TRUE, stringsAsFactors = FALSE, file = mutfile )
dim(mut)
mut[1:5, ]

## mutation load per sample

mut.load = ddply(mut, ~ sample, summarize, count = length(start) )
head(mut.load)

## gene level mutation

mut.genelevel = ddply(mut, ~ sample, function(s) {
    dcast(s, sample ~ gene, value.var = 'gene', fun.aggregate = length)
})
dim(mut.genelevel)

## NAs become 0 (that's really ok); muts are only 1
mut.genelevel[1:5, 1:5]
mut.genelevel[is.na(mut.genelevel)] = 0
mut.genelevel[1:5, 1:5]
sample.tmp = mut.genelevel$sample
mut.genelevel[mut.genelevel > 1] = 1
mut.genelevel$sample = sample.tmp
mut.genelevel[1:5, 1:10]

        

