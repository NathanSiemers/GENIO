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
