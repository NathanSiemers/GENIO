#!/bin/csh

cp biblio.xml biblio.medline

pandoc-citeproc -y biblio.medline > biblio.yaml
pandoc-citeproc -y Beck-2016.bib > Beck-2016.yaml
pandoc-citeproc -y R.citations.bib > R.citations.yaml

pandoc  --csl=csl.csl --reference-docx reference.docx --bibliography R.citations.bib --bibliography biblio.medline --bibliography Beck-2016.bib manuscript.md -o manuscript.docx

mkdir Final

cp manuscript.docx Final/manuscript.docx
# sometimes figure names not written correctly
cp figure/unnamed-chunk-1.tiff Final/Fig-1.tiff
cp figure/unnamed-chunk-2.tiff Final/Fig-2.tiff
cp figure/unnamed-chunk-3.tiff Final/Fig-3.tiff
cp figure/unnamed-chunk-4.tiff Final/Fig-4.tiff
cp figure/unnamed-chunk-5.tiff Final/Fig-5.tiff
cp figure/unnamed-chunk-6.tiff Final/Fig-6.tiff
cp figure/Fig-1.tiff Final
cp figure/Fig-2.tiff Final
cp figure/Fig-3.tiff Final
cp figure/Fig-4.tiff Final
cp figure/Fig-5.tiff Final
cp figure/Fig-6.tiff Final
#cp figure/unnamed-chunk-7.tiff Final/unnamed-chunk-7.tiff
##cp figure/unnamed-chunk-7.tiff Final
##cp figure/unnamed-chunk-8.tiff Final
##cp figure/unnamed-chunk-9.tiff Final
cp figure/SFig-1.tiff Final
cp Sdata1.csv.gz Final
cp manuscript.Rmd Final
cp one.set.vars.R Final
cp two.main.compute.R Final








