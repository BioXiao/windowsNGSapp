"""
exptf.py
A galaxy tool to find the top-n highest expressed TFs
"""
import subprocess
import os
import sys
import tempfile
#import htmlParser
import csv

galhtmlprefix = """<?xml version="1.0" encoding="utf-8" ?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<meta name="generator" content="Galaxy %s tool output - see http://galaxyproject.org/" />
<title>Conduct GO</title>
<link rel="stylesheet" href="/static/style/base.css" type="text/css" />
</head>
<body><br><br>

"""

galhtmlpostfix = """</body></html>"""

rprog="""rgo=function(tdir, title, txtfname, logfpath, anntIdentifier, rUtilName)
{
sink(file=file(logfpath, "at"), type="message")
sink(file=file(logfpath, "at"), type="output")
if (anntIdentifier == "fly.db0.db"){anntIdentifier = "fly.db0"}

newfiles = c()
newnames = c()
newtypes = c()

cdir=getwd()
setwd(tdir)

pkgTest <- function(x)
  {
    if (!require(x,character.only = TRUE))
    {
      #install.packages(x,dep=TRUE)
      source("http://bioconductor.org/biocLite.R")
      biocLite(x)
      if(!require(x,character.only = TRUE)) stop("Package not found")
    }
  }


 #source("http://bioconductor.org/biocLite.R")
 #biocLite("GO.db")
 
#source("http://bioconductor.org/biocLite.R")
#biocLite("GOstats")

#library("GOstats")
pkgTest("GOstats")

#library(GO.db)
pkgTest("GO.db")

pkgTest(anntIdentifier)

library(xtable)
#source(rUtilName)

print(txtfname)

dummy <- read.table(txtfname, sep = "\t", header = TRUE, as.is = TRUE, comment.char = "#", quote = "", check.names = FALSE)

selGenes <- unique(as.character(dummy[, "Gene"]))
selGenes <- as.character(dummy[, 1])
selGenes <- selGenes[!is.na(selGenes)]
selGenes <- selGenes[!selGenes == ""]
selGenes <- selGenes[!selGenes == "\n"]
selGenes <- selGenes[!selGenes == "\r"]
universeGenes <- NULL

# GO analysis, Cellular Component ontology
#annt <- "hgu133a"

#print(paste("THIS IS FROM OUTSIDE starting GO",annt))

#cellular component
outfname = "CC_RESULT"
outfname = paste(title, outfname, sep="_")
outfname = paste(outfname, "txt", sep=".")
outffname = file.path(tdir, outfname) #provides full file path

params <- new("GOHyperGParams", geneIds = selGenes, universeGeneIds = universeGenes,
annotation = anntIdentifier , ontology = "CC", pvalueCutoff = 0.05,conditional = FALSE, testDirection = "over")
#print(paste("annotation:", annotation)
GOResultCC <- hyperGTest(params)

#Filename <- checkExistantFile(outffname)
AllTerms <- summary(GOResultCC, pvalue = 2)
write.table(AllTerms, file = outffname, sep = "\t", row.names = FALSE, quote = FALSE)

#for output
newnames = c(newnames, outfname)
newfiles = c(newfiles, outffname)
newtypes = c(newtypes, 'txt')

#biological process
bpoutfname = "BP_RESULT"
bpoutfname = paste(title, bpoutfname, sep="_")
bpoutfname = paste(bpoutfname, "txt", sep=".")
bpoutffname = file.path(tdir, bpoutfname) #provides full file path

params <- new("GOHyperGParams", geneIds = selGenes, universeGeneIds = universeGenes,
annotation = anntIdentifier , ontology = "BP", pvalueCutoff = 0.05,conditional = FALSE, testDirection = "over")
GOResultBP <- hyperGTest(params)

#Filename <- checkExistantFile("GO-Mapping-resultBP.txt")
AllTerms <- summary(GOResultBP, pvalue = 2)
write.table(AllTerms, file = bpoutffname, sep = "\t", row.names = FALSE, quote = FALSE)

#for output
newnames = c(newnames, bpoutfname)
newfiles = c(newfiles, bpoutffname)
newtypes = c(newtypes, 'txt')

#Molecular Function
mfoutfname = "MF_RESULT"
mfoutfname = paste(title, mfoutfname, sep="_")
mfoutfname = paste(mfoutfname, "txt", sep=".")
mfoutffname = file.path(tdir, mfoutfname) #provides full file path

params <- new("GOHyperGParams", geneIds = selGenes, universeGeneIds = universeGenes,
annotation = anntIdentifier , ontology = "MF", pvalueCutoff = 0.05, conditional = FALSE, testDirection = "over")
GOResultMF <- hyperGTest(params)

#Saving the file of GO analysis, Molecular Function ontology
AllTerms <- summary(GOResultMF, pvalue = 2)
write.table(AllTerms, file = mfoutffname, sep = "\t", row.names = FALSE,quote = FALSE)

#for output
newnames = c(newnames, mfoutfname)
newfiles = c(newfiles, mfoutffname)
newtypes = c(newtypes, 'txt')

#for output

f1 = "Conduct_GO_using_David"
htmlname = paste(title, f1, sep="_")
htmlname = paste(htmlname, "html", sep=".")
htmlpath = file.path(tdir, htmlname) #provides full file path

newnames = c(newnames, htmlname)
newfiles = c(newfiles, htmlpath)
newtypes = c(newtypes, 'html')

sink()
return(list(newfiles=newfiles,newnames=newnames,newtypes=newtypes))
}
"""

mng = '### makenewgalaxy' # used by exec_after_process hook to parse out new file paths/names
tmp = '### tmp' # used by exec_after_process hook to clean up tmp directory

def getEntrezGeneIds(inputExprFilePath):

    #read expression index file
    f = file(inputExprFilePath, 'r')
    exprFilelines = f.readlines()
    f.close()

    #entrezIdPosn = 3
    entrezIdPosn = 4
    #drop header
    exprFilelines = exprFilelines[1:]
    entrezIdStr = ""
    count = 1
    #when reference genes are ALL
    for (line) in exprFilelines:
        if line.startswith('#'):
            continue
        columns = line.strip().split("\t")
        if(len(columns)>entrezIdPosn):
            value = columns[entrezIdPosn]
            if(entrezIdStr.find(value))==-1:
                if(count>200):
                    break
                entrezIdStr = entrezIdStr +","+ value
                count = count+1

    if(len(entrezIdStr)>0):
        entrezIdStr = entrezIdStr[1:]

    return entrezIdStr


def getDisplayGeneList(geneLinkList):
    gnDisplayStr = ""
    gnDisplayStrList = geneLinkList.split(',')

    for k in range(len(gnDisplayStrList)):
        entrezId = gnDisplayStrList[k]

        if(k%20==0):
            gnDisplayStr = gnDisplayStr+",<br>"+entrezId
        else:
            gnDisplayStr = gnDisplayStr+","+entrezId

    if(len(gnDisplayStr)>0):
        gnDisplayStr = gnDisplayStr[1:]

    return gnDisplayStr

def createGoGeneFormat(inFilename, outFilename):
  inFile = csv.reader(open(inFilename, 'rb'), delimiter="\t")
  outFile = file(outFilename, 'w')

  #geneCol = 0
  geneCol = 4
  hasGeneHead = False
  for rowNum, row in enumerate(inFile):
    if len(row) != 0 and row[0].startswith('#'):
      if not hasGeneHead:
        for colNum, header in enumerate(row):
          if header.lower() == 'na':
            geneCol = colNum
            outFile.write('Gene\n')
            hasGeneHead = True
            break
    else:
      if not hasGeneHead:
        outFile.write('Gene\n')
        hasGeneHead = True
      if len(row) == 0:
        continue

      if row[geneCol] != "":
        outFile.write(row[geneCol]+'\n')
  outFile.close()




def main():
    """called as
<command interpreter="python">
go_analysis.py  '$title' '$diff_expr_file' '$logmeta' '$diff_expr_file.dbkey'
</command>
    """
    nparm = 5
    appname = sys.argv[0]

    if (len(sys.argv) < nparm):
        print "Conduct GO 1.0.1"
        print '%s needs %d parameters - given %d as %s' % (appname,nparm,len(sys.argv),';'.join(sys.argv))
        print 'usage: %s title diff_expression_value out_file annotationKey' % (appname)
        sys.exit(1)
        # short command line error

    title = sys.argv[1].strip()
    inputExprFilePath = sys.argv[2].strip()
    logFileName = sys.argv[3].strip()
    #dbKey = sys.argv[4].strip()
    # dbKey has been replaced by compulsory annotation selection
    dbKey = sys.argv[5].strip()
    #outhtml = sys.argv[5].strip()
   
    gnStrList = getEntrezGeneIds(inputExprFilePath)

    #utilRscriptName=os.path.dirname(sys.argv[0]) + "/utils.R"
    utilRscriptName=os.path.dirname(sys.argv[0]) + "/utils.R"
    utilRscriptName = utilRscriptName.replace('\\', '/')

    print 'input anntkey:%s'% (dbKey)

    if dbKey == None or dbKey == "" or dbKey == "?":
        dbKey = "hgu133a"

    #create temp directory
    title = title.strip(' ').replace('  ',' ').replace(' ','_')
    
    #tdir = tempfile.mkdtemp(prefix=title)
    #tdir = tdir.replace('\\', '/')
    tdir = os.path.dirname(inputExprFilePath)

    print 'title:%s'% (title)
    print 'ExpressionFilePath:%s'% (inputExprFilePath)
    print 'logfilename:%s'% (logFileName)
    print 'anntkey:%s'% (dbKey)
    print 'tmp dir: %s'% (tdir)
    #print 'outhtml:%s'% (outhtml)

    # call R function & read return array
    tempExprFilePath = os.path.dirname(inputExprFilePath) + "/tempExprFile.txt"
    print 'tempExprFilePath: %s'% (tempExprFilePath)
    createGoGeneFormat(inputExprFilePath, tempExprFilePath)
    rprogname = 'rgo'
    if dbKey == 'fly.db0.db' : dbKey = 'fly.db0'
    rcall = "%s('%s','%s','%s','%s','%s','%s')" % (rprogname,tdir,title,tempExprFilePath,logFileName,dbKey,utilRscriptName)
    print '##rcall=',rcall

    #p = subprocess.Popen("R --vanilla", shell=True, executable="/bin/bash", stdin=subprocess.PIPE)
    p = subprocess.Popen("R --vanilla", shell=True, stdin=subprocess.PIPE)
    p.communicate(rprog + "\n" + rcall)
    makenew = getFilenames(tdir, title)

    newfiles = makenew[0]
    newnames = makenew[1]
    newtypes = makenew[2]

    #Write PDF file path,name and type to Log file for plot_code reading. i.e, exec_after_process()
    logf = file(logFileName,'a')
    s = 'R code and call\n%s\n\n\nCalled as\n%s\n' % (rprog, rcall)
    logf.write(s)
    try:
        outlist = ['%s\t%s\t%s\t%s' % (mng,newfiles[i],newnames[i],newtypes[i]) for i in range(len(newfiles))]
    except:
        outlist = ['%s\t%s\t%s\t%s' % (mng,newfiles,newnames,newtypes),] # only 1
    tmpString = '\n%s\t%s\n' % (tmp, tdir)
    logf.write(tmpString)
    logf.write('\n'.join(outlist))
    logf.write('\n')
    logf.close()


    for j in range(len(newfiles)):
        htmlfile = newfiles[j]
        ft = newtypes[j]
        if(ft=="html"):
            html = []
            html.append(galhtmlprefix)
            
            #html.append('<form method=post name="davidform" action="http://david.abcc.ncifcrf.gov/api.jsp">')
            #html.append('<input type=hidden name="type" value="ENTREZ_GENE_ID">')
            #html.append('<input type=hidden name="tool" value="summary">')
            #html.append('<input type=hidden name="ids" value="'+gnStrList+'">')
            #html.append('</form>')
            #html.append('<ul><li>')
            
            gnDisplayStr = getDisplayGeneList(gnStrList)

            html.append('<u><center><h3>Conduct GO using DAVID</h3><br></center></u>')
            html.append('<table width="95%" border="0"><tr><td><b>Target Entrez Gene Id List:</b></td></tr><tr><td>'+gnDisplayStr+'</td></tr></table><br><br>')

            html.append('<a href="http://david.abcc.ncifcrf.gov/api.jsp?type=ENTREZ_GENE_ID&ids='+gnStrList+'&tool=summary">Functional Annotation Summary Link</a><br><br>')
            html.append('<a href="http://david.abcc.ncifcrf.gov/api.jsp?type=ENTREZ_GENE_ID&ids='+gnStrList+'&tool=chartReport&annot=,GOTERM_BP_ALL,GOTERM_CC_ALL,GOTERM_MF_ALL,INTERPRO,PIR_SUPERFAMILY,SMART,BBID,BIOCARTA,KEGG_PATHWAY,COG_ONTOLOGY,SP_PIR_KEYWORDS,UP_SEQ_FEATURE,GENETIC_ASSOCIATION_DB_DISEASE,OMIM_DISEASE">Functional Annotation Chart Link</a><br><br>')
            html.append('<a href="http://david.abcc.ncifcrf.gov/api.jsp?type=ENTREZ_GENE_ID&ids='+gnStrList+'&tool=annotationReport&annot=,GOTERM_BP_ALL,GOTERM_CC_ALL,GOTERM_MF_ALL,INTERPRO,PIR_SUPERFAMILY,SMART,BBID,BIOCARTA,KEGG_PATHWAY,COG_ONTOLOGY,SP_PIR_KEYWORDS,UP_SEQ_FEATURE,GENETIC_ASSOCIATION_DB_DISEASE,OMIM_DISEASE">Functional Annotation Table Link</a><br><br>')
            html.append('<a href="http://david.abcc.ncifcrf.gov/api.jsp?type=ENTREZ_GENE_ID&ids='+gnStrList+'&tool=term2term&annot=,GOTERM_BP_ALL,GOTERM_CC_ALL,GOTERM_MF_ALL,INTERPRO,PIR_SUPERFAMILY,SMART,BBID,BIOCARTA,KEGG_PATHWAY,COG_ONTOLOGY,SP_PIR_KEYWORDS,UP_SEQ_FEATURE,GENETIC_ASSOCIATION_DB_DISEASE,OMIM_DISEASE">Functional Annotation Clustering Link</a><br><br>')
            html.append('<a href="http://david.abcc.ncifcrf.gov/api.jsp?type=ENTREZ_GENE_ID&ids='+gnStrList+'&tool=geneReportFull">Gene Full Report Link</a><br><br>')
            html.append('<a href="http://david.abcc.ncifcrf.gov/api.jsp?type=ENTREZ_GENE_ID&ids='+gnStrList+'&tool=geneReport">Gene Report Link</a><br><br>')
            html.append('<a href="http://david.abcc.ncifcrf.gov/api.jsp?type=ENTREZ_GENE_ID&ids='+gnStrList+'&tool=list">Show Gene List Names in Batch Link</a><br><br>')
            html.append('<a href="http://david.abcc.ncifcrf.gov/api.jsp?type=ENTREZ_GENE_ID&ids='+gnStrList+'&tool=gene2gene">Gene Functional Classfication Link</a><br><br>')
                      
            html.append("<br><i><b>Note:</b><P>Please be aware that maximum of 200 genes can be sent due to DAVID on-line internal limitation </P></i>")
            #html.append('</li></ul>')
            html.append(galhtmlpostfix)

            outf = file(htmlfile,'w')
            outf.write('\n'.join(html))
            outf.close()

def getFilenames(tdir, title):
    titles = ["_CC_RESULT.txt", "_BP_RESULT.txt", "_MF_RESULT.txt", "_Conduct_GO_using_David.html"]
    titles = [title + x for x in titles]
    paths = [os.path.join(tdir, x) for x in titles]
    return [
            paths,
            titles,
            ['txt', 'txt', 'txt', 'html']
           ]

if __name__ == "__main__":
    main()
