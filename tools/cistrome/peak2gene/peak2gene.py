#!/usr/bin/env python
"""
#ver_0.03
#1. output two file: reflist and detail.
#2. use flag_line, flag_skip to accelerate Efficiency of script.
#3. output genes up/down/in separately.
#4. input file support MACS result: .bed/.xls

#ver_1.00
#add option --fo
#change the cols of result file

#ver_1.01
#always calc distance from peak center to gene TSS

#ver_1.02
#output gene can be symbol or refseq
#delete option --fo, delete xls file support
#output 2 file: annotation for each gene, and annotation for peaks.
#delete column "pos" in output, use +/- for "TSS2pCenter" column.
#use "|" as separate instead of ","

"""
import sys, os, time, re
import sqlite3
from optparse import OptionParser

# print current time and information on screen
def Info(infoStr):
    print "[%s] %s" %(time.strftime('%H:%M:%S'), infoStr)

def prepare_optparser():
    """Prepare optparser object. New options will be added in this
    function first.
    """
    usage = "usage: %prog <-t FILE -o FILE -d NUMBER -g GENOME_FILE> [options]"
    description = "Input a peak file, and It will search each peak on UCSC to get the refGenes near the peak summit/center."

    optparser = OptionParser(version="%prog v1.02", description=description, usage=usage, add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")
    optparser.add_option("-t","--treat",dest="peakFile",type="string",
                         help="Input the MACS's result peak file(.xls/.bed), it will recognize it by extension.")
    optparser.add_option("-n","--name",dest="name",type="string",
                         help="this argument is used to name the result file, it will output two file 'name_peaks_annotation.txt', 'name_gene_annotation.txt'.")
    optparser.add_option("--op",dest="outpos",type="string",
                        help="""select which type of genes need to output.\n
                        'up' for genes upstream to peak summit,\n
                        'down' for genes downstream to peak summit,\n
                        'all' for both 'up' and 'down'""", default="all")
    optparser.add_option("--symbol", dest="symbol", action="store_true", 
                         help="output gene symbol instead of refseq name.", default=False)
    optparser.add_option("-d","--distance", dest="distance", type="int",
                         help="Set a number which unit is 'base'. It will get the refGenes in n bases from peak center. default:10000", default=10000)
    optparser.add_option("-g","--genome",dest="genome",type="string",
                         help="Select a genome file (sqlite3 file) to search refGenes.")
    optparser.add_option("-e","--entrez",dest="translation",type="string",
                         help="Select a translation file to map refGenes to Entrez Ids. If no file is proveded no thranslations are available")

    return optparser

def opt_validate(optparser):
    """Validate options from a OptParser object.

    Return: Validated options object.
    """
    (options,args) = optparser.parse_args()
    #quoted bug fix
    if options.peakFile and '"' in options.peakFile:
        options.peakFile = options.peakFile.replace('"','')
    if options.translation and '"' in options.translation:
        options.translation = options.translation.replace('"','')
    if options.genome and '"' in options.genome:
        options.genome = options.genome.replace('"','')
        
    if not options.peakFile and not options.genome and not options.name:
        optparser.print_help()
        sys.exit(1)
    if not os.path.isfile(options.peakFile):
        Info('Wanring: cannot find peak file, a tab-peak file must be given through -t (--treat).')
        sys.exit(1)

    if not os.path.isfile(options.genome):
        Info("Warning: Genome file not found! An annotation file must be given through -g (--genome).")
        sys.exit(1)
        
    if not os.path.isfile(options.translation):
        Info("Warning: Translation file not found! A translation file must be given through -t (--translations). No entrez ID translations will be availalbe.")
        options.translation = ""

    if not options.name:
        options.name = os.path.splitext(options.peakFile)[0] + "result"

    if options.outpos not in ("up", "down"):
        Info("use default: find up+down")
        options.outpos = ("up", "down")
    else:
        options.outpos = (options.outpos,)

    # print arguments
    Info("Argument List: ")
    Info("Name = " + options.name)
    Info("peak file = " + options.peakFile)
    Info("gene pos to peak = " + ",".join(options.outpos))
    Info("distance = %d bp" %options.distance)
    Info("genome = %s" %options.genome)
    Info("Entrez translation = %s" %options.translation)
    Info("Output symbol as gene name = %s" %str(options.symbol))
    print
    return options

class PCA:
    # connect to sqlite, select genome.
    def __init__(self, options):
        self.genome = options.genome
        self.translation = options.translation
        self.getGeneType = options.outpos # up, down
        self.symbol = options.symbol
        self.peakFile = options.peakFile
        self.db = sqlite3.connect(self.genome)
        self.c = self.db.cursor()
        self.opts_string = "# Argument List:\n" +\
                           "# Name = %s\n" %options.name +\
                           "# peak file = %s\n" %options.peakFile +\
                           "# gene pos to peak = %s\n" %(",".join(options.outpos),) +\
                           "# distance = %d bp\n" %options.distance +\
                           "# genome = %s\n" %options.genome +\
                           "# Entrez translations = %s\n" %options.translation +\
                           "# Output symbol as gene name = %s\n" %str(options.symbol)
        self.peakList = []
        
        self.geneEntresIds = {}

    # import peakFile
    def readfile(self): #reads the file and returns a vector: each element is a bed_row. 
        peakf = open(self.peakFile)
        count = 0
        for line in peakf:
            if line.startswith('track') or line.startswith('#') or line.startswith('browser') or not line.strip(): #skip "#" lines and empty lines
                continue
            line = line.split() #.bed-> 0:chrom 1:pStart 2:pEnd 3:peakName 4:-10*log10(pvalue)
            if len(line) < 5:
                line.extend(['NA', '0'][len(line)-5:])
            self.peakList.append(line)
            count += 1
        peakf.close()
        self.peakList.sort()
        Info("Read file<%s> OK! <%d> peaks." %(self.peakFile, count-1))
        
    # import entrez translations file
    def readTranslationsFile(self):
        if os.path.isfile(self.translation):
            entrezf = open(self.translation)
            count = 0
            for line in entrezf:
                if line.startswith('track') or line.startswith('#') or line.startswith('browser') or not line.strip(): #skip "#" lines and empty lines
                    continue
                line = line.split() #Format: tax_id GeneID status RNA_nucleotide_accession.version RNA_nucleotide_gi protein_accession.version protein_gi genomic_nucleotide_accession.version genomic_nucleotide_gi start_position_on_the_genomic_accession end_position_on_the_genomic_accession orientation assembly (tab is used as a separator, pound sign - start of a comment) 
                if len(line) < 13:
                    continue
                self.geneEntresIds[line[3][0:-2]] = line[1]
                count += 1
            entrezf.close()
            #self.peakList.sort()
            Info("Read translations file<%s> OK! <%d> IDS." %(self.translation, count-1))
    #prepares string out of a list for writing. "-" if the list is empty
    def prepareJoinedString(self, list):
        if len(list) == 0:
            return '-'
        else:
            return '|'.join(list)
    # get refgenes +/-nb from ChIP center
    def GetChipSummitRefgene(self, distance):
        # get data from sqlite file
        sql = "select chrom, name, txStart, txEnd, strand, name2 from GeneTable order by chrom;"
        self.c.execute(sql)
        geneInfo = self.c.fetchall() # (0:chrom, 1:name, 2:txStart, 3:txEnd, 4:strand)
        geneInfo.sort()
        geneInfo.append(("NA",)) # add a last line, otherwise it will miss the last hit

        # add result file head line
        self.annotated = [["#chrom", "pStart", "pEnd", "pName", "pScore", "NA", "gene", "strand", "TSS2pCenter"]]
        self.geneanno = [["#chrom", "gTSS", "gTTS", "gene", "NA", "strand", "pName", "TSS2pCenter"]]
        self.genelist = ["genelist"]

        count = 0
        flag_line = 0 # mark the position of genelist having read.
        for peak in self.peakList:
            #if not re.findall("chr\d+", peak[0]):
            #    continue
            geneList = []
            strandList = []
            disList = []
            chrom, pSummit, pName, pScore = peak[0], (int(peak[1])+int(peak[2]))/2, peak[3], peak[4]

            flag_skip = 0 # use to optimize script.
            for gene in geneInfo[flag_line:]:
                if peak[0] == gene[0]:
                    flag_skip = 1
                    if gene[4] == "+" and pSummit-distance < gene[2] < pSummit+distance: # it's the gene near peak we need
                        dist = gene[2]-pSummit
                    elif gene[4] == "-" and pSummit-distance < gene[3] < pSummit+distance: # gene at -
                        dist = pSummit-gene[3]
                    else:
                        continue

                    if ("up" in self.getGeneType and dist <= 0) or ("down" in self.getGeneType and dist >= 0):
                        if self.symbol:
                            geneList.append(gene[5])
                        else:
                            geneList.append(gene[1])
                        strandList.append(gene[4])
                        disList.append(str(dist))
                        if geneList[-1] not in self.genelist:
                            self.genelist.append(geneList[-1])
                            entrez = '0'
                            #print self.geneEntresIds[geneList[-1]]
                            if geneList[-1] in self.geneEntresIds:
                                
                                entrez = self.geneEntresIds[geneList[-1]]
                            self.geneanno.append([gene[0], str(gene[2]), str(gene[3]), geneList[-1], entrez, gene[4], pName, str(dist)])
                        else:
                            ig = self.genelist.index(geneList[-1])
                            self.geneanno[ig][6] +='|'+pName
                            self.geneanno[ig][7] +='|'+str(dist)    
                else:
                    if flag_skip == 0:
                        flag_line += 1
                    elif flag_skip == 1:
                        self.annotated.append([chrom, peak[1], peak[2], pName, pScore, '0', self.prepareJoinedString(geneList), self.prepareJoinedString(strandList), self.prepareJoinedString(disList)])
                        count += 1
                        if count%500 == 0:
                            Info("process peaks <%d>" %count)
                        break # ->do next peak

    # write peaks and refgenes to file
    def Output2File(self, name):
        outf = open("%s_peaks_annotation.txt"%name, "w")
        outf.write(self.opts_string)
        for each in self.annotated:
            print each
            outf.write("\t".join(each)+"\n")
        outf.close()

        outf = open("%s_gene_annotation.txt"%name, "w")
        outf.write(self.opts_string)
        for each in self.geneanno:
            #print each
            outf.write("\t".join(each)+"\n")
        outf.close()

        Info("Please Cheak outfile: <%s_peaks_annotation.txt> <%s_gene_annotation.txt>" %(name, name))
   
        

def main():
    opts=opt_validate(prepare_optparser())
    g = PCA(opts)
    g.readfile()
    g.readTranslationsFile()
    g.GetChipSummitRefgene(opts.distance)
    g.Output2File(opts.name)

if __name__ == "__main__":
    main()