#################################################################################################
##
#################################################################################################
getSkipFromGPR <- function(file){
        con <- file(file,"r")
        twoLines <- readLines(con,n=2)
    close(con)
        if(substring(twoLines,1,3)[1]=="ATF"){
                skip <- substring(twoLines,1,2)
        }
        else{
                skip=c(0,-2)
        }
        as.numeric(skip[2])+2
}

###########################################################################################
## creates a new GAL file from the information from the gpr files
###########################################################################################
createGALFromGPR <- function(file,galfile="dummy.gal",colnames=c("Block","Column","Row","Name","ID")){
        table <- read.table(file,as.is=TRUE,comment.char="",check.names=FALSE,header=TRUE,skip=getSkipFromGPR(file),sep="\t")
        names <- colnames(table)
        result <- data.frame(matrix(ncol=5,nrow=length(table[,1])))
        colnames(result) <- colnames
        for(i in 1:length(result[1,])){
                result[,i] <- table[,colnames[i]]
        }
        #result[,1:5] <- table[,1:5]
        #colnames(result) <- names[1:5]
        write.table(file=galfile,result,sep="\t",row.name=FALSE,quote=FALSE)
}

plotHistogram <- function(R,G,r.col=2,g.col=3,breaks=NULL,log.transform=FALSE,add=FALSE,weights=NULL,return.breaks=FALSE,ylim=NULL,xlim=NULL,title=NULL,relative=FALSE,...){
        if(is.null(weights)){
                weights <- matrix(ncol=1,nrow=length(R))
                weights[,] <- 1
        }

        if(log.transform==TRUE){
                r <- log2(R)
                g <- log2(G)
                by <- 0.1
                toadd <- 0.1
        }
        else{
                r <- R
                g <- G
                by <- 50
                toadd <- 50
        }
        r <- r[weights==1]
        g <- g[weights==1]
        r <- r[!is.na(r)]
        g <- g[!is.na(g)]
        r <- r[!is.infinite(r)]
        g <- g[!is.infinite(g)]

        if(is.null(breaks)){
                breaks <- seq(min(0,(min(r,g))),(max(max(r),max(g))+toadd),by=by)
        }
        
        hr <- hist(r,plot=FALSE,breaks=breaks)
        hg <- hist(g,plot=FALSE,breaks=breaks)
        if(is.null(xlim)){
                xlim <- range(breaks)
        }

        r.counts <- hr$counts
        g.counts <- hg$counts
        ylabb <- "Counts"
        if(relative){
                r.counts <- r.counts / length(r)
                g.counts <- g.counts / length(g)
                ylabb <- "relative Counts"
                yy <- -0.001
        }
        else{
                yy <- -10
        }

        if(is.null(ylim)){
                ylim <- c(0,max(r.counts,g.counts))
        }
                
        if(!add){
                plot(hr$mids,r.counts,col=r.col,xlab="Intensity",ylab=ylabb,type="l",ylim=ylim,xlim=xlim,main=title,...)
        }
        else{
                points(hr$mids,r.counts,col=r.col,type="l",...)
        }
        points(hg$mids,g.counts,col=g.col,type="l",...)
        points(medianWithoutNA(r),yy,col=r.col,pch=3)
        points(medianWithoutNA(g),yy,col=g.col,pch=3)
        
        if(return.breaks){
                breaks
        }
#       hr <- density(r,na.rm=TRUE)
#       hg <- density(g,na.rm=TRUE)
#       xlim <- c(min(hr$x,hg$x),max(hr$x,hg$x))
#       ylim <- c(min(hr$y,hg$y),max(hr$y,hg$y))
#       plot(hr$x,hr$y,col=2,type="l",ylim=ylim,...)
#       points(hg$x,hg$y,col=3,type="l")
#       hr <- hist(r,breaks=breaks,plot=FALSE)
#       hg <- hist(g,breaks=breaks,plot=FALSE)
}


#############################################################################
##
#############################################################################
plotHistogramFromSlide <- function(lObject,slide=1,use.weights=TRUE,log.transform=FALSE,breaks=NULL,bg.lty=3,...){
        if(slide=="all"){
        
        
        
        }
        if(is.null(lObject$R)){
                R <- MAtoR(lObject$M[,slide],lObject$A[,slide])
                G <- MAtoG(lObject$M[,slide],lObject$A[,slide])
                Rb <- NULL
                Gb <- NULL
                draw.bg <- FALSE
                names <- colnames(lObject$M)
        }
        else{
                R <- lObject$R[,slide]
                G <- lObject$G[,slide]
                Rb <- lObject$Rb[,slide]
                Gb <- lObject$Gb[,slide]
                draw.bg <- !is.null(Rb)
                names <- colnames(lObject$R)
        }
        weights <- matrix(ncol=1,nrow=length(R))
        weights[,] <- 1
        if(is.null(lObject$weights)){
        }
        else{
                if(use.weights){
                        weights <- lObject$weights[,slide]
                }
        }
#       r <- R[weights==1]
#       g <- G[weights==1]
#       rb <- Rb[weights==1]
#       gb <- Gb[weights==1]
        
        breaks <- plotHistogram(R=R,G=G,log.transform=log.transform,weights=weights,breaks=breaks,return.breaks=TRUE,title=names[slide],...)
        if(draw.bg){
                plotHistogram(R=Rb,G=Gb,log.transform=log.transform,breaks=NULL,weights=weights,lty=bg.lty,add=TRUE,...)
        }
}




checkExistantFile <- function(file,path=".",no.subdir=TRUE){
        file <- gsub(pattern="_",replacement=".",file)
        file <- gsub(pattern="/",replacement="-",file)
        counter <- 0
        suff <- NULL
        ending <- NULL
        ffile <- file
        splitted <- unlist(strsplit(file,split="\\."))
        if(length(splitted)>1){
                ending <- paste(".",splitted[length(splitted)],sep="")
                ffile <- paste(splitted[1:(length(splitted)-1)],collapse=".")
                #cat(ffile,"_",ending,"\n")
        }
        while(TRUE){
                newfile <- paste(ffile,suff,ending,sep="")
                if(file.exists(newfile)){
                        #file does exist...
                        counter <- counter +1
                        suff <- paste("(",counter,")",sep="")
                }
                else{
                        break
                }
        }
        newfile
}




saveMAList <- function(object,filename,remove.flagged=TRUE,na.value="NA",subset=NULL, allvalues=TRUE){
        genes <- NULL
        weights <- NULL
        g.length <- 0
        col.names <- NULL
        if(class(object)=="MAList"){
                if(is.null(dim(object$M))){
                        object$M <- as.matrix(object$M)
                        object$A <- as.matrix(object$A)
                        if(!is.null(object$weights)){
                                object$weights <- as.matrix(object$weights)
                        }
                }
                if(is.null(subset)){
                        subset <- rep(TRUE,nrow(object$M))
                }
                M <- as.matrix(object$M[subset,])
                CN <- colnames(object$M)
                if(is.null(CN)){
                        CN <- object$targets[1,1]
                }
                colnames(M) <- CN
#               colnames(M) <- colnames(object$M)
                A <- as.matrix(object$A[subset,])
                colnames(A) <- CN
#               colnames(A) <- colnames(object$A)
                if(!is.null(object$weights)){
                        weights <- as.matrix(object$weights[subset,])
                }
                if(!is.null(object$genes)){
                        genes <- object$genes[subset,]
                        if(!is.null(genes)){
                                g.length <- ncol(genes)
                        }
                }
        }
        else if(class(object)=="EexprSet" | class(object)=="MadbSet"){
                if(object@type=="TwoColor"){
                        if(is.null(subset)){
                                subset <- rep(TRUE,nrow(exprs(object)))
                        }
                        M <- getM(object)
                        Colnames <- colnames(M)
                        M <- as.matrix(M[subset,])
                        colnames(M) <- Colnames
#                       M <- as.matrix(getM(object)[subset,])
#                       cat(class(M),colnames(M),"\n")
#                       A <- as.matrix(getA(object)[subset,])
                        A <- getA(object)
                        A <- as.matrix(A[subset,])
                        colnames(A) <- Colnames
                        #cat("1")
                        weights <- as.matrix(getWeights(object)[subset,])
                        #cat("2")
                        genes <- as.matrix(object@genes[subset,])
                        if(!is.null(genes)){
                                g.length <- ncol(genes)
                                if(g.length==1){
                                        colnames(genes) <- "ID"
                                }
                        }
                        #cat("g.length",g.length,"dim M",dim(M),"\n")
                        #cat("3")
                }
        }
        else if(class(object)=="RGList"){
                if(is.null(subset)){
                        subset <- rep(TRUE,nrow(object$R))
                }
                if(!.is.log(object$R)){
                        object$R <- log2(object$R)
                }
                if(!.is.log(object$G)){
                        object$G <- log2(object$G)
                }
                M <- object$R[subset,] - object$G[subset,]
                A <- 0.5*(object$R[subset,] + object$G[subset,])
                if(!is.null(object$weights)){
                        weights <- object$weights[subset,]
                }
                if(!is.null(object$genes)){
                        genes <- object$genes[subset,]
                        if(!is.null(genes)){
                                g.length <- ncol(genes)
                        }
                }
        }
        else{
                stop("Objects from the class ",class(object)," can not be handled by this function\n")
        }
        #cat("ncol=",g.length+(ncol(M)*2),"nrow=",nrow=nrow(M),"\n")
        Table <- matrix(ncol=(g.length+(ncol(M)*2)),nrow=nrow(M))
        if(!is.null(genes)){
                Table[,1:ncol(genes)] <- as.matrix(genes)
                col.names <- colnames(genes)
        }

        #cat(colnames(M),"\n")

        for(i in 1:ncol(M)){
                col.names <- c(col.names,paste(colnames(M)[i],"M"),paste(colnames(A)[i],"A"))
                Table[,((i*2)-1+g.length)] <- M[,i]
                if(!is.null(weights) & remove.flagged){
                        Table[weights[,i]==0,((i*2)-1+g.length)] <- NA
                }
                Table[is.na(Table[,((i*2)-1+g.length)]),((i*2)-1+g.length)] <- NA
                
                Table[,((i*2)+g.length)] <- A[,i]
                if(!is.null(weights) & remove.flagged){
                        Table[weights[,i]==0,((i*2)+g.length)] <- NA
                }
                Table[is.na(Table[,((i*2)+g.length)]),((i*2)+g.length)] <- NA
        }
        colnames(Table) <- col.names
        
        if(allvalues){
                ## calculate the red and green intensities from the M and A values...
                R <- MAtoR(M, A)
                G <- MAtoG(M, A)
                if(is.null(ncol(R))){
                        R <- matrix(R, ncol=1)
                        G <- matrix(G, ncol=1)
                }
                colnames(R) <- paste(colnames(M), "R")
                colnames(G) <- paste(colnames(M), "G")
                Table <- cbind(Table, R, G)
        }
        
        write.table(Table,file=filename,sep="\t",na=na.value,row.names=FALSE,quote=FALSE)
}

## closes all connections of gpr files...
closeMyConnections <- function(pattern="gpr"){
    All <- showConnections(TRUE)
    what <- grep(All[,"description"],pattern=pattern)
    what <- what-1
    for(theCon in what){
        close(getConnection(theCon))
    }
}

checkGPRs <- function(files=NULL,overwrite=TRUE){
        ##
        # check for correct axon format...
        requiredcolumns <- c("Block","Column","Row","Name","ID","F635 Median","F635 Mean","B635 Median","B635 Mean","F532 Median","F532 Mean","B532 Median","B532 Mean","Flags")
        for(i in 1:length(files)){
                con <- file(files[i],"r")
                twoLines <- readLines(con,n=2)
                close(con)
#        twoLines <- scan(file=files[i], what="character",nlines=2, quiet=TRUE)
#        closeMyConnections()
                if(substring(twoLines,1,3)[1]!="ATF"){
                        # is not in correct Axon file format
                        cat("The file",files[i],"is not in correct Axon file format!\n")
                        header <- read.table(files[i],header=TRUE,nrows=1,sep="\t",check.names=FALSE)
                        if(sum(colnames(header) %in% requiredcolumns)==length(requiredcolumns)){
                                BigFile <- readLines(con=files[i])
                                writeLines(c("ATF\t1.0",paste("1\t",ncol(header),sep=""),"Provider=CARMAweb"),con=files[i])
                                con <- file(files[i],"at")
                                writeLines(con=con,BigFile)
                flush(con)
                                close(con)
                                cat("Successfully adjusted file",files[i],".\n")
#                closeMyConnections()
                        }
                        else{
                                stop("The file",files[i],"is not in correct Axon file format! Was this file scanned with an Axon GenePix Scanner?\n")
                        }
                }
        }

        skips <- NULL
        for(i in 1:length(files)){
                skips <- c(skips,getSkipFromGPR(files[i]))
        }
        if(length(unique(skips))!=1){
                cat("GPR files have different number of header rows!\nAdjusting GPR files\n")
                for(i in 1:length(files)){
                        dummy <- read.table(files[i],skip=skips[i],header=TRUE,comment.char="",sep="\t",check.names=FALSE)
                        #cat(colnames(dummy),"\n")
                        write.table(dummy,file=files[i],sep="\t",quote=FALSE,row.names=FALSE)
                }
        checkGPRs(files)    # just re-run the check once again...
        }
}

createM <- function(object,red.channels,green.channels){
        M <- matrix(ncol=length(red.channels),nrow=nrow(exprs(object)))
        weights <- getWeights(object)
        col.names <- NULL
        for(i in 1:ncol(M)){
                r.idx <- myGrep(colnames(exprs(object)),pattern=red.channels[i],exact.match=TRUE)
                g.idx <- myGrep(colnames(exprs(object)),pattern=green.channels[i],exact.match=TRUE)
                M[,i] <- getM(object,r=r.idx,g=g.idx)
                M[weights[,r.idx]==0 | weights[,g.idx]==0,i] <- NA
                col.names <- c(col.names,paste(red.channels[i]," vs ",green.channels[i],sep=""))
        }
        colnames(M) <- col.names
        M
}


createA <- function(object,red.channels,green.channels){
        A <- matrix(ncol=length(red.channels),nrow=nrow(exprs(object)))
        weights <- getWeights(object)
        col.names <- NULL
        for(i in 1:ncol(A)){
                r.idx <- myGrep(colnames(exprs(object)),pattern=red.channels[i],exact.match=TRUE)
                g.idx <- myGrep(colnames(exprs(object)),pattern=green.channels[i],exact.match=TRUE)
                A[,i] <- getA(object,r=r.idx,g=g.idx)
                A[weights[,r.idx]==0 | weights[,g.idx]==0,i] <- NA
                col.names <- c(col.names,paste(red.channels[i]," vs ",green.channels[i],sep=""))
        }
        colnames(A) <- col.names
        A
}


##
# to remove all those genes that have flags in 
excludeFromTest <- function(data,classlabels,weights=NULL,nr.present=2){
        if(!is.null(weights)){
                if(sum(dim(data)==dim(weights))==2){
                        data <- .removeFlagged(data,weights)
                }
        }
        result <- apply(data,MARGIN=1,.dummyFun,cl=classlabels,np=nr.present)
        result
}

.removeFlagged <- function(data,weights,flagged.value=NA){
        for(i in 1:ncol(data)){
                data[weights[,i]==0,i] <- flagged.value
        }
        data
}

.dummyFun <- function(x,cl,np){
        res <- sum(!is.na(x[cl==0])) < np | sum(!is.na(x[cl==1])) < np
        res
}

tukey.biweight.mod <- function(x,c=5,epsilon=1e-4,na.rm=TRUE){
        x <- x[!is.na(x)]
        ret <- tukey.biweight(x,c=c,epsilon=epsilon)
}

if( !isGeneric("filterOnVariance") )
        setGeneric("filterOnVariance", function(object,variance=0,array.names=NULL,remove.flagged=TRUE,v=TRUE,...)
        standardGeneric("filterOnVariance"))

setMethod("filterOnVariance","MadbSet",
        function(object,variance=0,array.names=NULL,remove.flagged=TRUE,v=TRUE,...){
                if(is.null(array.names)){
                        array.names <- colnames(exprs(object))
                }
                E <- exprs(object)[,array.names]
                if(remove.flagged){
                        weights <- getWeights(object)[,array.names]
                        for(i in 1:ncol(weights)){
                                E[weights[,i]==0,i] <- NA
                        }
                }
                sds <- apply(E,MARGIN=1,sd,na.rm=TRUE)
                quant <- quantile(sds,variance,na.rm=TRUE)
                sel <- sds >= quant
                sel[is.na(sel)] <- FALSE
                if(v){
                        cat(sum(sel)," features out of ",length(sel)," have a sd bigger than ",quant,"\n")
                }
                ret <- object
                exprs(ret) <- exprs(object)[sel,]
                ret@weights <- getWeights(object)[sel,]
                if(nrow(exprs(object))==nrow(object@genes)){
                        if(ncol(object@genes)==1){
                                G <- matrix(object@genes[sel,],ncol=1)
                                rownames(G) <- rownames(exprs(ret))
                                colnames(G) <- colnames(object@genes)
                                ret@genes <- G
                        }
                        else{
                                ret@genes <- object@genes[sel,]
                        }
                }
                ret
})


##
# this function is a generalization of Robert Gentlemens GOHyperG funciton from the GOstats package.
# it allows to submit all EntrezGenes present on the used microarray.
# in principal this function performs the same steps as the GOHyperG function, nevertheless the results
# are sligtly different (maybe because of the different annotation packages used:
#  mapping to GO terms using the affymetrix probe sets -> lokuslink -> GO terms (GOstats) or mapping the locuslink ids directly
#  to GO terms using the GO package like it is done here.)
calculateGOHyperG <- function(x,allLLs=NULL,chip=NULL,what="MF",use.all=FALSE,v=TRUE){
        require(GO)
        require(GOstats)
        if (chip == "fly.db0.db"){chip = "fly.db0"}
        if(is.null(allLLs) & is.null(chip) & !use.all){
                stop("You have to submit either all EntrezGene IDs present on the array used, or the Affymetrix GeneChip used, or select use.all=TRUE to use all EntrezGenes available in GO!\n\n")
        }
        if(use.all){
                ##
                # just that the code below throws no exception...
                allLLs <- x
        }
        ##
        # if chip was submitted retrieve all LLs from the specific chip
        if(!is.null(chip)){
                require(chip,character.only=TRUE) || stop("need data package",chip)
                allLLs <- as.character(unlist(as.list(get(paste(chip,"LOCUSID",sep=""),mode="environment"))))
        }
        ##
        # prepare the x (myLLs) and allLLs. remove duplicates, NAs, use only those x that are present in allLLs
        # remove possible LL. or LL strings (old way of writing LocusLink IDs)
        x <- sub("LL.","",x)
        x <- sub("LL","",x)
        allLLs <- sub("LL.","",allLLs)
        allLLs <- sub("LL","",allLLs)
        allLLs <- unique(allLLs)
        allLLs <- allLLs[!is.na(allLLs) & allLLs!=""]
        x <- unique(x)
        x <- x[!is.na(x) & x!=""]
        ##
        # just interested in those x that are present in allLLs
        x <- allLLs[match(x,allLLs)]
        x <- x[!is.na(x) & x!=""]
        
        ##
        # do some checking...
        if(length(x)==0 || is.na(x)){
                stop("None of the submitted interesting genes are present in the list of all EntrezGene (LocusLink) IDs submitted!\nPossible causes:\n-wrong list of genes (wrong IDs) submitted\n-EntrezGene IDs of the interesting genes and list of EntrezGene IDs of all genes on the array are not present on the same microarray\n")
        }
        
        
        
        ##
        # ok, now map x to GO, retrieve all GO terms where one of the x EntrezGene IDs is present (using GOALLLOCUSID from the GO package)
        GOLL <- as.list(get("GOALLLOCUSID",mode="environment"))     ## warning, GOALLLOCUSID is deprecated in newer versions, use GOALLENTREZID
        GOLL <- GOLL[!is.na(GOLL)] # just removing all the GO ids that are not mapped to any LL
        PresentGO <- sapply(GOLL,function(z){
                if(is.na(z) || length(z)==0)
                        return(FALSE)
                any(x %in% z)
                }
        )
        ##
        # another checking...
        if(max(PresentGO)==0){
                stop("The submitted list of EntrezGene IDs of the interesting genes can not be mapped to any of the GO terms.\nMaybe the list contained no EntrezGene (LocusLink) IDs?\n")
        }
        if(v){
                cat("The",length(x),"submitted EntrezGene IDs map to",sum(PresentGO),"GO terms in all tree categories.\n")
        }
        GOLL <- GOLL[PresentGO]
        ##
        # next restrict the GO terms to the selected category (specified by "what")
        goCat <- unlist(getGOOntology(names(GOLL)))
        goodGO <- goCat == what
        goodGO[is.na(goodGO)] <- FALSE
        if(v){
                cat(sum(goodGO),"of them are in the GO category we are interested in.\n")
        }
        mm <- match(names(goCat),names(GOLL))
        mm <- mm[goodGO]
        GOLL <- GOLL[mm]
        
        ##
        # next we have to remove all the EntrezGene IDs from each GO term, that was not present on the microarray used.
        if(!use.all){
                myGOLL <- sapply(GOLL,function(z){
                        if(length(z)==0 || is.na(z))
                                return(NA)
                        mylls <- unique(z[match(allLLs,z)])
                        mylls <- mylls[!is.na(mylls)]
                        mylls
                }
                )
                ##
                # removing all those GO terms that have no entries left (just in the case there are some...)
                empties <- sapply(myGOLL,function(z){length(z)==1 && is.na(z)})
                myGOLL <- myGOLL[!empties]
                GOLL <- GOLL[!empties]
        }
        else{
                ##
                # if use.all was selected, the proportion of the submitted EntrezGene IDs per GO term is compared to all
                # EntrezGenes mapped to that GO term.
                myGOLL <- GOLL
        }
        

        ##
        # the rest is the same as in the GOHyperG function...
        cLLs <- unique(unlist(myGOLL))
        nLL <- length(cLLs)                    # nr of unique LL IDs present in all GO terms
        goCounts <- sapply(myGOLL,length)      # nr of LL IDs present on the array per GO term
        #allCounts <- sapply(GOLL,length)       # nr of all LL IDs per GO term
        ##
        # some LL ids are represented more than once in a GO term...
        uniqueAllCounts <- sapply(GOLL,
                function(z){
                        if(length(z)==0 || is.na(z))
                                return(0)
                        lls <- unique(z)
                        lls <- lls[!is.na(lls)]
                        return(length(lls))
                })
        whGood <- x[x %in% cLLs]               # those LLs from the interesting ones that are present in the GO terms selected
        nInt <- length(whGood)
        if(nInt==0){
                warning("No interesting genes left")
        }
        
        MappingGOLL <- lapply(myGOLL,function(z){whGood[whGood %in% z]})
        useCts <- sapply(myGOLL,function(z){sum(whGood %in% z)})   # nr of interesting LLs per GO term
        ##
        # calculate the hypergeometric p values...
        pvs <- phyper(useCts-1,nInt,nLL-nInt,goCounts,lower.tail=FALSE)
        ord <- order(pvs)

        return(
                list(pvalues=pvs[ord],goCounts=goCounts[ord],intCounts=useCts[ord],allCounts=uniqueAllCounts[ord],GOLL=myGOLL[ord],intMapping=MappingGOLL[ord],numLL=nLL,numInt=nInt,intLLs=x)
        )
        
}


################
## get EntrezGenes for the specified GO term
getEntrezIDsFromGO <- function(GO, package="GO", entrezgenes){
        match.arg(package, c("GO", "humanLLMappings", "mouseLLMappings", "ratLLMappings"))
        require(package, character.only=TRUE)
        if(package=="GO"){
                ALLIDs <- mget(GO, env=GOALLENTREZID, ifnotfound=NA)
                if(!missing(entrezgenes)){
                        ## subset to the specified entrezgene ids.
                        ALLIDs <- lapply(ALLIDs, function(x){ x[ x %in% entrezgenes] })
                }
                return(ALLIDs)
        }
        if(package=="humanLLMappings"){
                ALLIDs <- mget(GO, env=humanLLMappingsGO2LL, ifnotfound=NA)
                if(!missing(entrezgenes)){
                        ## subset to the specified entrezgene ids.
                        ALLIDs <- lapply(ALLIDs, function(x){ x[ x %in% entrezgenes] })
                }
                return(ALLIDs)          
        }
        if(package=="mouseLLMappings"){
                ALLIDs <- mget(GO, env=mouseLLMappingsGO2LL, ifnotfound=NA)
                if(!missing(entrezgenes)){
                        ## subset to the specified entrezgene ids.
                        ALLIDs <- lapply(ALLIDs, function(x){ x[ x %in% entrezgenes] })
                }
                return(ALLIDs)          
        }
        if(package=="ratLLMappings"){
                ALLIDs <- mget(GO, env=ratLLMappingsGO2LL, ifnotfound=NA)
                if(!missing(entrezgenes)){
                        ## subset to the specified entrezgene ids.
                        ALLIDs <- lapply(ALLIDs, function(x){ x[ x %in% entrezgenes] })
                }
                return(ALLIDs)          
        }
}


#getLLGOmapping <- function(x,what="MF"){
#       gos <- mget(x, env=GOLOCUSID2GO,ifnotfound=NA)
#       if(length(gos)==1)
#               bd = is.na(gos[[1]])
#       else
#               bd = is.na(gos)
#       gos <- gos[!bd]
#       mymapping <- lapply(gos,getOntology,ontology=what)
#       emptygos <- sapply(mymapping,function(x) length(x)==0 )
#       mymapping <- mymapping[!emptygos]
#       mymapping
#}

#getGOFromLL <- function(x,Ontology="MF"){
#       require(GO) || stop("no GO library")
#       match.arg(Ontology, c("MF", "BP", "CC"))
#       goterms <- mget(x, env = GOLOCUSID2GO, ifnotfound=NA)
#       if(length(goterms)==1){
#               goterms <- goterms[!is.na(goterms[[1]])]
#       }
#       else{
#               goterms <- goterms[!is.na(goterms)]
#       }
#       goterms <- lapply(goterms, function(x) x[sapply(x,
#                      function(x) {if (is.na(x$Ontology) )
#                                          return(FALSE)
#                                      else
#                                          x$Ontology == Ontology})])
#       goids <- NULL
#       if(length(goterms)==0){
#               stop("The submitted list of EntrezGene IDs of the interesting genes can not be mapped to any of the GO terms.\nMaybe the list contained no EntrezGene (LocusLink) IDs?\n")
#       }
#       for(i in 1:length(goterms)){
#               goids <- c(goids,unique(unlist(sapply(goterms[[i]], function(x) x$GOID))))
#       }
#       goids
#}
#

##
# check if package is available
isAnnotationAvailable <- function(chip,develOK=FALSE){
        ok = FALSE
        if (chip == "fly.db0.db"){chip  = "fly.db0"}
       if(!require("reposTools",character.only=TRUE)){
               warning("reposTools package is not installed.\n")
               if(!require(chip,character.only=TRUE)){
                       ok = FALSE
               }
               else{
                       ok = TRUE
               }
       }
       else{
                if(!require(chip,character.only=TRUE)){
             source("http://www.bioconductor.org/getBioC.R")
                      library(reposTools)
                      # have to install the package
                      tryCatch(install.packages2(pkg=chip,develOK=develOK,getAllDeps=TRUE),error=return(FALSE))
                      try(getBioC(chip))
                      if(!require(chip,character.only=TRUE)){
                          ok=FALSE
                      }
                       else{
                           ok=TRUE
                       }
                }
                else{
                        ok = TRUE
                }
       }
        ok
}

saveSAMResults <- function(SAM,table,delta,filename){
        dummy <- print(SAM)
        for(i in 1:nrow(dummy)){
                saveSAMResult(SAM,table,delta=dummy[i,"Delta"],filename=filename)
        }
}

saveSAMResult <- function(SAM, table,delta,filename){
        sum.SAM <- (summary(SAM,delta))@mat.sig
        if(!is.null(sum.SAM) & nrow(sum.SAM)>0){
                table.sub <- table[sum.SAM[,"Row"],]
                if(class(table)=="aafTable"){
                        for(i in 2:(ncol(sum.SAM))){
                                table.sub <- merge(table.sub,aafTable(sum.SAM[,i],colnames=colnames(sum.SAM)[i]))
                        }
                        Filename <- checkExistantFile(paste(filename,"-delta-",delta,".txt",sep=""))
                        saveText(table.sub,filename=Filename)
                }
                else{
                        if(nrow(sum.SAM)==1){
                                table.sub <- matrix(table.sub,nrow=1)
                                colnames(table.sub) <- colnames(table)
                        }
                        table.sub <- cbind(as.matrix(table.sub),as.matrix(sum.SAM[,2:ncol(sum.SAM)]))
                        Filename <- checkExistantFile(paste(filename,"-delta-",delta,".txt",sep=""))
                        write.table(table.sub,file=Filename,sep="\t",quote=FALSE,row.names=FALSE)
                }
                cat(nrow(table.sub),"significant genes (delta=",delta,") were saved to the file:",Filename,"\n")
        }
        else{
                cat("No significant genes with a delta of",delta,"\n")
        }
}

###
# fills missing spots
fillMissingSpots <- function(object){
        ret <- object
        if(class(ret)=="RGList"){
                if(!is.null(ret$printer)){
                        layout <- object$printer
                        nr.row <- layout$ngrid.r * layout$ngrid.c * layout$nspot.r * layout$nspot.c
                        missing <- (nr.row) - nrow(object$R)
                        if(missing!=0){
                                cat("Have to fill in",missing,"missing spots\n")
                                index <- as.numeric(apply(object$genes,MARGIN=1,FUN=.getIndexFromBCR,nspot.r=layout$nspot.r,nspot.c=layout$nspot.c))
                                R <- matrix(ncol=ncol(object$R),nrow=nr.row)
                                colnames(R) <- colnames(object$R)
                                rownames(R) <- 1:nrow(R)

                                ## check if length of index and nrow of R are the same (<- problem with custom GAL files...)
                                if(length(index)!=nrow(object$R)){
                                        stop("ERROR: number of data values does not match number of annotations!\nPlease try again without specifying a GAL file!")
                                }
                                R[index,] <- object$R
                                ret$R <- R
                                R[index,] <- object$G
                                ret$G <- R
                                if(!is.null(object$Rb)){
                                        R[index,] <- object$Rb
                                        ret$Rb <- R
                                        R[index,] <- object$Gb
                                        ret$Gb <- R
                                }
                                if(!is.null(object$weights)){
                                        R[,] <- 0
                                        R[index,] <- object$weights
                                        ret$weights <- R
                                }
                                else{
                                        R[,] <- 0
                                        R[index,] <- 1
                                        ret$weights <- R
                                }
                                rm(R)
                                g <- gc()
                                Genes = data.frame(matrix(ncol=ncol(object$genes),nrow=nr.row))
                                colnames(Genes) <- colnames(object$genes)
                                rownames(Genes) <- 1:nrow(Genes)
#                               cat("nrow Genes",nrow(Genes),"ncol Genes",ncol(Genes),"length index",length(index),"\n")
                                Genes[index,1:ncol(Genes)] = object$genes[,1:ncol(Genes)]
#                               cat("genes filled\n","ncol= ",ncol(Genes)," nrow: ",nrow(Genes),"\n")
                                ret$genes = NULL
                                ret$genes = Genes
#                               ret$genes = Genes
#                               cat("genes set\n")
                                rm(Genes)
                                g <- gc()
                        }
                }
        }
        ret
}

.getIndexFromBCR <- function(x,nspot.c,nspot.r){
        ((as.numeric(x["Block"])-1)*nspot.c*nspot.r) + ((as.numeric(x["Row"])-1)*nspot.c) + as.numeric(x["Column"])
}

if( !isGeneric("drawBoxplot") )
        setGeneric("drawBoxplot", function(x,new.plot=TRUE,xlim=NULL,at=0,exclude.flagged=FALSE,...)
        standardGeneric("drawBoxplot"))

setMethod("drawBoxplot","MadbSet",
        function(x,new.plot=TRUE,xlim=NULL,at=0,exclude.flagged=FALSE,...){
                if(exclude.flagged){
                        W <- getWeights(x)
                        for(i in 1:ncol(exprs(x))){
                                exprs(x)[W[,i]==0,i] <- NA
                        }
                }
                boxplots(x=exprs(x),new.plot=new.plot,xlim=xlim,at=at,...)
        }
)

boxplots <- function(x,new.plot=TRUE,xlim=NULL,at=0,col=NULL,log2.transform=FALSE,ylim=NULL,...){
        data <- x
#       if(class(x)=="exprSet" | class(x)=="EexprSet"){
#               data <- exprs(x)
#       }
        if(log2.transform){
                data <- log2(data)
        }
        if(is.null(ncol(data))){
                data <- as.matrix(data)
        }
        ## remove those nasty infinite values
        for(i in 1:ncol(data)){
                data[is.infinite(data[,i]),i] <- NA
        }

        CN <- colnames(data)
        if(is.null(CN)){
                CN <- 1:ncol(data)
        }
        if(is.null(col)){
                col=rep(0,ncol(data))
        }
        if(length(col)!=ncol(data)){
                col=rep(col[1],ncol(data))
        }
        if(new.plot){
        if(is.null(xlim)){
                xlim=c(0.5,(ncol(data)+0.5))
        }
        par(omi=0.7*c(min(2,max(strwidth(CN,units="i"))),0,0,0),cex.axis=0.7)
        if(is.null(ylim)){
                ylim=c(min(data,na.rm=TRUE),max(data,na.rm=TRUE))
        }
        plot(1,1,pch=NA,ylab="expression",xaxt="n",ylim=ylim,xlim=xlim,xlab=NA)
        }
        axis(side=1,labels=CN,at=(1+at):(ncol(data)+at),las=3)
        for(i in 1:ncol(data)){
                boxplot(data[,i],at=(i+at),add=TRUE,col=col[i],range=0,...)
        }
}


plotRankedVariance <- function(data,subset=NULL,cut){
        if(!is.null(subset)){
                data <- data[,subset]
        }
        if(!.is.log(data)){
                data <- log2(data)
        }
        SD <- apply(data,MARGIN=1,sd,na.rm=TRUE)
        CUT <- quantile(SD,cut,na.rm=TRUE)
        plot(sort(SD),pch=16,cex=0.33,ylab="sorted variance")
        abline(h=CUT,col=2)
        text((length(SD)/20),(CUT+0.2),label=format(CUT,digits=4),col=2)
}

plotSortedPValues <- function(x,top=NULL,supported.cols=c("rawp","adjp","minP","maxT","BH","BY","Bonferroni","Hochberg","Holm","SidakSS","SidakSD"),...){
        if(is.null(nrow(x))){
                stop("Please submit a matrix, data.frame or instance of aafTable!\n")
        }
        if(!is.null(top)){
                top <- min(c(nrow(x),top))
        }
        else{
                top <- nrow(x)
        }
        Columns <- supported.cols[supported.cols %in% colnames(x)]
        if(length(Columns)==0){
                stop("The submitted matrix must contain at least one column with p-values that is named according to one of the supported column names!")
        }
        if(class(x)=="aafTable"){
                require(annaffy)
                #have to do some transformations...
                x.new <- matrix(ncol=length(Columns),nrow=nrow(x))
                colnames(x.new) <- Columns
                for(i in Columns){
                        x.new[,i] <- as.numeric(getText(x[[i]]))
                }
                x <- x.new
        }
        #ok, now i can start drawing...
        x <- apply(x[,Columns],MARGIN=2,sort,na.last=TRUE)
        Max <- max(x[1:top,],na.rm=TRUE)

        plot(x[1:top,Columns[1]],xlab="number of significant genes",ylab="sorted p-value",col=1,ylim=c(0,Max),main=paste("Top",top,"sorted p-values"),...)
        if(length(Columns)>1){
                for(i in 2:length(Columns)){
                        points(x[1:top,Columns[i]],col=i,...)
                }
        }
        legend(top-top/4,Max,Columns,col=1:length(Columns),pch=16)
}

.is.log <- function(data,maxi=100){
        data <- data[!is.na(data)]
        data <- data[!is.infinite(data)]
#       cat("Try to check wheter or not the values are in log scale...")
        if(max(data>maxi)){
                log.values = FALSE
#               cat("signals are not in log2 scale, transforming...\n")
        }
        else{
                log.values = TRUE
#               cat("signals are already in log2 scale \n")
        }
        log.values
}

toNumeric <- function(x,Sub=","){
        if(!is.null(dim(x))){
                RN <- rownames(x)
                if(class(x[1,])=="numeric"){
                        return(x)
                }
                if(!is.null(Sub)){
                        for(s in Sub){
                                x <- apply(x,MARGIN=2,sub,pattern=s,replacement="")
                        }
                }
                x <- tryCatch(apply(x,MARGIN=2,as.numeric),error=function(error){stop("ERROR: Can not perform calculations on non-numeric data! Please submit only numeric values.")},warning=function(warning){stop("ERROR: Can not perform calculations on non-numeric data! Please submit only numeric values.")})
                rownames(x) <- RN
                x
        }
}

####
## function similar to PWamat function from the annotate package,
## generates a matrix with rows unique LL ids, columns KEGG pathways
## containing 1 if gene can be associated with the pathway, 0 otherwise.
makeLL2KEGGMappingMatrix <- function(x,nrGenes=0,verbose=TRUE){
        if(is.null(x)){
                stop("No EntrezGene (LocusLink) IDs have been submitted!\n")
        }
        require(KEGG)
        x <- unique(as.character(x))
        paths <- mget(x,KEGGEXTID2PATHID,ifnotfound=NA)
        not_found <- sapply(paths,function(x){if(length(x)==1 && is.na(x)) TRUE else FALSE})
        paths <- paths[!not_found]
        if(verbose){
                cat(length(paths),"EntrezGene IDs from the inputlist can be associated with at least one KEGG pathway.\n",sum(not_found),"IDs can not be associated with KEGG pathways.\n")
        }
        unique_paths <- unique(unlist(paths))
        Test <- sapply(paths,function(x){
                mtch <- match(x,unique_paths)
                res <- rep(0,length(unique_paths))
                res[mtch] <- 1
                res
        })
        Test <- t(Test)
        rownames(Test) <- x
        colnames(Test) <- unique_paths
        if(nrGenes>0){
                to_exclude <- colSums(Test) < nrGenes
                if(verbose){
                        cat(sum(to_exclude),"Pathways are excluded, because less than",nrGenes,"genes can be mapped to them.\n")
                }
                Test <- Test[,!to_exclude]
        }
        return(Test)
}
