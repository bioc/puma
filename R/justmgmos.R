justmgMOS <- function(
  ...
, filenames=character(0)
, widget=getOption("BioC")$affy$use.widgets
, compress=getOption("BioC")$affy$compress.cel
, celfile.path=getwd()
, sampleNames=NULL
, phenoData=NULL
, description=NULL
, notes=""
, background=TRUE
, gsnorm=c("median", "none", "mean", "meanlog")
, savepar=FALSE
, eps=1.0e-6
)
{

  l <- AllButCelsForReadAffy(..., filenames=filenames,
                             widget=widget,
                             celfile.path=celfile.path,
                             sampleNames=sampleNames,
                             phenoData=phenoData,
                             description=description)


  ##and now we are ready to read cel files
 ret<- just.mgmos(filenames=l$filenames,
                  phenoData=l$phenoData,
                  description=l$description,
                  notes=notes,
                  compress=compress,
                  background=background,
                  gsnorm=gsnorm,
                  savepar=savepar,
                  eps=eps)
  ##sampleNames(ret) <- l$sampleNames
  return(ret)

}

just.mgmos <- function(
  ...
, filenames=character(0)
, phenoData=new("AnnotatedDataFrame")
, description=NULL
, notes=""
, compress=getOption("BioC")$affy$compress.cel
, background=TRUE
, gsnorm=c("median", "none", "mean", "meanlog")
, savepar=FALSE
, eps=1.0e-6
)
{

  auxnames <- as.list(substitute(list(...)))[-1]
  filenames <- .Primitive("c")(filenames, auxnames)

  n <- length(filenames)

  ## error if no file name !
  if (n == 0)
    stop("No file name given !")

  pdata <- pData(phenoData)
  ##try to read sample names form phenoData. if not there use CEL filenames
  if(dim(pdata)[1]!=n){#if empty pdata filename are samplenames
    warning("Incompatible phenoData object. Created a new one.\n")

    samplenames <- gsub("^/?([^/]*/)*", "", unlist(filenames))
    pdata <- data.frame(sample=1:n,row.names=samplenames)
    phenoData <- new(
      "phenoData"
     , pData=pdata
     , varLabels=list(sample="arbitrary numbering")
     )
  }
  else samplenames <- rownames(pdata)

  if (is.null(description))
    {
      description <- new("MIAME")
      description@preprocessing$filenames <- filenames
      description@preprocessing$affyversion <- library(
                                                 help=affy
                                               )$info[[2]][[2]][2]
	  description@other <- list(notes)
    }
  ## read the first file to see what we have
  ##if (verbose) cat(1, "reading",filenames[[1]],"...")

  ## get information from cdf environment

  headdetails <- read.celfile.header(filenames[[1]])
  dim.intensity <- headdetails[[2]]
  cdfname <- headdetails[[1]]

  tmp <- new("AffyBatch",
             cdfName=cdfname,
             annotation=cleancdfname(cdfname, addcdf=FALSE))
  pmIndex <- pmindex(tmp)
  probenames <- rep(names(pmIndex), unlist(lapply(pmIndex,length)))
  pmIndex <- unlist(pmIndex)

  ## read pm data into matrix

  pm <- read.probematrix(filenames=unlist(filenames),which="pm")$pm
  mm <- read.probematrix(filenames=unlist(filenames),which="mm")$mm

  ## pass matrix of probe values to mgmos
  ## call mgmos
  conds <- n
  genes <- length(featureNames(tmp))

  phis <- c(0,0,0)  
  prctiles <- 0.01*c(5, 25, 50, 75, 95);
  
  if (background == TRUE)
  {
    for (i in c(1:conds)){
      m<-min(c(min(pm[,i]),min(mm[,i])))
      pm[,i]<-pm[,i]-m+1
      mm[,i]<-mm[,i]-m+1
    }
  }

  res <-
    .Call(
      "mgmos_c"
     , pm
     , mm
     , genes
     , probenames
     , phis
     , prctiles
     , length(prctiles)
     , savepar
     , eps
     , PACKAGE="puma"
     )

  expr <- matrix(res[c(1:genes),],genes,conds)
  se <- matrix(res[c((genes+1):(2*genes)),],genes,conds)
  prc5 <- matrix(res[c((2*genes+1):(3*genes)),],genes,conds)
  prc25 <- matrix(res[c((3*genes+1):(4*genes)),],genes,conds)
  prc50 <- matrix(res[c((4*genes+1):(5*genes)),],genes,conds)
  prc75 <- matrix(res[c((5*genes+1):(6*genes)),],genes,conds)
  prc95 <- matrix(res[c((6*genes+1):(7*genes)),],genes,conds)

  rm(res)

  if (gsnorm[1]=="mean")
  {
    expr <- as.data.frame(2^expr)
    
    chipm <- apply(expr,2,mean)
    chipm <- chipm/chipm[1]

    expr <- as.matrix(log2(expr))
    for (i in 1:conds)
    {
      expr[,i] <- expr[,i]-log2(chipm[i])
      prc5[,i] <- prc5[,i]-log2(chipm[i])
      prc25[,i] <- prc25[,i]-log2(chipm[i])
      prc50[,i] <- prc50[,i]-log2(chipm[i])
      prc75[,i] <- prc75[,i]-log2(chipm[i])
      prc95[,i] <- prc95[,i]-log2(chipm[i])
    }
  }
  else if (gsnorm[1]=="median")
  {
    expr <- as.data.frame(2^expr)
    
    chipm <- apply(expr,2,median)
    chipm <- chipm/chipm[1]

    expr <- as.matrix(log2(expr))
    for (i in 1:conds)
    {
      expr[,i] <- expr[,i]-log2(chipm[i])
      prc5[,i] <- prc5[,i]-log2(chipm[i])
      prc25[,i] <- prc25[,i]-log2(chipm[i])
      prc50[,i] <- prc50[,i]-log2(chipm[i])
      prc75[,i] <- prc75[,i]-log2(chipm[i])
      prc95[,i] <- prc95[,i]-log2(chipm[i])
    }
  }
  else if (gsnorm[1]=="meanlog")
  {
    chipm <- apply(expr,2,mean)
    chipm <- chipm-chipm[1]

    for (i in 1:conds)
    {
      expr[,i] <- expr[,i]-chipm[i]
      prc5[,i] <- prc5[,i]-chipm[i]
      prc25[,i] <- prc25[,i]-chipm[i]
      prc50[,i] <- prc50[,i]-chipm[i]
      prc75[,i] <- prc75[,i]-chipm[i]
      prc95[,i] <- prc95[,i]-chipm[i]
    }
  }

  rownames(expr) <- featureNames(tmp)
  colnames(expr) <- samplenames
  rownames(se) <- featureNames(tmp)
  colnames(se) <- samplenames
  rownames(prc5) <- featureNames(tmp)
  colnames(prc5) <- samplenames
  rownames(prc25) <- featureNames(tmp)
  colnames(prc25) <- samplenames
  rownames(prc50) <- featureNames(tmp)
  colnames(prc50) <- samplenames
  rownames(prc75) <- featureNames(tmp)
  colnames(prc75) <- samplenames
  rownames(prc95) <- featureNames(tmp)
  colnames(prc95) <- samplenames

  annotation <- annotation(tmp)

  new(
    "exprReslt"
    , prcfive=prc5
    , prctwfive=prc25
    , prcfifty=prc50
    , prcsevfive=prc75
    , prcninfive=prc95
    , exprs=expr
    , se.exprs=se
    , phenoData=phenoData
    , annotation=annotation
    , experimentData=description
    )
}
