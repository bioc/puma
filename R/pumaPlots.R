plotHistTwoClasses <- function (
	scores			# A set of scores (e.g. output of pplr function)
,	class1Elements	# Logical vector, TRUE if element is in first class
,	class2Elements	# Logical vector, TRUE if element is in second class
,	space=0			# no space between bars
,	col=c("white", "grey40")
,	xlab="PPLR"
,	ylab="Number of genes"
# ,	ylim=c(0, 5000)
,	ylim=NULL
,	las=0 # axis labels all perpendicular to axes
,	legend=c("non-spike-in genes", "spike-in genes")
,	inset=0.05
,	minScore=0
,	maxScore=1
,	numOfBars=20
,	main=NULL
)
{
	class1Bars <- hist(
						scores[class1Elements]
					,	breaks=seq(
							minScore
						,	maxScore
						,	(maxScore-minScore)/numOfBars
						)
					,	plot=FALSE
					)
	class2Bars <- hist(
						scores[class2Elements]
					,	breaks=seq(
							minScore
						,	maxScore
						,	(maxScore-minScore)/numOfBars
						)
					,	plot=FALSE
					)
	barsToPlot <- rbind(class2Bars$counts, class1Bars$counts)
	barplot(
		barsToPlot
	,	space=space
	,	col=col
	,	xlab=xlab
	,	ylab=ylab
	,	ylim=ylim
	,	main=main
	)
	axis(1, at=0:numOfBars, labels=class1Bars$breaks)
	box()
	legend(
		"topright"
	,	legend=legend
	,	fill=col
	,	inset=inset
	)
}

plotWhiskers <- function (
	eset
,	comparisons=c(1,2)
,	sortMethod = c("logRatio", "PPLR")
,	numGenes=50
,	xlim
,	main = "PUMA Whiskers plot"
,	highlightedGenes=NULL
)
{
	if(length(comparisons) != 2)
		stop("Must supply a 2-vector for comparisons (columns to compare)")
	if(max(comparisons) > dim(exprs(eset))[2])
		stop("Maximum comparison must be less than number of samples")
	logRatio <- abs(exprs(eset[,comparisons[1]])-exprs(eset[,comparisons[2]]))
	if(sortMethod[1] == "logRatio")
	{
		sortLogRatio <- sort(logRatio, decreasing=TRUE, method="quick", index.return=TRUE)[[2]]
		topGenes <- sortLogRatio[1:numGenes]
	}
	else
	{
		topGenes <- pumaDE(eset[,comparisons])[[2]][1:numGenes]
	}
	topLogRatio <- logRatio[topGenes]
	topPlotPoints <- cbind(topLogRatio, -1:-numGenes)
	SEs <- sqrt(assayDataElement(eset[,comparisons[1]],"se.exprs")^2 + assayDataElement(eset[,comparisons[2]],"se.exprs")^2)
	topSE <- SEs[topGenes]
	plot(topPlotPoints, xlim=xlim, pch=8, main=main)
	segments(topLogRatio + qnorm(0.05)*topSE, -1:-numGenes, topLogRatio + qnorm(0.95)*topSE, -1:-numGenes)
	topGenesSpikeIn <- which(!is.na(match(topGenes,highlightedGenes)))
	points(topPlotPoints[topGenesSpikeIn,1], topPlotPoints[topGenesSpikeIn,2],pch=22,cex=2)
	abline(v=0)
}

plotROC <- function (
	scoresList
,	truthValues
,	includedProbesets=1:length(truthValues)
,	legendTitles=1:length(scoresList)
,	main = "PUMA ROC plot"
,	lty = 1:length(scoresList)
,	col = rep(1,length(scoresList))
,	lwd = rep(1,length(scoresList))
,	yaxisStat = "tpr"
,	xaxisStat = "fpr"
,	downsampling = 100
,	showLegend = TRUE
,	showAUC = TRUE
,	...
)
{
	require(ROCR)
	
	predictions <- list()
	ROCs <- list()
	AUCs <- list()
	legend <- NULL
	for(i in 1:length(scoresList))
	{
		predictions[[i]] <- prediction(scoresList[[i]][includedProbesets]
			,	truthValues[includedProbesets])
		ROCs[[i]] <- performance(predictions[[i]], yaxisStat, xaxisStat)
		AUCs[[i]] <- performance(predictions[[i]], "auc")@y.values[[1]]
		if(i==1)
			plot(ROCs[[i]], lty=lty[i], main=main, col=col[i], lwd=lwd[i], downsampling=downsampling, ...)
		else
			plot(ROCs[[i]], lty=lty[i], add=TRUE, col=col[i], lwd=lwd[i], downsampling=downsampling, ...)
		if(showAUC)
			legend = c(legend, paste(legendTitles[i],"AUC = ", round(AUCs[[i]],3)))
		if(!showAUC)
			legend = c(legend, paste(legendTitles[i]))
	}
	if(showLegend)
	{
		legend2(
			"bottomright"
		,	legend=legend
		,	lty=lty
		,	col=col
		,	inset=0.05
		,	seg.len=10
		# ,	...
		)
	}
}

calcAUC <- function (
	scores
,	truthValues
,	includedProbesets=1:length(truthValues)
)
{
	require(ROCR)
	
	predictions <- prediction(scores[includedProbesets]
		,	truthValues[includedProbesets])
	AUC <- performance(predictions, "auc")@y.values[[1]]
}

plot.pumaPCARes <- function(
	x
,	...
,	firstComponent=1
,	secondComponent=2
,	useFilenames=FALSE
,	phenotype = pData(x@phenoData)
,	legend1pos = "topright"
,	legend2pos = "bottomright"
)
{
	plot(
		x=x@model@W[,firstComponent]
	,	y=x@model@W[,secondComponent]
	,	type=if(useFilenames) "n" else "p"
	,	xlab=paste("Component", firstComponent)
	,	ylab=paste("Component", secondComponent)
	,	pch=if(dim(phenotype)[2] >= 1)
				unclass(as.factor(phenotype[,1]))
			else 1
	,	col=if(dim(phenotype)[2] >= 2)
				unclass(as.factor(phenotype[,2]))
			else if(dim(phenotype)[2] == 1)
				unclass(as.factor(phenotype[,1]))
			else
				1
	,	...
	)
	if(useFilenames)
		text(
			x=x@model@W[,1]
		,	y=x@model@W[,2]
		,	labels=rownames(phenotype)
		)
    if(dim(phenotype)[2] >= 1)
		legend(
	      x=legend1pos
	    , legend=levels(as.factor(phenotype[,1]))
	    , pch=1:nlevels(as.factor(phenotype[,1]))
		, col=1:nlevels(as.factor(phenotype[,1]))
		, title=colnames(phenotype)[1]
	    )
	if(dim(phenotype)[2] >= 2)
	    legend(
	      x=legend2pos
	    , legend=levels(as.factor(phenotype[,2]))
	    , text.col=1:nlevels(as.factor(phenotype[,2]))
		, title=colnames(phenotype)[2]
	    )
}

legend2 <- function (x, y = NULL, legend, fill = NULL, col = par("col"), 
    lty, lwd, pch, angle = 45, density = NULL, bty = "o", bg = par("bg"), 
    box.lwd = par("lwd"), box.lty = par("lty"), pt.bg = NA, cex = 1, 
    pt.cex = cex, pt.lwd = lwd, xjust = 0, yjust = 1, x.intersp = 1, 
    y.intersp = 1, adj = c(0, 0.5), text.width = NULL, text.col = par("col"), 
    merge = do.lines && has.pch, trace = FALSE, plot = TRUE, 
    ncol = 1, horiz = FALSE, title = NULL, inset = 0, seg.len = 2) 
{
    if (missing(legend) && !missing(y) && (is.character(y) || 
        is.expression(y))) {
        legend <- y
        y <- NULL
    }
    mfill <- !missing(fill) || !missing(density)
    if (length(title) > 1) 
        stop("invalid title")
    n.leg <- if (is.call(legend)) 
        1
    else length(legend)
    if (n.leg == 0) 
        stop("'legend' is of length 0")
    auto <- if (is.character(x)) 
        match.arg(x, c("bottomright", "bottom", "bottomleft", 
            "left", "topleft", "top", "topright", "right", "center"))
    else NA
    if (is.na(auto)) {
        xy <- xy.coords(x, y)
        x <- xy$x
        y <- xy$y
        nx <- length(x)
        if (nx < 1 || nx > 2) 
            stop("invalid coordinate lengths")
    }
    else nx <- 0
    xlog <- par("xlog")
    ylog <- par("ylog")
    rect2 <- function(left, top, dx, dy, density = NULL, angle, 
        ...) {
        r <- left + dx
        if (xlog) {
            left <- 10^left
            r <- 10^r
        }
        b <- top - dy
        if (ylog) {
            top <- 10^top
            b <- 10^b
        }
        rect(left, top, r, b, angle = angle, density = density, 
            ...)
    }
    segments2 <- function(x1, y1, dx, dy, ...) {
        x2 <- x1 + dx
        if (xlog) {
            x1 <- 10^x1
            x2 <- 10^x2
        }
        y2 <- y1 + dy
        if (ylog) {
            y1 <- 10^y1
            y2 <- 10^y2
        }
        segments(x1, y1, x2, y2, ...)
    }
    points2 <- function(x, y, ...) {
        if (xlog) 
            x <- 10^x
        if (ylog) 
            y <- 10^y
        points(x, y, ...)
    }
    text2 <- function(x, y, ...) {
        if (xlog) 
            x <- 10^x
        if (ylog) 
            y <- 10^y
        text(x, y, ...)
    }
    if (trace) 
        catn <- function(...) do.call("cat", c(lapply(list(...), 
            formatC), list("\n")))
    cin <- par("cin")
    Cex <- cex * par("cex")
    if (is.null(text.width)) 
        text.width <- max(abs(strwidth(legend, units = "user", 
            cex = cex)))
    else if (!is.numeric(text.width) || text.width < 0) 
        stop("'text.width' must be numeric, >= 0")
    xc <- Cex * xinch(cin[1], warn.log = FALSE)
    yc <- Cex * yinch(cin[2], warn.log = FALSE)
    if (xc < 0) 
        text.width <- -text.width
    xchar <- xc
    xextra <- 0
    yextra <- yc * (y.intersp - 1)
    ymax <- yc * max(1, strheight(legend, units = "user", cex = cex)/yc)
    ychar <- yextra + ymax
    if (trace) 
        catn("  xchar=", xchar, "; (yextra,ychar)=", c(yextra, 
            ychar))
    if (mfill) {
        xbox <- xc * 0.8
        ybox <- yc * 0.5
        dx.fill <- xbox
    }
    do.lines <- (!missing(lty) && (is.character(lty) || any(lty > 
        0))) || !missing(lwd)
    n.legpercol <- if (horiz) {
        if (ncol != 1) 
            warning("horizontal specification overrides: Number of columns := ", 
                n.leg)
        ncol <- n.leg
        1
    }
    else ceiling(n.leg/ncol)
    if (has.pch <- !missing(pch) && length(pch) > 0) {
        if (is.character(pch) && !is.na(pch[1]) && nchar(pch[1], 
            type = "c") > 1) {
            if (length(pch) > 1) 
                warning("not using pch[2..] since pch[1] has multiple chars")
            np <- nchar(pch[1], type = "c")
            pch <- substr(rep.int(pch[1], np), 1:np, 1:np)
        }
        if (!merge) 
            dx.pch <- x.intersp/2 * xchar
    }
    x.off <- if (merge) 
        -0.7
    else 0
    if (is.na(auto)) {
        if (xlog) 
            x <- log10(x)
        if (ylog) 
            y <- log10(y)
    }
    if (nx == 2) {
        x <- sort(x)
        y <- sort(y)
        left <- x[1]
        top <- y[2]
        w <- diff(x)
        h <- diff(y)
        w0 <- w/ncol
        x <- mean(x)
        y <- mean(y)
        if (missing(xjust)) 
            xjust <- 0.5
        if (missing(yjust)) 
            yjust <- 0.5
    }
    else {
        h <- (n.legpercol + (!is.null(title))) * ychar + yc
        w0 <- text.width + (x.intersp + 1) * xchar
        if (mfill) 
            w0 <- w0 + dx.fill
        if (has.pch && !merge) 
            w0 <- w0 + dx.pch
        if (do.lines) 
            w0 <- w0 + (seg.len + x.off) * xchar
        w <- ncol * w0 + 0.5 * xchar
        if (!is.null(title) && (tw <- strwidth(title, units = "user", 
            cex = cex) + 0.5 * xchar) > w) {
            xextra <- (tw - w)/2
            w <- tw
        }
        if (is.na(auto)) {
            left <- x - xjust * w
            top <- y + (1 - yjust) * h
        }
        else {
            usr <- par("usr")
            inset <- rep(inset, length.out = 2)
            insetx <- inset[1] * (usr[2] - usr[1])
            left <- switch(auto, bottomright = , topright = , 
                right = usr[2] - w - insetx, bottomleft = , left = , 
                topleft = usr[1] + insetx, bottom = , top = , 
                center = (usr[1] + usr[2] - w)/2)
            insety <- inset[2] * (usr[4] - usr[3])
            top <- switch(auto, bottomright = , bottom = , bottomleft = usr[3] + 
                h + insety, topleft = , top = , topright = usr[4] - 
                insety, left = , right = , center = (usr[3] + 
                usr[4] + h)/2)
        }
    }
    if (plot && bty != "n") {
        if (trace) 
            catn("  rect2(", left, ",", top, ", w=", w, ", h=", 
                h, ", ...)", sep = "")
        rect2(left, top, dx = w, dy = h, col = bg, density = NULL, 
            lwd = box.lwd, lty = box.lty)
    }
    xt <- left + xchar + xextra + (w0 * rep.int(0:(ncol - 1), 
        rep.int(n.legpercol, ncol)))[1:n.leg]
    yt <- top - 0.5 * yextra - ymax - (rep.int(1:n.legpercol, 
        ncol)[1:n.leg] - 1 + (!is.null(title))) * ychar
    if (mfill) {
        if (plot) {
            fill <- rep(fill, length.out = n.leg)
            rect2(left = xt, top = yt + ybox/2, dx = xbox, dy = ybox, 
                col = fill, density = density, angle = angle, 
                border = "black")
        }
        xt <- xt + dx.fill
    }
    if (plot && (has.pch || do.lines)) 
        col <- rep(col, length.out = n.leg)
    if (missing(lwd)) 
        lwd <- par("lwd")
    if (do.lines) {
        if (missing(lty)) 
            lty <- 1
        lty <- rep(lty, length.out = n.leg)
        lwd <- rep(lwd, length.out = n.leg)
        ok.l <- !is.na(lty) & (is.character(lty) | lty > 0)
        if (trace) 
            catn("  segments2(", xt[ok.l] + x.off * xchar, ",", 
                yt[ok.l], ", dx=", seg.len * xchar, ", dy=0, ...)")
        if (plot) 
            segments2(xt[ok.l] + x.off * xchar, yt[ok.l], dx = seg.len * 
                xchar, dy = 0, lty = lty[ok.l], lwd = lwd[ok.l], 
                col = col[ok.l])
        xt <- xt + (seg.len + x.off) * xchar
    }
    if (has.pch) {
        pch <- rep(pch, length.out = n.leg)
        pt.bg <- rep(pt.bg, length.out = n.leg)
        pt.cex <- rep(pt.cex, length.out = n.leg)
        pt.lwd <- rep(pt.lwd, length.out = n.leg)
        ok <- !is.na(pch) & (is.character(pch) | pch >= 0)
        x1 <- (if (merge) 
            xt - (seg.len/2) * xchar
        else xt)[ok]
        y1 <- yt[ok]
        if (trace) 
            catn("  points2(", x1, ",", y1, ", pch=", pch[ok], 
                ", ...)")
        if (plot) 
            points2(x1, y1, pch = pch[ok], col = col[ok], cex = pt.cex[ok], 
                bg = pt.bg[ok], lwd = pt.lwd[ok])
        if (!merge) 
            xt <- xt + dx.pch
    }
    xt <- xt + x.intersp * xchar
    if (plot) {
        if (!is.null(title)) 
            text2(left + w/2, top - ymax, labels = title, adj = c(0.5, 
                0), cex = cex, col = text.col)
        text2(xt, yt, labels = legend, adj = adj, cex = cex, 
            col = text.col)
    }
    invisible(list(rect = list(w = w, h = h, left = left, top = top), 
        text = list(x = xt, y = yt)))
}
