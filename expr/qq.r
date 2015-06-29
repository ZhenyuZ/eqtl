# Old qqman.r qq plot code from github stephenturner/qqman
# This version support confidence interview and a few more functions compared to later versions

qq = function(pvector,gridlines=F,gridlines.col='gray83',gridlines.lwd=1,gridlines.lty=1,confidence=T,confidence.col='gray81',
    pt.cex=0.5,pt.col='black',pt.bg='black',pch=21,abline.col='red',abline.lwd=1.8,abline.lty=1,ymax=8,ymax.soft=T,
    highlight=NULL,highlight.col=c('green3','magenta'),highlight.bg=c('green3','magenta'),
    annotate=NULL,annotate.cex=0.7,annotate.font=3,cex.axis=0.95,...) {
    #======================================================================================================
    ######## Check data and arguments; create observed and expected distributions
    d = suppressWarnings(as.numeric(pvector))
    names(d) = names(pvector)
    d = d[!is.na(d)] # remove NA, and non-numeric [which were converted to NA during as.numeric()]
    d = d[d>0 & d<1] # only Ps between 0 and 1
    
    
    if (!is.null(highlight) | !is.null(annotate)){
        if (is.null(names(d))) stop("P-value vector must have names to use highlight or annotate features.")
        d = d[!is.na(names(d))]
        if (!is.null(highlight) & FALSE %in% (highlight %in% names(d))) stop ("D'oh! Highlight vector must be a subset of names(pvector).")
        if (!is.null(annotate) & FALSE %in% (annotate %in% names(d))) stop ("D'oh! Annotate vector must be a subset of names(pvector).")
    }
    
    d = d[order(d,decreasing=F)] # sort
    o = -log10(d)
    e = -log10( ppoints(length(d) ))
    if (!is.null(highlight) | !is.null(annotate)) names(e) = names(o) = names(d)
    
    if (!is.numeric(ymax) | ymax<max(o)) ymax <- max(o) 
    
    if (!is.numeric(pt.cex) | pt.cex<0) pt.cex=0.5
    if (!is.numeric(annotate.cex) | annotate.cex<0) annotate.cex=0.7
    if (!is.numeric(annotate.font)) annotate.font=3
    
    if (is.character(gridlines.col[1]) & !(gridlines.col[1] %in% colors())) gridlines.col = 'gray83'
    if (is.character(confidence.col[1]) & !(confidence.col[1] %in% colors())) confidence.col = 'gray81'
    if (is.character(abline.col[1]) & !(abline.col[1] %in% colors())) abline.col = 'red'
    
    if (FALSE %in% (pt.col %in% colors() | !is.na(suppressWarnings(as.numeric(pt.col))) )){
        pt.col = 'black'; warning("pt.col argument(s) not recognized. Setting to default: 'black'.")
    }

    if (FALSE %in% (pt.bg %in% colors() | !is.na(suppressWarnings(as.numeric(pt.bg))) )){
        pt.bg = 'black'; warning("pt.bg argument(s) not recognized. Setting to default: 'black'.")
    }
    
    if (FALSE %in% (highlight.col %in% colors() | !is.na(suppressWarnings(as.numeric(highlight.col))) )){
        highlight.col = 'blue'; warning("highlight.col argument(s) not recognized. Setting to default: 'blue'.")
    }

    if (FALSE %in% (highlight.bg %in% colors() | !is.na(suppressWarnings(as.numeric(highlight.bg))) )){
        highlight.bg = 'blue'; warning("highlight.bg argument(s) not recognized. Setting to default: 'blue'.")
    }
    
    # Ymax
    if(is.na(suppressWarnings(as.numeric(ymax)))){  # not numeric
        ymax = ceiling(max(o))
        warning('non-numeric ymax argument.')
    } else if (as.numeric(ymax) < 0){           # negative
        ymax = ceiling(max(o))
        warning('negative ymax argument.')
    }
    if (ymax.soft==T){ #if soft, ymax is just the lower limit for ymax
        ymax = max(ymax, ceiling(max(o)))
    } #else, ymax = ymax
            
    
    ################################
    
    # Initialize plot
    #print('Setting up plot.')
    #print(ymax)
    xspace = 0.078
    xmax = max(e) * 1.019
    xmin = max(e) * -0.035
    #ymax = ceiling(ymax * 1.03)
    ymin = -ymax*0.03
    plot(0,xlab=expression(Expected~~-log[10](italic(p))),ylab=expression(Observed~~-log[10](italic(p))),
            col=F,las=1,xaxt='n',xlim=c(xmin,xmax),ylim=c(ymin,ymax),bty='n',xaxs='i',yaxs='i',cex.axis=cex.axis, ...)
    axis(side=1,labels=seq(0,max(e),1),at=seq(0,max(e),1),cex.axis=cex.axis,lwd=0,lwd.ticks=1)
    
    # Grid lines
    if (isTRUE(gridlines)){
        yvals = par('yaxp')
        yticks = seq(yvals[1],yvals[2],yvals[2]/yvals[3])
        abline(v=seq(0,max(e),1),col=gridlines.col[1],lwd=gridlines.lwd,lty=gridlines.lty)
        abline(h=yticks,col=gridlines.col[1],lwd=gridlines.lwd,lty=gridlines.lty)
    }
    
     #Confidence intervals
     find_conf_intervals = function(row){
        i = row[1]
        len = row[2]
        if (i < 10000 | i %% 100 == 0){
            return(c(-log10(qbeta(0.95,i,len-i+1)), -log10(qbeta(0.05,i,len-i+1))))
        } else { # Speed up
            return(c(NA,NA))
        }
     }

     # Find approximate confidence intervals
    if (isTRUE(confidence)){
        #print('Plotting confidence intervals.')
        ci = apply(cbind( 1:length(e), rep(length(e),length(e))), MARGIN=1, FUN=find_conf_intervals)
        bks = append(seq(10000,length(e),100),length(e)+1)
        for (i in 1:(length(bks)-1)){
            ci[1, bks[i]:(bks[i+1]-1)] = ci[1, bks[i]]
            ci[2, bks[i]:(bks[i+1]-1)] = ci[2, bks[i]]
        }
        colnames(ci) = names(e)
        # Extrapolate to make plotting prettier (doesn't affect intepretation at data points)
        slopes = c((ci[1,1] - ci[1,2]) / (e[1] - e[2]), (ci[2,1] - ci[2,2]) / (e[1] - e[2]))
        extrap_x = append(e[1]+xspace,e) #extrapolate slightly for plotting purposes only
        extrap_y = cbind( c(ci[1,1] + slopes[1]*xspace, ci[2,1] + slopes[2]*xspace), ci)
        
        polygon(c(extrap_x, rev(extrap_x)), c(extrap_y[1,], rev(extrap_y[2,])),col = confidence.col[1], border = confidence.col[1]) 
    }
    
    # Points (with optional highlighting)
    #print('Plotting data points.')
    fills = rep(pt.bg,length(o))
    borders = rep(pt.col,length(o))
    names(fills) = names(borders) = names(o)
    if (!is.null(highlight)){   
        borders[highlight] = rep(NA,length(highlight))
        fills[highlight] = rep(NA,length(highlight))
    }
    points(e,o,pch=pch,cex=pt.cex,col=borders,bg=fills)
    
    if (!is.null(highlight)){
        points(e[highlight],o[highlight],pch=pch,cex=pt.cex,col=highlight.col,bg=highlight.bg)
    }
    
    #Abline
    abline(0,1,col=abline.col,lwd=abline.lwd,lty=abline.lty)
    
    # Annotate SNPs
    if (!is.null(annotate)){
        x = e[annotate] # x will definitely be the same
        y = -0.1 + apply(rbind(o[annotate],ci[1,annotate]),2,min)
        text(x,y,labels=annotate,srt=90,cex=annotate.cex,adj=c(1,0.48),font=annotate.font)      
    }
    # Box
    box()
}







