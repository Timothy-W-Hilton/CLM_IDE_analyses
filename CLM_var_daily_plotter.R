library(ggplot2)
library(ggthemes)
library(dplyr)
library(boot)

DEBUGFLAG <- TRUE

bootThetaMean <- function(x,i) {
    mean(x[i])
}

boot_5_95 <- function(vals, R=1000) {
    if (length(unique(vals)) > 1) {
        b <- boot(vals, R=1000, statistic=bootThetaMean)
        bci <- boot.ci(b, type='basic')
        ## see boot.ci documentation: basic[4:5] contain the
        ## confidence interval bounds (cilo = "CI low", cihi = "CI
        ## high")
        cilo <- bci[['basic']][4]
        cihi <- bci[['basic']][5]
        ci <- c(cilo=cilo, cihi=cihi)
    } else {
        unique_val <- unique(vals)
        ci <- c(cilo=unique_val, cihi=unique_val)
    }
    return(ci)
}

plotter <- function(varname,
                    plot_min=NA, plot_max=NA,
                    units="", units_factor=1.0) {
    ## plot daily CLM variable with 95% CI envelope, one panel per site
    ##
    ## ARGS:
    ##  varname (string): name of variable to be plotted
    ##  plot_min (float): vertical axis minimum (default is min(95%
    ##      confidence interval)
    ##  plot_max (float): vertical axis maximum (default is max(95%
    ##      confidence interval)
    ##  units (string): units, to be appended to vertical axis label
    ##  units_factor (float): units conversion factor to be multiplied
    ##      into requested variable.  Default is 1.0.
    ##
    ## RETURNS:
    ##  data frame containing daily values, labeled by site

    ## TODO: min, max vals should be arguments
    ## TODO: where to put scale factor?  Argument here, or in shell
    ## script that builds daily netcdf files?
    fname_csv_base <- paste(varname, '_daily_all.csv.gz', sep='')
    df <- read.csv(file.path('/', 'global', 'cscratch1', 'sd',
                             'twhilton', 'daily_CLM_output', 'output',
                             fname_csv_base),
                   header=TRUE)
    df[['value']] <- df[['value']] * units_factor
    by_run <- group_by(select(df, case, loc, doy, value),
                       case, doy, loc)
    if (DEBUGFLAG){
        ibootstrap <- 5  ## bootstrap iterations
    } else {
        ibootstrap <- 1000
    }
    s <- summarize(by_run,
                   count=n(),
                   val=mean(value),
                   minval=min(value),
                   maxval=max(value),
                   ci=list(boot_5_95(value, R=ibootstrap)))
    s[['cilo']] <- unlist(lapply(s[['ci']], function(x) x[['cilo']]))
    s[['cilo']][s[['cilo']] < plot_min] <- plot_min
    s[['cihi']] <- unlist(lapply(s[['ci']], function(x) x[['cihi']]))
    s[['cihi']][s[['cihi']] > plot_max] <- plot_max
    if (is.na(plot_min)) {
        plot_min <- min(s[['cilo']], na.rm=TRUE)
        cat(paste('plot_min', plot_min))
    }
    if (is.na(plot_max)) {
        plot_max <- max(s[['cihi']], na.rm=TRUE)
        cat(paste('plot_min', plot_min))
    }
    s <- select(s, loc, case, doy, count, val, minval, maxval, cilo, cihi)
    levels(s$case) <- c('control', 'drought')
    h <- ggplot(s, aes(doy, val, group=case)) +
        labs(y=paste(varname, units)) +
        geom_ribbon(aes(ymin = cilo, ymax = cihi, linetype=case),
                    fill="grey50", alpha=0.4) +
        geom_line(aes(y = val, linetype=case)) +
        theme_few() +  ## https://www.r-bloggers.com/ggplot2-themes-examples/
        ylim(plot_min, plot_max) + ## BTRAN varies in [0.0, 1.0]
        facet_wrap(~ loc, ncol=3 )
    fname <- paste(varname, '_daily_sites.pdf', sep='')
    cat(paste('saving', fname, '...'))
    ggsave(filename=fname)
    cat('done\n')
    return(s)
}

s_per_day <- 24*60*60
rain <- plotter('RAIN', plot_min=0.0, units='(mm/d)', units_factor=s_per_day)
btran <- plotter('BTRAN', plot_min=0.0, plot_max=1.0)
wt <- plotter('WT', plot_min=0.0, units='(mm)')
fpsn <- plotter('FPSN', units='(umol/m2/s)')
