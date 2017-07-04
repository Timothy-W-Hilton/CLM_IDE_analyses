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

sitenames_reorder_factor <- function(f) {
    ## abbreviate long site names so they fit above their plots
    ## group sites geographically rather than alphabetically
    f <- recode(
        f,
        "Sierra Foothill Research Extension Center"="Sierra Foothill",
        "Loma Ridge Global Change Experiment"="Loma Ridge",
        "ARM Southern Great Plains"="ARM S. Great Plains",
        "McLaughlin NRS"="McLaughlin",
        "Sedgewick NRS"="Sedgwick",
        "Mammoth Lakes"="SNARL",
        .default=levels(f))
    f <- factor(f, levels=c("Harvard Forest",
                            "ARM S. Great Plains",
                            "WLEF",
                            "Sierra Foothill",
                            "McLaughlin",
                            "Younger Lagoon",
                            "Box Springs",
                            "Loma Ridge",
                            "Sedgwick",
                            "SNARL",
                            "Carrizo Plain"))
    return(f)
    }

calc_annual_difference <- function(df, units_factor=1.0) {
    ## calculate annual difference between control run, drought run
    ##
    ann_sum <- summarize(group_by(df, case, loc),
                         annual_sum=sum(val * units_factor, na.rm=TRUE))
    ann_diff <- spread(ann_sum, case, annual_sum)
    ann_diff <- within(ann_diff, d <- control - drought)
    ann_diff <- within(ann_diff, pct <- d / control)
    return(ann_diff)
}

bootstrapper <- function(varname, units_factor=1.0) {
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
    return(s)
}

vert_axis_label_builder <- function(varname, units) {
    ## build axis label for vertical axis
    if (varname == 'FPSN') {
        return(parse(text='GPP~(gC~m^-2~d^-1)'))
    } else if (varname == 'BTRAN') {
        return(parse(text='beta[t]'))
    } else {
        return(paste(varname, units))
    }
}

plotter <- function(s,
                    varname,
                    plot_min=NA, plot_max=NA,
                    units="", units_factor=1.0) {
    ## plot daily CLM variable with 95% CI envelope, one panel per site
    ##
    ## ARGS:
    ##  s: summary of bootstrap returned by bootstrapper()
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

    s[['cilo']] <- unlist(lapply(s[['ci']], function(x) x[['cilo']]))
    s[['cilo']][s[['cilo']] < plot_min] <- plot_min
    s[['cihi']] <- unlist(lapply(s[['ci']], function(x) x[['cihi']]))
    s[['cihi']][s[['cihi']] > plot_max] <- plot_max
    if (is.na(plot_min)) {
        plot_min <- min(s[['cilo']], na.rm=TRUE)
    }
    if (is.na(plot_max)) {
        plot_max <- max(s[['cihi']], na.rm=TRUE)
    }
    s <- select(s, loc, case, doy, count, val, minval, maxval, cilo, cihi)
    levels(s$case) <- c('control', 'drought')
    s[['loc']] <- sitenames_reorder_factor(s[['loc']])
    ann_diff <- calc_annual_difference(s)
    ann_diff[['label']] <- paste0('Delta==', round(ann_diff[['d']]),
                                 '~', units)
    h <- ggplot(s, aes(doy, val, group=case)) +
        labs(y=vert_axis_label_builder(varname, units), x='Day of Year') +
        geom_ribbon(aes(ymin = cilo, ymax = cihi, linetype=case),
                    fill="grey50", alpha=0.4) +
        geom_line(aes(y = val, linetype=case)) +
        theme_few() +  ## https://www.r-bloggers.com/ggplot2-themes-examples/
        ylim(plot_min, plot_max) + ## BTRAN varies in [0.0, 1.0]
        facet_wrap(~ loc, ncol=3 ) +
        theme(legend.position=c(0.8, 0.1)) ## http://www.cookbook-r.com/Graphs/Legends <- (ggplot2)/
    fname <- paste(varname, '_daily_sites.pdf', sep='')
    cat(paste('saving', fname, '...'))
    ggsave(filename=fname)
    cat('done\n')
    return(list(annsum=s, anndiff=ann_diff))
}

do_bootstrap <- FALSE
if (do_bootstrap) {
    s_btran <- bootstrapper('BTRAN')
    s_FPSN <- bootstrapper('FPSN')
    s_rain <- bootstrapper('RAIN')
    s_wt <- bootstrapper('WT')
    s_h2osoi_lev1 <- bootstrapper('H2OSOIlev00')
    h2osoi_sum <- bootstrapper('H2OSOIsum')
}

s_per_day <- 24*60*60
mw_C <- 12.0107
umol_per_mol <- 1e-6
umol_m2_s_2_gC_m2_d <- s_per_day * mw_C * umol_per_mol

btran <- plotter(s_btran, 'BTRAN', plot_min=0.0, plot_max=1.0)
fpsn <- plotter(s_FPSN,
                'FPSN',
                units_factor=umol_m2_s_to_gC_m2_d,
                units='(gC/m2/d)')

rain <- plotter(s_rain, 'RAIN', plot_min=0.0, units='(mm/d)',
                units_factor=s_per_day)
wt <- plotter(s_wt, 'WT', plot_min=0.0, units='(mm)')
h2osoi_lev1 <- plotter(s_h2osoi_lev1, 'H2OSOIlev0', units='(mm3/mm3)')
h2osoi_sum <- plotter(h2osoi_sum, 'H2OSOIsum', units='(mm3/mm3)')
