library(ggplot2)
library(dplyr)
library(boot)

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

df <- read.csv('./btran_daily.csv', header=TRUE)
by_run <- group_by(select(df, case, loc, doy, value),
                   case, doy)
s <- summarize(by_run,
               count=n(),
               BTRAN=mean(value),
               minval=min(value),
               maxval=max(value),
               ci=list(boot_5_95(value, R=5)))  ## needs a single value
s[['cilo']] <- unlist(lapply(s[['ci']], function(x) x[['cilo']]))
s[['cihi']] <- unlist(lapply(s[['ci']], function(x) x[['cihi']]))
s <- select(s, case, doy, count, BTRAN, minval, maxval, cilo, cihi)


group_colors <- c(control='#8c510a', drought='#01665e')
levels(s$case) <- c('control', 'drought')
h <- ggplot(s, aes(doy, BTRAN, group=case, col=case)) +
    geom_ribbon(aes(ymin = cilo, ymax = cihi),
                fill="grey50", alpha=0.4) +
    geom_line(aes(y = BTRAN), linetype=2) +
    scale_colour_manual(values=group_colors)
