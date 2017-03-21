library(ggplot2)
library(dplyr)
library(boot)

df <- read.csv('./btran_daily.csv', header=TRUE)
by_run <- group_by(select(df, case, loc, doy, value),
                   case, doy)
s <- summarize(by_run,
               count=n(),
               BTRAN=mean(value),
               minval=min(value),
               maxval=max(value),
               pctl05=(quantile(value, 0.05)),
               pctl95=(quantile(value, 0.95)))
group_colors <- c(control='#8c510a', drought='#01665e')
levels(s$case) <- c('control', 'drought')
h <- ggplot(s, aes(doy, BTRAN, group=case, col=case)) +
    geom_ribbon(aes(ymin = pctl05, ymax = pctl95), fill="grey50", alpha=0.4) +
    geom_line(aes(y = BTRAN), linetype=2) +
    scale_colour_manual(values=group_colors)
h2 <- ggplot(filter(s, case=='IDE_redpcp'), aes(doy, BTRAN)) +
    geom_ribbon(aes(ymin = minval, ymax = maxval), fill="grey70") +
        geom_line(aes(y = BTRAN))
## question: how to plot cases separately?




## huron <- data.frame(year = 1875:1972, level = as.vector(LakeHuron))
## h <- ggplot(huron, aes(year))
## h +
##   geom_ribbon(aes(ymin = level - 1, ymax = level + 1), fill = "grey70") +
##   geom_line(aes(y = level))
