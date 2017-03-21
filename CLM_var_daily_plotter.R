library(ggplot2)
library(dplyr)

df <- read.csv('./btran_daily.csv', header=TRUE)
by_run <- group_by(select(df, case, loc, doy, value),
                   case, doy)
s <- summarize(by_run,
               count=n(),
               BTRAN=mean(value),
               minval=min(value),
               maxval=max(value))
h <- ggplot(s, aes(doy, BTRAN, group=case, col=case)) +
    geom_ribbon(aes(ymin = minval, ymax = maxval), fill="grey50", alpha=0.4) +
        geom_line(aes(y = BTRAN), linetype=2)
h2 <- ggplot(filter(s, case=='IDE_redpcp'), aes(doy, BTRAN)) +
    geom_ribbon(aes(ymin = minval, ymax = maxval), fill="grey70") +
        geom_line(aes(y = BTRAN))
## question: how to plot cases separately?




## huron <- data.frame(year = 1875:1972, level = as.vector(LakeHuron))
## h <- ggplot(huron, aes(year))
## h +
##   geom_ribbon(aes(ymin = level - 1, ymax = level + 1), fill = "grey70") +
##   geom_line(aes(y = level))
