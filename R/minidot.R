#!/usr/bin/env Rscript
## Read minimap .paf files and generate a-v-a dotplots

## * minimap OUTPUT FORMAT
##
## Minimap outputs mapping positions in the Pairwise mApping Format (PAF). PAF
## is a TAB-delimited text format with each line consisting of at least 12
## fields as are described in the following table:
##
## ┌────┬────────┬─────────────────────────────────────────────────────────────┐
## │Col │  Type  │                         Description                         │
## ├────┼────────┼─────────────────────────────────────────────────────────────┤
## │  1 │ string │ Query sequence name                                         │
## │  2 │  int   │ Query sequence length                                       │
## │  3 │  int   │ Query start coordinate (0-based)                            │
## │  4 │  int   │ Query end coordinate (0-based)                              │
## │  5 │  char  │ `+' if query and target on the same strand; `-' if opposite │
## │  6 │ string │ Target sequence name                                        │
## │  7 │  int   │ Target sequence length                                      │
## │  8 │  int   │ Target start coordinate on the original strand              │
## │  9 │  int   │ Target end coordinate on the original strand                │
## │ 10 │  int   │ Number of matching bases in the mapping                     │
## │ 11 │  int   │ Number bases, including gaps, in the mapping                │
## │ 12 │  int   │ Mapping quality (0-255 with 255 for missing)                │
## └────┴────────┴─────────────────────────────────────────────────────────────┘

## * read args/input

library(ggplot2)
library(scales)
library(stringr)

args = commandArgs(trailingOnly=TRUE)
if (length(args)<1) {
  stop("Usage: ./minidot.R PAF LEN [OUT]", call.=FALSE)
}

## debug
#setwd("/home/thackl/projects/coding/sandbox/R-minimap-dotplots")
#args <- c("minidot.paf", "minidot.len");

## * read paf and len
paf <- read.table(args[1])
len <- read.table(args[2]) #, stringAsFactor=FALSE)
pdf.file <- if (length(args)==3) args[3] else "minidot.pdf"

len$V1 <- factor(len$V1, levels=len$V1)

paf$V1 <- factor(paf$V1, levels=len$V1)
paf$V6 <- factor(paf$V6, levels=len$V1)

paf$ava <- paf$V1:paf$V6
paf$strand <- ifelse(paf$V5=='+', 1, -1)
paf[paf$strand==-1,8:9] <- paf[paf$strand==-1,9:8]
paf$idy <- paf$V10 / paf$V11 * paf$strand

## * map contig boundaries to gglayer

len.cum <- cbind(len, cumsum=c(lapply(split(len, len$V1),
               function(x) cumsum(x$V2)), recursive=T))
yava <- data.frame(ava=character(0), yi=numeric(0))
xava <- data.frame(ava=character(0), xi=numeric(0))
yava.rt <- data.frame(ava=character(0), xmin=numeric(0), xmax=numeric(0), ymin=numeric(0), ymax=numeric(0))
xava.rt <- data.frame(ava=character(0), xmin=numeric(0), xmax=numeric(0), ymin=numeric(0), ymax=numeric(0))

for(i in unique(paf$ava)){
    r <- str_split_fixed(i, ":", 2)
    yi <- c(0,len.cum$cumsum[len.cum$V1==r[2]])
    xi <- c(0,len.cum$cumsum[len.cum$V1==r[1]])

    yi.even <- if (length(yi) %% 2) yi[-length(yi)] else yi
    xi.even <- if (length(xi) %% 2) xi[-length(xi)] else xi
    
    ## odd or even
    yi.mt <- matrix(yi.even, ncol=2, byrow=T)
    xi.mt <- matrix(xi.even, ncol=2, byrow=T)

   
    yava <- rbind(yava, data.frame(ava=i, yi=yi))
    yava.rt <- rbind(yava.rt, data.frame(ava=i,
                         xmin=-Inf, xmax=Inf,
                         ymin=yi.mt[,1], ymax=yi.mt[,2]))
    xava <- rbind(xava, data.frame(ava=i, xi=xi))
    xava.rt <- rbind(xava.rt, data.frame(ava=i,
                         xmin=xi.mt[,1], xmax=xi.mt[,2],
                         ymin=-Inf, ymax=Inf))
}

gg <- ggplot(paf)
gg <- gg + geom_rect(data=xava.rt, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey75", alpha=0.5)
gg <- gg + geom_rect(data=yava.rt, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey75", alpha=0.5)
gg <- gg + geom_segment(aes(x=V3, xend=V4, y=V8, yend=V9, color=idy), size=1.5)
#gg <- gg + geom_vline(xava, aes(xintercept=xi), size=.1)
#gg <- gg + geom_hline(yava, aes(yintercept=yi), size=.1)
gg <- gg + scale_colour_distiller(palette="Spectral", direction=1)
gg <- gg + coord_fixed(1) + theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "grey55")
)

gg <- gg + scale_x_continuous(label=scientific_format(digits=0), expand=c(0,0))
gg <- gg + scale_y_continuous(label=scientific_format(digits=0), expand=c(0,0))

samples.n <- length(unique(paf$V1))
#gg <- gg + facet_wrap(~ava, ncol=samples.n, drop=FALSE)
gg <- gg + facet_grid(V6~V1, drop=FALSE, as.table=FALSE)

ggsave(pdf.file, plot=gg, width=10, height=10)
ggsave(paste(pdf.file, ".png", sep=""), plot=gg, width=10, height=10)
