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

library(ggplot2)
library(scales)
library(stringr)
library(argparse)
library(proto)

## * aux
humanize_format <- function(...) {
    function(x) humanize(x, ...)
}

humanize <- function(x, digits=2, sep=" ", unit="bp"){
  pre <- c("","k","M","G","T");
  px <- as.integer(ifelse(x > 0, log10(x)/3, 0))
  v <- signif(x/(10^(px*3)), digits=digits)
  pre[!pre==''] <- paste(sep, pre[!pre==''], sep='')
  return(paste(v, pre[px+1], unit, sep=""))
}


## * read args/input
ARGPARSE_TRUE <- TRUE
parser <- ArgumentParser()

## ** dryrun some args - allows interface with bash master script
args.cmd <- commandArgs(TRUE)
if (args.cmd[1] == "--argparse"){
    ARGPARSE_TRUE <- "" # FALSE does not work - argparse bug
    args.cmd <- args.cmd[-1]
}

## ** argparse

## -s (hort), --long ...
parser$add_argument("-i", required=ARGPARSE_TRUE, metavar="PAF", help="supported formats: paf")
parser$add_argument("-l", required=ARGPARSE_TRUE, metavar="LEN", help="per set sequence lengths")
parser$add_argument("-o", metavar="OUT", default="minidot.pdf", help="output file, .pdf/.png")
parser$add_argument("-S", "--no-self", action="store_true", default=FALSE, help="exclude plots of set against itself")
parser$add_argument("--title", help="plot title")
parser$add_argument("--theme", default="dark", help="themes: dark, light. [dark]")
parser$add_argument("--width", default=20, help="plot width (cm)", type="integer")


args <- parser$parse_args(args.cmd)

if (ARGPARSE_TRUE == "") quit();

## debug
#setwd("/home/thackl/projects/coding/sandbox/R-minimap-dotplots")
#args <- c("minidot.paf", "minidot.len");

## * read paf and len
paf <- read.table(args$i)
len <- read.table(args$l) #, stringAsFactor=FALSE)

sets <- unique(len$V1)

len$V1 <- factor(len$V1, levels=sets)
paf$V1 <- factor(paf$V1, levels=sets)
paf$V6 <- factor(paf$V6, levels=sets)

sets.n <- length(unique(paf$V1))

if(args$no_self){
    if (sets.n > 1){
        paf<-paf[!paf$V1==paf$V6,]
        if(sets.n == 2){
            ## make sure that for two sets, they don't appear in same dimension
            paf<-paf[paf$V1==paf$V1[1],]
        }
    }
}

if(dim(paf)[1]==0){
    write("no matches between sets, nothing to plot", stderr())
    quit(status=1);
}

paf$ava <- paf$V1:paf$V6
paf$strand <- ifelse(paf$V5=='+', 1, -1)
paf[paf$strand==-1,8:9] <- paf[paf$strand==-1,9:8]
paf$idy <- paf$V10 / paf$V11 * paf$strand

## * map contig boundaries to gglayer

len.cum <- cbind(len, cumsum=c(lapply(split(len, len$V1),
               function(x) cumsum(x$V2)), recursive=T))
yava <- data.frame(V1=character(0), V6=character(0), yi=numeric(0))
ava.se <- data.frame(V1=character(0), V6=character(0), x=numeric(0), xend=numeric(0), y=numeric(0), yend=numeric(0))
xava <- data.frame(V1=character(0), V6=character(0), xi=numeric(0))
yava.rt <- data.frame(V1=character(0), V6=character(0), xmin=numeric(0), xmax=numeric(0), ymin=numeric(0), ymax=numeric(0))
xava.rt <- data.frame(V1=character(0), V6=character(0), xmin=numeric(0), xmax=numeric(0), ymin=numeric(0), ymax=numeric(0))

ava.bg <- data.frame(V1=character(0), V6=character(0), xmin=numeric(0), xmax=numeric(0), ymin=numeric(0), ymax=numeric(0))

for(i in unique(paf$ava)){
    r <- str_split_fixed(i, ":", 2)
    yi <- c(0,len.cum$cumsum[len.cum$V1==r[2]])
    yi.max <- yi[length(yi)]
    xi <- c(0,len.cum$cumsum[len.cum$V1==r[1]])
    xi.max <- xi[length(xi)]

    ava.bg <- rbind(ava.bg, data.frame(V1=r[1], V6=r[2],
                         xmin=xi.max * -0.01, xmax=xi.max*1.01,
                         ymin=yi.max * -0.01, ymax=yi.max*1.01))

    ava.se <- rbind(ava.se, data.frame(V1=r[1], V6=r[2],  y=yi, yend=yi, x=0, xend=xi.max))
    ava.se <- rbind(ava.se, data.frame(V1=r[1], V6=r[2],  y=0, yend=yi.max, x=xi, xend=xi))
}

## * plot

gg <- ggplot(paf)
#gg <- gg + geom_rect(data=xava.rt, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="green", alpha=0.5)
#gg <- gg + geom_rect(data=yava.rt, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="red", alpha=0.5)

if (args$theme=="dark"){
    col.fill <- "grey20"
    col.line <- "grey50"
}else if (args$theme=="light"){
    col.fill <- "white"
    col.line <- "grey20"
}

gg <- gg + geom_rect(data=ava.bg, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill=col.fill)
gg <- gg + geom_segment(data=ava.se, aes(x=x,xend=xend,y=y,yend=yend), size=.1, color=col.line, linetype = 3)
gg <- gg + geom_segment(data=paf, aes(x=V3, xend=V4, y=V8, yend=V9, color=idy), size=.4, lineend = "round")

if (args$theme=="dark") gg <- gg + scale_colour_distiller("Identity", palette="Spectral", direction=1)
if (args$theme=="light") gg <- gg + scale_colour_gradientn("Identity", colours = c("#d60004", "#e8ae00", "#666666", "#666666", "#19bf5e", "#1701d2"))

gg <- gg + coord_fixed(1) + theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "grey92"),
    axis.title.x=element_blank(),
    axis.title.y=element_blank()
)

gg <- gg + scale_x_continuous(label=humanize_format(unit='',sep=''), expand=c(0,0))
gg <- gg + scale_y_continuous(label=humanize_format(unit='',sep=''), expand=c(0,0))
gg <- gg + facet_grid(V6~V1, drop=TRUE, as.table=FALSE)

if (!is.null(args$title)) gg <- gg + ggtitle(args$title)

# adjust width for direct compare of two set with different genome length
plot.height = args$width -2 # assume square plot, but leave some extra room for legend
if (sets.n == 2 && args$no_self) {
    x <- ava.bg$xmax
    y <- ava.bg$ymax
    bp.width <- args$width/x
    plot.height <- ava.bg$ymax * bp.width - ifelse(x > y, 0, 8)
}

ggsave(args$o, plot=gg, width=args$width, height=plot.height, units="cm")
