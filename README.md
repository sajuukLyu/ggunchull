# gg-unc-hull
A ggplot extension for drawing smooth non-convex irregular hulls around groups of points, which can fit the shape of boundaries perfectly!

## Install

`ggunchull` can be installed from github.

```R
devtools::install_github("sajuukLyu/ggunchull", type = "source")
```

## Usage

Take a pancreas dataset for example:

```R
library(ggplot2)
library(ggunchull)

data("pancExample")
plotData <- pancExample

# setup a proper distance to extend (and simplify) the circle (here use 3% of the range of x axis)
delta_r <- 0.03
th_r <- 0.03
delta <- diff(range(plotData$UMAP_1)) * delta_r
th <- diff(range(plotData$UMAP_1)) * th_r

ggplot(plotData, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = celltype), size = 0.3, shape = 16) +
  stat_unchull(aes(color = celltype, fill = celltype), alpha = 0.3, n = 5, delta = delta, th = th) +
  theme(
    aspect.ratio = 1,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line()
  )
```

<img src="plot\circleClusters.png"/>

The detailed algorithm will be explained later soon... (when I have free time)

