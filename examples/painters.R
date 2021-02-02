library("MASS")
data(painters)
source("../cig.R")

X.painters = painters[, 1:4]
cig <- get_cig(X.painters, 0.05)
plot(cig$graph, layout=layout_with_kk, label.cex=.25, label.dist=10)
print(cig$r)
View(cig$mat)
