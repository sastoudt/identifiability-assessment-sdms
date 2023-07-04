## theta labels done manually

### Figure 1 ###

p1 <- function(x, y) {
  -(x^2 + y^2)
}

x <- -10:10
y <- -10:10
z <- outer(x, y, p1)

### Figure A ###
x2 <- -20:10
jpeg("paper_figures/fig1a.jpeg", width = 824, height = 634, units = "px")
res <- persp(x, y, z, theta = 50, xlab = "x", ylab = "y", zlab = "z", border = "black", col = "white")
points(trans3d(0, y = 0, z = 0, pmat = res), col = "red", cex = 2, pch = 19)
dev.off()
### Figure C ###

x2 <- -20:10
jpeg("paper_figures/fig1c.jpeg", width = 824, height = 634, units = "px")
res <- persp(x, y, z, theta = 50, xlab = "x", ylab = "y", zlab = "z", border = "black", col = "white")
points(trans3d(0, y = 0, z = 0, pmat = res), col = "red", cex = 2, pch = 19)

lines(trans3d(x2, y = 5, z = -(25 + x2^2), pmat = res), col = "black", lwd = 4)
points(trans3d(0, y = 5, z = -25, pmat = res), col = "black", cex = 2, pch = 19)
dev.off()

### Figure E ###

jpeg("paper_figures/fig1e.jpeg", width = 824, height = 634, units = "px")
contour(x, y, z)
lines(x, rep(5, length(x)), lwd = 2)
points(0, 0, col = "red", pch = 19, cex = 2)
points(0, 5, col = "black", pch = 19, cex = 2)
dev.off()

### Figure B ###
p2 <- function(x, y) {
  -(x^2)
}

x <- -10:10
y <- -10:10
z <- outer(x, y, p2)

jpeg("paper_figures/fig1b.jpeg", width = 824, height = 634, units = "px")
res <- persp(x, y, z, theta = 50, xlab = "x", ylab = "y", zlab = "z", border = "black", col = "white")
points(trans3d(0, y = -7, z = 0, pmat = res), col = "red", cex = 2, pch = 19)
dev.off()

### Figure D ###

jpeg("paper_figures/fig1d.jpeg", width = 824, height = 634, units = "px")
res <- persp(x, y, z, theta = 50, xlab = "x", ylab = "y", zlab = "z", border = "black", col = "white")
points(trans3d(0, y = -7, z = 0, pmat = res), col = "red", cex = 2, pch = 19)
lines(trans3d(x, 5, z = -x^2, pmat = res), col = "black", lwd = 4)
points(trans3d(0, y = 5, z = 0, pmat = res), col = "black", cex = 2, pch = 19)
dev.off()

### Figure F ###

jpeg("paper_figures/fig1f.jpeg", width = 824, height = 634, units = "px")
contour(x, y, z)
lines(x, rep(5, length(x)), lwd = 2)
points(0, -7, col = "red", pch = 19, cex = 2)
points(0, 5, col = "black", pch = 19, cex = 2)
dev.off()
