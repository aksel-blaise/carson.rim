# load
devtools::install_github("mlcollyer/RRPP")
devtools::install_github("geomorphR/geomorph",
                         ref = "Stable", 
                         build_vignettes = TRUE)

library(here)
library(StereoMorph)
library(geomorph)
library(ggplot2)
library(wesanderson)

# digitize dataset
digitizeImage(image.file = 'rims', 
              shapes.file = 'shapes',
              landmarks.ref = 'lm.txt', 
              curves.ref = 'curve.txt')

# read shape data and define number of sLMs
shapes <- readShapes("shapes")
shapesGM <- readland.shapes(shapes, nCurvePts = 20)

# read qualitative data
qdata <- read.csv("qdata.csv",
                  header = TRUE,
                  row.names = 1)

# generalised procrustes analysis
Y.gpa <- gpagen(shapesGM, print.progress = FALSE)
plot(Y.gpa)

# dataframe
gdf <- geomorph.data.frame(shape = Y.gpa$coords,
                           size = Y.gpa$Csize,
                           surf.treat = qdata$surf.treat,
                           cont = qdata$ context)

# add centroid size to qdata
qdata$csz <- Y.gpa$Csize

# attributes
csz <- qdata$csz
surf.treat <- qdata$surf.treat
cont <- qdata$context

# palette
pal = wes_palette("Moonrise2")

# boxplot of Perdiz arrow points by context
csz.temp <- ggplot(qdata, aes(x = context, y = csz, color = context)) +
  geom_boxplot() +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.3) +
  scale_color_manual(values = pal) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = 'Raw Material', y = 'Centroid Size')

# render plot
csz.temp

# pca
pca <- gm.prcomp(Y.gpa$coords)
summary(pca)

# set plot parameters
pch.gps <- c(15,17)[as.factor(qdata$context)]
col.gps <- pal[as.factor(qdata$context)]
col.hull <- c("#798E87","#C27D38")

# pca plot
pc.plot <- plot(pca,
                asp = 1,
                pch = pch.gps,
                col = col.gps)
shapeHulls(pc.plot,
           groups = qdata$surf.treat,
           group.cols = col.hull)

# general allometry
fit.size <- procD.lm(shape ~ size,
                     data = gdf,
                     print.progress = FALSE,
                     iter = 9999)

# general allometry
anova(fit.size)

# size
fit.sz.reg <- procD.lm(size ~ cont,
                       data = gdf,
                       print.progress = FALSE,
                       iter = 9999)

# size
anova(fit.sz.reg)

# shape
fit.sh.reg <- procD.lm(shape ~ cont,
                       data = gdf,
                       print.progress = FALSE,
                       iter = 9999)

# shape
anova(fit.sh.reg)

land.gps <- c("A","A","B","B","B","A","A","A","A","A","A","A","A","A",
              "B","B","B","B","B","B","B","B","A","A","A","A","A","A",
              "A","A")
it <- integration.test(Y.gpa$coords,
                       partition.gp = land.gps,
                       iter = 9999,
                       print.progress = FALSE)

summary(it)
plot(it)

morphol.disparity(shape ~ cont,
                  groups = qdata$cont,
                  data = gdf,
                  print.progress = FALSE,
                  iter = 9999)

## viz ----
# find mean specimen
findMeanSpec(Y.gpa$coords)

# create outline object warped to the mean shape
warpRefOutline(556,)

# subset landmark coordinates to produce mean shapes
new.coords <- coords.subset(A = Y.gpa$coords,
                            group = qdata$context)
names(new.coords)

# group shape means
mean <- lapply(new.coords, mshape)

# plot mean shapes
plot(mean$`Pit 434`)
plot(mean$`Structure 23`)

# comparison plots
plotRefToTarget(mean$`Pit 434`,
                mean$`Structure 23`,
                method = "points",
                mag = 1)

plotRefToTarget(mean$burnished,
                mean$slipped,
                method = "points",
                mag = 1)

plotRefToTarget(mean$plain,
                mean$slipped,
                method = "points",
                mag = 1)
