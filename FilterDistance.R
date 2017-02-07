setwd("~/Desktop/SDM")
PRES <- read.csv("Dipo_presence.csv", header = T)[, -1] #54 presence points
ABSV <- read.csv("Dipo_background.csv", header = T)[, -1] #background points developed from a bias file
PRED <- list.files(path = "/Volumes/CHAETO/SDM/bio_30s_bil", pattern = "bil", full.names = T) #30s Bioclim data
PRED <- stack(PRED)
MCP <- readOGR(dsn = path.expand("~/Desktop/SDM"), layer = "MCP") #minimum convex polygon of presence points
EXT <- extent(c(MCP@bbox[, 1] - 1, MCP@bbox[, 2] + 1)[c(1, 3, 2, 4)]) #add one degree to the polygon to create the extent


#Grabbing PCA scores of environmental variables for study extent 
SUBR <- crop(PRED, EXT)
#alt <- lonlatFromCell(SUBR[[1]], 1:ncell(SUBR))
##do we need geographic distance for anything (e.g., determining an appropriate distance to use for spatial filtering)? If so, this would give lat/long for each cell
VALR <- getValues(SUBR)
VALR <- na.omit(VALR)
PCA <- prcomp(VALR)
summary(PCA)
DISTE <- dist(PCA$x[, 1:2])

#Grabbing PCA scores of environmental variables for presence points
PCB <- prcomp(PRES[, -1:-2])
summary(PCB)
PRESR <- cbind(PRES, PCB$x[, 1:2])
DISTP <- dist(PRESR[, 22:23]) #euclidean distance for all presence points

###At this point, where do we go? We have distances in environmental space for each presence point and same for each cell in the study extent. The question is what is an appropriate distance to filter by and can we get it from the study extent in some way.