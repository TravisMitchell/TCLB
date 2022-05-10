current_path <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path ))

load("../output/Densities.RDATA")

dim(densities)

velocities.x = c(0,1,0,-1,0,1,-1,-1,1)
velocities.y = c(0,0,1,0,-1,1,1,-1,-1)

phasefield <- rowSums(densities[,,9:18],dims = 2)
image(phasefield)

# Collapse to 2d data set for python
data.out <- densities
dim(densities)
dim(data.out) <- c(31*2000,18)
head(data.out)
write.csv(x = data.out,file = "data_for_python.csv", row.names = FALSE)
