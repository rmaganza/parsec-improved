library(parsec)
source("/Users/XelmagaX/Desktop/UNI/Tesi/code/computePolarizationMeasure_improved.r")
source("/Users/XelmagaX/Desktop/UNI/Tesi/code/runC.pre_parsec_polarization.r")
source("/Users/XelmagaX/Desktop/UNI/Tesi/code/summary.parsec_polarization.r")
dyn.load("/Users/XelmagaX/Desktop/UNI/Tesi/code/working/parsec_improved.so")

congo <- read.table("/Users/XelmagaX/Desktop/UNI/Tesi/code/data/CONGO.txt", header=T)

profiles <- congo[,1:4]
frequencies <- congo[, 5:16]
congo[, 5:16] <- sweep(frequencies,2, colSums(frequencies), "/")

profileslist <- lapply(5:16, function(i) pop2prof(profiles, weights = congo[,i]))
names(profileslist) <- names(congo[,5:16])

lambda <- getlambda(Water > Sanitation, Health > Sanitation, Sanitation > Shelter)

reduceLE <- function(x,y) {
  ord <- rownames(x)
  return(x * y[ord, ord])
}
lst <- LE(lambda)
lstZeta <- LE2incidence(lst, varmod=list(Shelter=0:1, Sanitation=0:1, Water=0:1, Health=0:1))
zeta <- Reduce(reduceLE, lstZeta)
class(zeta) <- "incidence"

kss <- computePolarizationMeasure(profileslist$KSS, measure = "leti", lambda = lambda)
measures <- lapply(profileslist, computePolarizationMeasure, measure="leti", lambda=lambda)

lambda2 <- getlambda(Health,Water,Shelter,Sanitation)
measures2 <- lapply(profileslist, computePolarizationMeasure, measure="leti", lambda=lambda2)
