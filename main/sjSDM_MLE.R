library(sjSDM)
set.seed(42)
load("data/data.Rdata")
env <- scale(temp)

env <- as.matrix(env)
colnames(env)[1] <- "V1"

ixp <- as.matrix(ixp)
X = model.matrix(~poly(V1, degree = 2 ), data.frame(env))
Y = ixp

sampling = 4000L
model = sjSDM(Y = ixp, env = linear(data= env, formula = ~poly(V1, degree = 2 )), biotic = bioticStruct(df = 6L), 
              se = FALSE, iter = 600L, family=multinomial(), sampling = 4000L, step_size = 1L, device = 0)

saveRDS(model, "results/mle_sjSDM.RDS")
