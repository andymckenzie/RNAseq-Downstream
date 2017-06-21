
library(bnlearn)

a = runif(1000, 0, 1)
b = runif(1000, 0, 1) + a
c = runif(1000, 0, 1) + b

plot(a, b)
plot(b, c)
plot(a, c)

#plot the conditional independence relationship
c1 = c - b
plot(a, c1)

ci.test(a, b, c, test = "mc-cor")
ci.test(b, c, a, test = "mc-cor")
ci.test(a, c, b, test = "mc-cor")
