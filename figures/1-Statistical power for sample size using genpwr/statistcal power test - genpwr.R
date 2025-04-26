library(genpwr)

case = 9903
control = 2099293
total_N = case + control
ratio = case/total_N

pw <- genpwr.calc(calc = "power", model = "logistic", ge.interaction = NULL,
                  N=total_N, Case.Rate=ratio, k=NULL,
                  MAF=seq(0.01, 0.89, 0.02), OR=c(1.4),Alpha=5e-8,
                  True.Model='Additive', 
                  Test.Model="All")

head(pw)

power.plot(pw)
