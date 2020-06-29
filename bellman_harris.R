#these libraries need to be loaded
library(utils)
library(rstan)
library(matrixStats)
library(ggplot2)
library(gridExtra)

d <- readRDS('data/south-korea.rds')
serial_interval <- readRDS('data/serial-interval.rds')
pad_serial.interval <- data.frame(
  "X" = (length(serial_interval$fit)+1):200,
  "fit" = rep(1e-17, 200)
)
serial_interval = rbind(serial_interval, pad_serial.interval)

stan_data<-list()
stan_data$N  <- nrow(d)
stan_data$cases <- d$cases
stan_data$N2 <- 60

stan_data$SI <- serial_interval[1:stan_data$N,2]

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
m <- stan_model(file.path('stan-models','bellman.stan'))

fit <- sampling(m,data=stan_data,iter=2000,warmup=1000,chains=5,thin=5,control = list(adapt_delta = 0.99, max_treedepth = 30))

out <- rstan::extract(fit)

data<-data.frame(
  x=d$t,
  y_real=stan_data$cases,
  y = colMeans(out$prediction),
  yl = colQuantiles(out$prediction,probs=c(0.025)),
  yh = colQuantiles(out$prediction,probs=c(0.975)),
  Rt = colMeans(out$Rt),
  Rtl = colQuantiles(out$Rt,probs=c(0.025)),
  Rth = colQuantiles(out$Rt,probs=c(0.975))
)
data2<-data.frame(
  x=d$t[1:stan_data$N2],
  ex = colMeans(out$mu)[1:stan_data$N2],
  exl = colQuantiles(out$mu,probs=c(0.025))[1:stan_data$N2],
  exh = colQuantiles(out$mu,probs=c(0.975))[1:stan_data$N2]
)

g1 <- ggplot(data, aes(x = x)) +
  geom_line(data,mapping=aes(y=y,x=x),col='deepskyblue4',size=1.2) +
  geom_bar(data = data, aes(y = y_real),
           fill = "coral4", stat='identity', alpha=0.5) +
  geom_ribbon(data,mapping=aes(x=x,ymin=yl,ymax=yh),fill='deepskyblue4',alpha=0.3) + 
  theme_classic() +  xlab("Days") + ylab('Infections')

g2 <- ggplot(data, aes(x = x)) +
  geom_line(data,mapping=aes(y=Rt,x=x),col='seagreen',size=1.2) +
  geom_ribbon(data,mapping=aes(x=x,ymin=Rtl,ymax=Rth),fill='seagreen',alpha=0.4) + 
  geom_hline(yintercept = 1,color = "black", size=1,alpha=0.5) +
  theme_classic() +  xlab("Days") + ylab('Reproduction Number') + scale_y_continuous(trans = "log10")

g3 <- ggplot(data2, aes(x = x)) +
  geom_line(data2,mapping=aes(y=ex,x=x),col='tan3',size=1.2) +
  geom_ribbon(data2,mapping=aes(x=x,ymin=exl,ymax=exh),fill='tan3',alpha=0.4) + 
  theme_classic() +  xlab("Days") + ylab('Imported Infections')

grid.arrange(g1,g2,g3)
