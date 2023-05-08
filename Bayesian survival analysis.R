set.seed(1)
library("rjags")
library("tidyverse")
library("tidybayes")
library("survminer")
library("survival")
library('bayesplot')
library("ggthemes")
library("coxed")

#simulate data. Ignore covariates. 
#data in reality is cut off at follow-up at the time of the interim analysis, but not in this example
#specify own hazard fun
my_hazard <- function(t){ 
  dnorm((log(t) - log(50))/log(40)) /
    (log(10)*t*(1 - pnorm((log(t) - log(50))/log(40))))
}

#interim analysis 1: 150 HO132 + 300 HO102 control treatment (matched with control treatment arm of HO132)
data_150 <- sim.survdata(N=450, T=10, xvars = 1, num.data.frames=1, hazard.fun = my_hazard)$data %>% select(y, failed) %>% 
  mutate(
  treatment = c(rep(0:1, length.out = 150), rep(0, length.out = 300)),
  trialnr = c(rep("HO132", length.out = 150), rep("HO102", length.out = 300)),
  failed = ifelse(y >= 24, FALSE, failed) # right censor all time>=24
)

#interim analysis 2: 300 HO132 + 300 HO102 control treatment (matched with control treatment arm of HO132)
data_300 <- sim.survdata(N=600, T=16, xvars = 1, num.data.frames=1, hazard.fun = my_hazard)$data %>% select(y, failed)%>% mutate(
  treatment = c(rep(0:1, length.out = 300), rep(0, length.out = 300)),
  trialnr = c(rep("HO132", length.out = 300), rep("HO102", length.out = 300)),
  failed = ifelse(y >= 24, FALSE, failed)
)

#interim analysis 3: 450 HO132 + 300 HO102 control treatment (matched with control treatment arm of HO132)
data_450 <- sim.survdata(N=750, T=22, xvars = 1, num.data.frames=1, hazard.fun = my_hazard)$data %>% select(y, failed)%>% mutate(
  treatment = c(rep(0:1, length.out = 450), rep(0, length.out = 300)),
  trialnr = c(rep("HO132", length.out = 450), rep("HO102", length.out = 300)),
  failed = ifelse(y >= 24, FALSE, failed)
)

#interim analysis 4: 600 HO132 + 300 HO102 control treatment (matched with control treatment arm of HO132)
data_600 <- sim.survdata(N=900, T=28, xvars = 1, num.data.frames=1, hazard.fun = my_hazard)$data %>% select(y, failed)%>% mutate(
  treatment = c(rep(0:1, length.out = 600), rep(0, length.out = 300)),
  trialnr = c(rep("HO132", length.out = 600), rep("HO102", length.out = 300)),
  failed = ifelse(y >= 24, FALSE, failed)
)

#get HR
cph <- coxph(Surv(y, failed) ~ treatment, data=data_150)
summary(cph)

#plot the curves of EFS
fit1 <- survfit(Surv(y, failed) ~ trialnr + treatment, data=data_150)

ggsurvplot(
  fit1
)

##############################################
#Bayesian interim analysis event free survival
##############################################
# The model specification
file.show("Bayesian_survival_model.txt")

#make a list of all dataframes
interim_list <- list(data_150 = data_150 %>% arrange(trialnr),
                     data_300 = data_300 %>% arrange(trialnr),
                     data_450 = data_450 %>% arrange(trialnr),
                     data_600 = data_600 %>% arrange(trialnr))

#create list with censored times indicator and time to event data, make time to event NA if censored
censored <- list()
is_censored <- list()
time <- list()
for (i in 1:length(interim_list)) {
censored[[i]] <- interim_list[[i]]["failed"] == FALSE
is_censored[[i]] <- censored[[i]]*1
time[[i]] <- interim_list[[i]]["y"]
time[[i]][which(time[[i]] == censored[[i]]),] <- NA
}

#also get the censored times for separate distribution in the JAGS model
t_cen <- list()
group <- list()
n1 <- list()
n2 <- list()
for (i in 1:length(censored)){
t_cen[[i]] <- interim_list[[i]]["y"]
t_cen[[i]][censored[[i]] == FALSE] <- 0
group[[i]] <- ifelse(interim_list[[i]]["treatment"]== 1, TRUE, FALSE)
n1[[i]] <- nrow(interim_list[[i]] %>% filter(trialnr == "HO102"))
n2[[i]] <- nrow(interim_list[[i]])
}

aml_data <- list()
mcmc_samples <- list()
results <- list()

for(i in 1:length(interim_list)){
  # put the data together for JAGS
  aml_data[[i]] <- list(time = time[[i]][,1],
                        t_cen = t_cen[[i]][,1],
                        is_censored = is_censored[[i]][,1],
                        n1 = n1[[i]],
                        n2 = n2[[i]],
                        group = group[[i]][,1],
                        #"Spike and Slab" hyperparameters
                        l.slab= 1.0E-4, #S_l (lower bound of slab) 
                        u.slab= 2, #S_u (upper bound of slab)
                        p0= 0.5, #p_0 (prior probability of being in the slab)
                        spike=5000 #entry is Spike for logit_theta1
  )
}



for(i in 1:length(interim_list)){
#which paramaters to monitor
params_to_monitor <- c("beta_0", "beta_1" , "hr_1", "p1", "p2", "p3")

# Running the model
model <- jags.model(file = "Bayesian_survival_model.txt", aml_data[[i]], 
                    n.chains = 3, n.adapt= 1000)
update(model, 1000); # Burning 1000 samples to the MCMC gods...
mcmc_samples[[i]] <- coda.samples(model, params_to_monitor, n.iter= 50000,
                                  inits = list(.RNG.name = "base::Wichmann-Hill",
                                               .RNG.seed = 2022,
                                               thin = 50))

results[[i]] <- mcmc_samples[[i]][[1]]
}

#check convergence
plot(results[[1]])
gelman.plot(mcmc_samples[[1]][,1])
plot(results[[2]])
gelman.plot(mcmc_samples[[2]][,1])
plot(results[[3]])
gelman.plot(mcmc_samples[[3]][,1])
plot(results[[4]])
gelman.plot(mcmc_samples[[4]][,1])

summarise_draws(mcmc_samples[[1]], mean, sd, ~quantile(.x, probs = c(0.5, 0.025,0.975))) #n = 300 
summarise_draws(mcmc_samples[[2]], mean, sd, ~quantile(.x, probs = c(0.5, 0.025,0.975))) #n = 450
summarise_draws(mcmc_samples[[3]], mean, sd, ~quantile(.x, probs = c(0.5, 0.025,0.975))) #n = 600
summarise_draws(mcmc_samples[[4]], mean, sd, ~quantile(.x, probs = c(0.5, 0.025,0.975))) #n = 600

#function for auc posterior
calc_auc_post <- function(data, threshold) {
  prob <- sum(data<threshold)/length(data)
  prb <- paste0(sapply(prob, round, 2)*100,"%")
  return(prb)
}

# results 
hr_1 <- list()
prob <- list()
for (hr in 1:length(results)){
  hr_1[[hr]] <- results[[hr]][,3]
  
  prob[[hr]] <- calc_auc_post(hr_1[[hr]], 0.76)
}

# Make some nice density plots
prob <- unlist(prob)
post <- data.frame(sample(hr_1[[1]], 1000), sample(hr_1[[2]], 1000), 
                   sample(hr_1[[3]], 1000), sample(hr_1[[4]], 1000))

colnames(post) <- c("Interim analysis 1", "Interim analysis 2", "Interim analysis 3", "Interim analysis 4")

dat_text <- data.frame(
  label = prob,
  key   = c("Interim analysis 1", "Interim analysis 2", "Interim analysis 3", "Interim analysis 4")
)

# plot posteriors with HDI
p <- post %>% 
  gather() %>% 
  
  mutate(key = factor(key)) %>% 
  
  ggplot() +
  geom_density(aes(x = value),
               fill = "grey92", adjust=2) +
  geom_vline(aes(xintercept = 0), linetype = "dashed") +
  stat_pointinterval(aes(x = value, y = 0),
                     point_interval = mode_hdi, .width = c(.95, .5)) +
  scale_y_continuous(NULL, breaks = NULL) +
  xlab(NULL) +
  facet_wrap(~key, scales = "free", ncol = 2)

#extract values for fill
to_fill <- data_frame(
  x = ggplot_build(p)$data[[1]]$x,
  y = ggplot_build(p)$data[[1]]$y,
  key = factor(ggplot_build(p)$data[[1]]$PANEL, levels = c(1,2,3,4), 
               labels = c("Interim analysis 1", "Interim analysis 2", "Interim analysis 3", "Interim analysis 4")))

post %>% 
  gather() %>% 
  
  mutate(key = factor(key)) %>% 
  
  ggplot() +
  geom_density(aes(x = value), trim = TRUE,
               fill = "#86D2ED", size = 0.8, adjust=2, alpha = 0.3) +
  geom_vline(aes(xintercept = 0.76), linetype = "dashed") +
  stat_pointinterval(aes(x = value, y = 0),
                     point_interval = mode_hdi, .width = c(.95, .5)) +
  scale_y_continuous(NULL, breaks = NULL) +
  scale_x_continuous(limits = c(0, 2), breaks = c(0,0.5,0.76,1,1.5,2)) + 
  xlab(NULL) + geom_text(
    data    = dat_text,
    mapping = aes(x = 0.75, y = 0.5, label = label),
    hjust = 1, vjust = -1.5) +
  geom_area(data = to_fill[to_fill$x <= 0.76, ], 
            aes(x=x, y=y), fill = "#0C2074",
            alpha = 0.3) +
  facet_wrap(~key, scales = "fixed", ncol = 1, strip.position = "left") + 
  theme_fivethirtyeight() + 
  theme(
    strip.text.y.left = element_text(angle = 90),
    plot.background = element_rect(fill = " white"),
    panel.background = element_rect(fill = "white"),
    legend.background = element_rect(fill = "white"),
    strip.background = element_rect(fill = "white"),
    axis.title.x =  element_text()
  ) +
  xlab("Hazard ratio")
