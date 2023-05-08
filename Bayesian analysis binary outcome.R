set.seed(123)
library("rjags")
library("tidyverse")
library("tidybayes")
library("ggthemes")
library("scales")

#data
simulated_outcome <- c(rbinom(300, 1, .80), rbinom(300, 1, .70), rbinom(300, 1, .80))

#interim analysis 1: 150 HO132 + 300 HO102 control treatment (matched with control treatment arm of HO132)
data_150 <- simulated_outcome[c(1:75, 301:375, 601:900)] %>% as_tibble() %>% mutate(
  treatment = c(rep(0, length.out = 75), rep(1, length.out = 75), rep(0, length.out = 300)),
  trialnr = c(rep("HO132", length.out = 150), rep("HO102", length.out = 300)),
  group = paste(treatment, trialnr)
)

#interim analysis 2: 300 HO132 + 300 HO102 control treatment (matched with control treatment arm of HO132)
data_300 <- simulated_outcome[c(1:150, 301:450, 601:900)] %>% as_tibble() %>% mutate(
  treatment = c(rep(0, length.out = 150), rep(1, length.out = 150), rep(0, length.out = 300)),
  trialnr = c(rep("HO132", length.out = 300), rep("HO102", length.out = 300)),
  group = paste(treatment, trialnr)
)

#interim analysis 3: 450 HO132 + 300 HO102 control treatment (matched with control treatment arm of HO132)
data_450 <- simulated_outcome[c(1:225, 301:525, 601:900)] %>% as_tibble() %>% mutate(
  treatment = c(rep(0, length.out = 225), rep(1, length.out = 225), rep(0, length.out = 300)),
  trialnr = c(rep("HO132", length.out = 450), rep("HO102", length.out = 300)),
  group = paste(treatment, trialnr)
)

#interim analysis 4: 600 HO132 + 300 HO102 control treatment (matched with control treatment arm of HO132)
data_600 <- simulated_outcome %>% as_tibble() %>% mutate(
  treatment = c(rep(0, length.out = 300), rep(1, length.out = 300), rep(0, length.out = 300)),
  trialnr = c(rep("HO132", length.out = 600), rep("HO102", length.out = 300)),
  group = paste(treatment, trialnr)
)

######################################
#Bayesian interim analysis
######################################
# The model specification
file.show("Bayesian_analysis_binary_outcome.txt")

#make a list of all dataframes
interim_list <- list(data_150 = data_150, 
                     data_300 = data_300,
                     data_450 = data_450,
                     data_600 = data_600)
success <- list()
g <- list()
y1 <- vector() # cr for historical data
y2 <- vector() # cr for control
y3 <- vector() # cr for treatment
for (i in 1:length(interim_list)) {
  success[[i]] <- with(interim_list[[i]], tapply(interim_list[[i]][["value"]], group, sum, na.rm = TRUE))
  g[[i]] <- interim_list[[i]] %>% count(group) %>% select(n)
}

for (j in 1:length(success)){
  y1[j] <- success[[j]][1]
  y2[j] <- success[[j]][2]
  y3[j] <- success[[j]][3]
}

h <- vector()# total number of historical treatment
n <- vector()# total number of control
m <- vector()# total number of treatment

for (k in 1:length(g)){
  h[k] <- g[[k]]$n[1]
  n[k] <- g[[k]]$n[2]
  m[k] <- g[[k]]$n[3]
}

data <- data.frame(
  y1 = y1, # cr for historical data
  y2 = y2, # cr for control
  y3 = y3, # cr for treatment
  n = n, # number of trials control
  m = m, # number of trials treatment
  h = h # number of trials historical control
)

data_list <- list()
results <- list()
mcmc_samples <- list()
for (l in 1:nrow(data)){
  data_list[[l]] <- list(
    y = c(data$y1[l], data$y2[l], data$y3[l], data$y2[l], data$y2[l]),
    n = data$n[l],
    m = data$m[l],
    h = data$h[l],
    #"Spike and Slab" hyperparameters
    l.slab= 1.0E-4, #S_l (lower bound of slab) 
    u.slab= 2, #S_u (upper bound of slab)
    p0= 0.5, #p_0 (prior probability of being in the slab)
    spike= 5000 #entry is Spike for logit_theta1
  )

  params_to_monitor <- c("diff", "theta2", "theta3", "p")
  
  # Running the model
  model <- jags.model(file = "Bayesian_analysis_binary_outcome.txt", data_list[[l]], n.chains = 3, n.adapt= 1000)
  update(model, 1000); # Burning 1000 samples to the MCMC gods...
  mcmc_samples[[l]] <- coda.samples(model, params_to_monitor, n.iter=50000,
                                    inits = list(.RNG.name = "base::Wichmann-Hill",
                                                 .RNG.seed = 2022))
  
  results[[l]] <- mcmc_samples[[l]][[1]]
}

#check for convergence
#plot(mcmc_samples[[1]]) #seems like convergence
#gelman.plot(mcmc_samples[[1]])

summarise_draws(mcmc_samples[[1]], mean, sd, ~quantile(.x, probs = c(0.5, 0.025,0.975))) #n = 150
summarise_draws(mcmc_samples[[2]], mean, sd, ~quantile(.x, probs = c(0.5, 0.025,0.975))) #n = 300
summarise_draws(mcmc_samples[[3]], mean, sd, ~quantile(.x, probs = c(0.5, 0.025,0.975))) #n = 450
summarise_draws(mcmc_samples[[4]], mean, sd, ~quantile(.x, probs = c(0.5, 0.025,0.975))) #n = 600


# results and treatment effect
theta_diff <- list()
prob <- list()
for (cr in 1:length(results)){
theta_diff[[cr]] <- results[[cr]][,1]
prob[[cr]] <- sum(theta_diff[[cr]]>0)/length(theta_diff[[cr]])
}

# Make some nice density plots
prob <- unlist(prob)
post <- data.frame(theta_diff[[1]], theta_diff[[2]],
theta_diff[[3]], theta_diff[[4]])

colnames(post) <- c("Interim analysis 1", "Interim analysis 2", 
                    "Interim analysis 3", "Interim analysis 4")
dat_text <- data.frame(
  label = paste0(sapply(prob[1:4], round, 2)*100, "%"),
  key   = c("Interim analysis 1", "Interim analysis 2", 
            "Interim analysis 3", "Interim analysis 4"),
  mean = sapply(theta_diff, mean)
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
  xlab(NULL) + geom_text(
    data    = dat_text,
    mapping = aes(x = 0, y = 0.5, label = label),
    hjust = 1, vjust = -0.5) +
  facet_wrap(~key, scales = "free", ncol = 1)

#extract values for fill
to_fill <- data_frame(
  x = ggplot_build(p)$data[[1]]$x,
  y = ggplot_build(p)$data[[1]]$y,
  key = factor(ggplot_build(p)$data[[1]]$PANEL, levels = c(1,2,3,4), 
               labels = c("Interim analysis 1", "Interim analysis 2", 
                          "Interim analysis 3", "Interim analysis 4")))

#plot
post %>% 
  gather() %>% 
  
  mutate(key = factor(key)) %>% 
  
  ggplot() +
  geom_density(aes(x = value), trim = TRUE,
               fill = "#86D2ED", size = 0.8, adjust=2, alpha = 0.3) +
  geom_vline(aes(xintercept = 0), linetype = "dashed") +
  stat_pointinterval(aes(x = value, y = 0),
                     point_interval = mode_hdi, .width = c(.95, .5)) +
  scale_y_continuous(NULL, breaks = NULL) + 
  scale_x_continuous(limits = c(-0.3,0.1), labels = c("-30", "-20", "-10", "0", "10")) +
  geom_text(
    data    = dat_text,
    mapping = aes(x = mean, y = 0, label = paste0(round(mean*100,1), "%")), vjust = -1) +
  geom_text(
    data    = dat_text,
    mapping = aes(x = 0, y = 0.5, label = label),
    hjust = -1, vjust = -1.5) +
  geom_area(data = to_fill[to_fill$x >= 0, ], 
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
  xlab("Treatment difference %")
