library(knitr)
opts_chunk$set(cache=TRUE)
library(ggplot2)
library(grid)
library(plyr)
library(RColorBrewer)
library(truncnorm)
library(dplyr)
library(reshape2)
library(compiler)
setCompilerOptions(optimize=3)
NagelkerkeR2 <- function (rr) 
{
    n <- nrow(rr$model)
    R2 <- (1 - exp((rr$dev - rr$null)/n))/(1 - exp(-rr$null/n))
    RVAL <- list(N = n, R2 = R2)
    return(RVAL)
}
result_h2l <- function(k, r2n, p) {
    x <- qnorm(1 - k)
    z <- dnorm(x)
    i <- z / k
    cc <- k * (1 - k) * k * (1 - k) / (z^2 * p * (1 - p))
    theta <- i * ((p - k)/(1 - k)) * (i * ((p - k) / ( 1 - k)) - x)
    e <- 1 - p^(2 * p) * (1 - p)^(2 * (1 - p))
    h2l <- cc * e * r2n / (1 + cc * e * theta * r2n)
    return(h2l)
}


# Variable initialization -------------------------------------------------

nselect <- 5000
sample.size <- 10000000
r2 <- seq(0.01,0.5,0.02)
prevalence <- seq(0.01,0.35,0.05)
sampling <- c("Normal", "Healthy", "Disease")
control.ratio <- seq(1.0,2,0.5)
perm <- 1

simulated.result <- expand.grid(Sampling=sampling, Prevalence=prevalence, R2=r2, Ratio=control.ratio, Perm=1:perm, MeanSA=NA, MeanSU=NA,VarSA=NA,VarSU=NA,MeanPA=NA, MeanPU=NA,VarPA=NA,VarPU=NA,MeanSA3=NA, MeanSU3=NA,VarSA3=NA,VarSU3=NA, obsR2=NA, leeR2=NA, pearR2a=NA, pearR2b=NA)

for(p in 1:perm)
{
    selectPerm <- simulated.result$Perm %in% p
    for(i in r2)
    {
        pheno <- rnorm(sample.size, mean=0,sd=1)
        score <- pheno*sqrt(i)+rnorm(sample.size, mean=0,sd=sqrt(1-i))
        score <- score/sd(score)*rnorm(1)+rnorm(1)
        #score2 <- pheno*sqrt(i)+rnorm(sample.size, mean=0,sd=sqrt(1-i))
        #score2 <- score2/sd(score2)*rnorm(1)+rnorm(1)
        #score3 <- pheno*sqrt(i)+rnorm(sample.size, mean=0,sd=sqrt(1-i))
        #score3 <- score3/sd(score3)*rnorm(1)+rnorm(1)
        info <- data.frame(score,pheno) %>% arrange(pheno)
        selectR2 <- simulated.result$R2 %in% i & selectPerm
        for(k in prevalence)
        {
            selectK <- selectR2 & simulated.result$Prevalence%in%k
            threshold <- qnorm(k, lower.tail = F)
            info$cc <- as.numeric(info$pheno > threshold)
            case <- info$cc == 1
            control <- !case
            for(r in control.ratio)
            {
                cat(paste("\r",p,i,k,r))
                selectR <- selectK & simulated.result$Ratio %in% r
                for(s in sampling)
                {
                    controls <- sort(sample(which(control), nselect*r))
                    ncase <- sum(case)
                    if(s=="Normal")
                    {
                        cases <- sort(sample(which(case), nselect))
                    }else if(s=="Healthy")
                    {
                        cases <- sort(sample(head(which(case), 0.2*ncase), nselect))
                    }else if(s=="Disease"){
                        cases <- sort(sample(tail(which(case), 0.2*ncase), nselect))
                    }
                    obsR2 <- NagelkerkeR2(suppressWarnings(glm(info$cc[c(controls,cases)]~info$score[c(controls,cases)], family=binomial)))$R2
                    
                    
                    selected <- selectR &simulated.result$Sampling %in% s
                    simulated.result$MeanSA[selected] <- mean(info$score[cases])
                    simulated.result$MeanSU[selected] <- mean(info$score[controls])
                    simulated.result$VarSA[selected] <- var(info$score[cases])
                    simulated.result$VarSU[selected] <- var(info$score[controls])
                    simulated.result$MeanPA[selected] <- mean(info$pheno[cases])
                    simulated.result$MeanPU[selected] <- mean(info$pheno[controls])
                    simulated.result$VarPA[selected] <- var(info$pheno[cases])
                    simulated.result$VarPU[selected] <- var(info$pheno[controls])
                    
                    #simulated.result$MeanSA2[selected] <- mean(info$score2[cases])
                    #simulated.result$MeanSU2[selected] <- mean(info$score2[controls])
                    #simulated.result$VarSA2[selected] <- var(info$score2[cases])
                    #simulated.result$VarSU2[selected] <- var(info$score2[controls])
                    
                    #simulated.result$MeanSA3[selected] <- mean(info$score3[cases])
                    #simulated.result$MeanSU3[selected] <- mean(info$score3[controls])
                    #simulated.result$VarSA3[selected] <- var(info$score3[cases])
                    #simulated.result$VarSU3[selected] <- var(info$score3[controls])
                    
                    simulated.result$obsR2[selected] <- obsR2
                    simulated.result$leeR2[selected] <- result_h2l(k, obsR2, nselect/(nselect*r+nselect))
                    
                    pearR <-  matrix(c(var(info$score[controls]),-1,1-vtruncnorm(b=qnorm(k,lower.tail=F))), ncol=1) %>%polyroot
                    simulated.result$pearR2a[selected] <- pearR[1]
                    simulated.result$pearR2b[selected] <- pearR[2]
                }
            }
        }
    }
}

sim.result <- melt(simulated.result[,-c(4:8,12)], id=c("Sampling", "Prevalence","R2"))
sim.result$value <- Re(sim.result$value)
sim.result.sum <- sim.result %>%
    group_by(Sampling, Prevalence, R2, variable) %>%
    summarise(mean=mean(value), sd=sd(value),mse=mean((value-R2)^2))
levels(sim.result$variable) <- c("Observed", "Lee's Adjustment", "PA's Adjustment")
levels(sim.result$Sampling)<- c("Normal", "Healthy Bias", "Berkson's bias")
mean.full <- ggplot(sim.result)+
    geom_abline()+
    theme_bw()+
    geom_point(aes(x=R2, y=mean, color=Prevalence))+
    scale_color_distiller(palette = "Spectral")+
    facet_grid(Sampling~variable)+
    xlab("Empirical R2")+
    ylab("Mean Estimated R2")
sd.full <- ggplot(sim.result)+
    geom_abline()+
    theme_bw()+
    geom_point(aes(x=R2, y=sd, color=Prevalence))+
    scale_color_distiller(palette = "Spectral")+
    facet_grid(Sampling~variable)+
    xlab("Empirical R2")+
    ylab("SD of Estimated R2")
mse.full <- ggplot(sim.result)+
    geom_abline()+
    theme_bw()+
    geom_point(aes(x=R2, y=mse, color=Prevalence))+
    scale_color_distiller(palette = "Spectral")+
    facet_grid(Sampling~variable)+
    xlab("Empirical R2")+
    ylab("MSE Estimated R2")

mean.sub <- ggplot(subset(sim.result, R2<=0.3))+
    geom_abline()+
    theme_bw()+
    geom_point(aes(x=R2, y=mean, color=Prevalence))+
    scale_color_distiller(palette = "Spectral")+
    facet_grid(Sampling~variable)+
    xlab("Empirical R2")+
    ylab("Mean Estimated R2")
sd.sub <- ggplot(subset(sim.result, R2<=0.3))+
    geom_abline()+
    theme_bw()+
    geom_point(aes(x=R2, y=sd, color=Prevalence))+
    scale_color_distiller(palette = "Spectral")+
    facet_grid(Sampling~variable)+
    xlab("Empirical R2")+
    ylab("SD of Estimated R2")
mse.sub <- ggplot(subset(sim.result, R2<=0.3))+
    geom_abline()+
    theme_bw()+
    geom_point(aes(x=R2, y=mse, color=Prevalence))+
    scale_color_distiller(palette = "Spectral")+
    facet_grid(Sampling~variable)+
    xlab("Empirical R2")+
    ylab("MSE Estimated R2")

ggsave("Full.mean.png", mean.full)
ggsave("Full.sd.png", sd.full)
ggsave("Full.mse.png", mse.full)
ggsave("Sub.mean.png", mean.sub)
ggsave("Sub.sd.png", sd.sub)
ggsave("Sub.mse.png", mse.sub)


simulated.result$cov <- (simulated.result$MeanSA-simulated.result$MeanSU)/(etruncnorm(a=qnorm(simulated.result$Prevalence, lower.tail=F), b=qnorm(simulated.result$Prevalence+simulated.result$Prevalence*0.5))- etruncnorm(b=qnorm(simulated.result$Prevalence, lower.tail=F)))

simulated.result$varS <- simulated.result$VarSU+simulated.result$cov^2*(1- vtruncnorm(b=qnorm(simulated.result$Prevalence, lower.tail=F)))
