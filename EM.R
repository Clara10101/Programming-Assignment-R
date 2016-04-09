data <- c(rnorm(20,mean=1,sd=1),rnorm(20,mean=16,sd=1))

EM_algorithm <- function(data){
        u <- sample(data,2)
        #d <- c(var(data),var(data))
        pi <- 0.5
        for (i in 1:50){
                #Step E
                g <- pi*dnorm(data,mean=u[2])/(pi*dnorm(data,mean=u[2]) + (1-pi)*dnorm(data,mean=u[1]))
                
                #Step M
                u[1] <- sum((1-g)*data)/sum(1-g)
                u[2] <- sum(g*data)/sum(g)
                #d[1] <- sum((1-g)*(data-u[1])^2)/sum(1-g)
                #d[2] <- sum((g*(data-u[2])^2))/sum(g)
                pi <- mean(g)
                print(c(u,pi))
        }
}

EM_algorithm(data)
