#zadanie zaliczeniowe
setwd("C:/Users/Klara/Desktop/projekt sad")
#Reads a file in table format and creates a data frame from it, with cases corresponding to lines and variables to fields in the file.
expression_data_24h <- read.csv(file="PJ34_C2_expression_data_24h.csv",header=TRUE,sep=",")
expression_data_6h <- read.csv(file="PJ34_C2_expression_data_6h.csv",header=TRUE,sep=",")
expression_data_3h <- read.csv(file="PJ34_C2_expression_data_3h.csv",header=TRUE,sep=",")

#ile genow analizowanych w kazdym punkcie czasowym (unikalne symbole wierszy w plikach z danymi)
length(unique(as.vector(expression_data_24h[,1])))
length(unique(as.vector(expression_data_6h[,1])))
length(unique(as.vector(expression_data_3h[,1])))

dim(expression_data_24h)
dim(expression_data_6h)
dim(expression_data_3h)

#sprawdzenie czy kolumny symbol sa takie same dla wszystkich powtorzen w czasie
all(expression_data_24h[,1] == expression_data_6h[,1])
all(expression_data_24h[,1] == expression_data_3h[,1])

#jedna duza data frame
full_data <- cbind(expression_data_3h,expression_data_6h[,-c(1,2)],expression_data_24h[,-c(1,2)])

all_values <- c()
for(i in 3:38){
        all_values <- c(all_values,full_data[,i])
}

ggplot() + aes(all_values)+ geom_histogram(binwidth=0.3, colour="black", fill="white",alpha=0.75) +
labs(title="Histogram wszystkich wyników", x="poziom ekspresji", y="liczba wyników")

#podsumowanie
summary(all_values)
var(all_values)
sd(all_values)

#statystyki dla ka¿dej kolumny
for(i in 3:38){
        print("kolumna o nazwie")
        print(names(full_data)[i])
        print(summary(full_data[,i]))
        print(var(full_data[,i]))
        print(sd(full_data[,i]))
}

len <- length(full_data[,1])

#box plot podsumowujacy
x_1 <- data.frame(
        values = all_values,
        powtorzenia = rep(c(rep(c('1'), len),rep(c('2'), len),rep(c('3'), len)),12),
        proba = rep(c(rep(c('C2'),len*3),rep(c('K'),len*3),rep(c('PJ34'),len*3),rep(c('PJ34C2'),len*3)),3),
        punkty_czasowe = c(rep(c('3h'),len*12),rep(c('6h'),len*12),rep(c('24h'),len*12))
)

ggplot(x_1, aes(x = punkty_czasowe, y = values, fill = powtorzenia)) +
        geom_boxplot() +
        facet_wrap(~ proba) + labs(x="wartoœci", y="punkty czasowe")

#wybor genow do dalszej analizy
#t test podzial na dwie grupy - kontrola i pozostale

#wyznaczenie genow do dalszej analizy osobno dla roznych srodowisk

#indeksy kolumn z kontrola 3h
kolumny_kontrola <- c(6,7,8,18,19,20,30,31,32)

#C2
kolumna_c2 <- c(3,4,5,15,16,17)
pvals_c2 <- apply(full_data,1,function(x) {t.test(as.numeric(x[kolumny_kontrola]),as.numeric(x[kolumna_c2]))$p.value})
sum(pvals_c2 < 0.05)
#poprawka zwiazana z testowanyem wielu hipotez
pad_c2 <- p.adjust(pvals_c2,method="fdr")
sum(pad_c2 < 0.05)
geny_c2 <- full_data$SYMBOL[pad_c2 < 0.05]

#PJ34
kolumna_pj34 <- c(9,10,11,21,22,23,33,34,35)
pvals_pj34 <- apply(full_data,1,function(x) {t.test(as.numeric(x[kolumny_kontrola]),as.numeric(x[kolumna_pj34]))$p.value})
sum(pvals_pj34 < 0.05)
pad_pj34 <- p.adjust(pvals_pj34,method="fdr")
sum(pad_pj34 < 0.05)
geny_pj34 <- full_data$SYMBOL[pad_pj34 < 0.05]

#PJ34C2
kolumna_pj34c2 <- c(12,13,14,24,25,26,36,37,38)
pvals_pj34c2 <- apply(full_data,1,function(x) {t.test(as.numeric(x[kolumny_kontrola]),as.numeric(x[kolumna_pj34c2]))$p.value})
sum(pvals_pj34c2 < 0.05)
pad_pj34c2 <- p.adjust(pvals_pj34c2,method="fdr")
sum(pad_pj34c2 < 0.05)
geny_pj34c2 <- full_data$SYMBOL[pad_pj34c2 < 0.05]

all_genes <- c(as.character(geny_c2), as.character(geny_pj34), as.character(geny_pj34c2))
#geny wybrane do dalszej analizy
length(unique(all_genes))

#wektor odpowiadajacy analizowanym genom
geny <- pad_pj34 < 0.05 | pad_pj34c2 < 0.05 | pad_c2 < 0.05

#zakladam model liniowy, gdzie poziom ekspresji zalezy od zmiennych:
#obecnosc_c2
#obecnosc_pj34
#liczba_godzin
#czyli postaæ: y = a * obecnosc_c2 + b * obecnosc_pj34 + c * liczba_godzin + d + N(0,sd)
y <- c()
apply(full_data[,-c(1,2)],2,function(x){y <<- c(y,x[geny])})
#length(y) = 49788 czyli 36*1383
sum_geny <- sum(geny)
#trzy punkty czasowe
obecnosc_c2 <- rep(c(rep(1,sum_geny*3),rep(0,sum_geny*6),rep(1,sum_geny*3)),3)
obecnosc_pj34 <- rep(c(rep(0,sum_geny*6),rep(1,sum_geny*6)),3)
liczba_godzin <- c(rep(3,sum_geny*12),rep(6,sum_geny*12),rep(24,sum_geny*12))

length(y) == length(obecnosc_c2)
length(y) == length(obecnosc_pj34)
length(y) == length(liczba_godzin)

#zrodlo kodu https://theoreticalecology.wordpress.com/2010/09/17/metropolis-hastings-mcmc-in-r/
likelihood <- function(param){
        a = param[1]
        b = param[2]
        c = param[3]
        d = param[4]
        sd = param[5]
        
        pred = a * obecnosc_c2 + b * obecnosc_pj34 + c * liczba_godzin + d
        singlelikelihoods = dnorm(y, mean = pred, sd = sd, log = T)
        sumll = sum(singlelikelihoods)
        return(sumll)   
}

# Prior distribution
prior <- function(param){
        a = param[1]
        b = param[2]
        c = param[3]
        d = param[4]
        sd = param[5]
        
        aprior = dunif(a, min=0, max=20, log = T)
        bprior = dunif(b, min=0, max=20, log = T)
        cprior = dunif(b, min=3, max=24, log = T)
        dprior = dnorm(b, sd = 5, log = T)
        sdprior = dunif(sd, min=0, max=30, log = T)
        return(aprior+bprior+sdprior)
}

posterior <- function(param){
        return (likelihood(param) + prior(param))
}

######## Metropolis algorithm ################

proposalfunction <- function(param){
        return(rnorm(5,mean = param, sd= c(0.1,0.5,0.3,0.2,0.3)))
}

run_metropolis_MCMC <- function(startvalue, iterations){
        chain = array(dim = c(iterations+1,5))
        chain[1,] = startvalue
        for (i in 1:iterations){
                proposal = proposalfunction(chain[i,])
                
                probab = exp(posterior(proposal) - posterior(chain[i,]))
                if (runif(1) < probab){
                        chain[i+1,] = proposal
                }else{
                        chain[i+1,] = chain[i,]
                }
        }
        return(chain)
}

startvalue = c(2,2,1,0,10)
chain = run_metropolis_MCMC(startvalue, 10000)

burnIn = 5000
acceptance = 1-mean(duplicated(chain[-(1:burnIn),]))

#wyznaczone wartosci parametrow
mean(chain[-(1:burnIn),1])
mean(chain[-(1:burnIn),2])
mean(chain[-(1:burnIn),3])

#wczytanie identyfikatorow genow zwiazanych ze smiercia komorki
dat = readLines("idpantherFromUniprot.txt")
#ile takich genow w analizowanych danych
sum(apply(full_data[1],1,function(x){x %in% dat}))
#512
d <- apply(full_data[1],1,function(x){x %in% dat})

#model liniowy tylko dla genow zwiazanych ze smiercia komorki
geny_d <- pad_pj34 < 0.05 | pad_pj34c2 < 0.05 | pad_c2 < 0.05 | d
sum_geny_d <- sum(geny_d)

y_2 <- c()
apply(full_data[,-c(1,2)],2,function(x){y_2 <<- c(y_2,x[geny_d])})

#Alanogiczna analiza dla genow zwiazanych ze smiercia komorki
#trzy punkty czasowe
obecnosc_c2 <- rep(c(rep(1,sum_geny_d*3),rep(0,sum_geny_d*6),rep(1,sum_geny_d*3)),3)
obecnosc_pj34 <- rep(c(rep(0,sum_geny_d*6),rep(1,sum_geny_d*6)),3)
liczba_godzin <- c(rep(3,sum_geny_d*12),rep(6,sum_geny_d*12),rep(24,sum_geny_d*12))

length(y_2) == length(obecnosc_c2)
length(y_2) == length(obecnosc_pj34)
length(y_2) == length(liczba_godzin)

likelihood <- function(param){
        a = param[1]
        b = param[2]
        c = param[3]
        d = param[4]
        sd = param[5]
        
        pred = a * obecnosc_c2 + b * obecnosc_pj34 + c * liczba_godzin + d
        singlelikelihoods = dnorm(y_2, mean = pred, sd = sd, log = T)
        sumll = sum(singlelikelihoods)
        return(sumll)   
}

# Prior distribution
prior <- function(param){
        a = param[1]
        b = param[2]
        c = param[3]
        d = param[4]
        sd = param[5]
        
        aprior = dunif(a, min=0, max=20, log = T)
        bprior = dunif(b, min=0, max=20, log = T)
        cprior = dunif(b, min=3, max=24, log = T)
        dprior = dnorm(b, sd = 5, log = T)
        sdprior = dunif(sd, min=0, max=30, log = T)
        return(aprior+bprior+sdprior)
}

posterior <- function(param){
        return (likelihood(param) + prior(param))
}

######## Metropolis algorithm ################

proposalfunction <- function(param){
        return(rnorm(5,mean = param, sd= c(0.1,0.5,0.3,0.2,0.3)))
}

run_metropolis_MCMC <- function(startvalue, iterations){
        chain = array(dim = c(iterations+1,5))
        chain[1,] = startvalue
        for (i in 1:iterations){
                proposal = proposalfunction(chain[i,])
                
                probab = exp(posterior(proposal) - posterior(chain[i,]))
                if (runif(1) < probab){
                        chain[i+1,] = proposal
                }else{
                        chain[i+1,] = chain[i,]
                }
        }
        return(chain)
}

startvalue = c(2,2,1,0,10)
chain = run_metropolis_MCMC(startvalue, 10000)

burnIn = 5000
acceptance = 1-mean(duplicated(chain[-(1:burnIn),]))

#wyznaczone wartosci parametrow
mean(chain[-(1:burnIn),1])
mean(chain[-(1:burnIn),2])
mean(chain[-(1:burnIn),3])
