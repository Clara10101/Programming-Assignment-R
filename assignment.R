library(qvalue)
dane <- read.csv("http://bioputer.mimuw.edu.pl/~bartek/SAD-Z1/trial-6.csv",sep="\t")
#dim(dane) 200 106
zdrowi <- dane$severity==0
chorzy <- !zdrowi

f <- function(x){
      #zlicza ile z mutacja posrod zdrowych i chorych
      #interpretujac mut jako sukces jakie prawdop wylosowania tylu lub wiecej chorych z mutacja
      #sposrod z mutacja i bez mutacji losujac 100 chorych
      z <- sum(x[zdrowi])
      c <- sum(x[chorzy])
      p <- phyper(c-1,c+z,200-(c+z),100,FALSE)
      return(c(z,c,p))
}

do.call("rbind",lapply(7:106,function(x){f(dane[[x]])}))->wynik
wynik <- data.frame(wynik,sapply(wynik[,3],function(x){ifelse(x<0.1,"!","")}))
row.names(wynik)<- names(dane)[7:106]
names(wynik) <- c("zdrowi","chorzy","p-value","")
wynik
#mutacje ktore wplywaja na powstanie choroby
#przy poziomie istotnosci 0.05
row.names(wynik)[wynik[[3]]<0.05]->mutacje
mutacje

#oszacowanie FDR z wykorzystaniem pakietu qvalue
#przy progu istotnosci 0.05
qobj <- qvalue(wynik[,3])
max(qobj$qvalues[qobj$pvalues<=0.05])->fdr
fdr

m <- function(x){
  #zwraca wektor severity dla mutacji x
  #x powinno nalezec do przedzialu [7,106]
  mutacje <- x == 1
  dane$severity[chorzy & mutacje]
}

#srednie severity wsrod wszystkich chorych
mean(dane$severity[chorzy])->srednia

do.call("rbind",lapply(7:106,function(x){
  m(dane[[x]])->severity
  t.test(severity, mu = srednia, alternative="greater")$p.value -> pvalue
  return(c(length(severity),mean(severity),pvalue))
}))->wynik2
wynik2 <- data.frame(wynik2,sapply(wynik2[,3],function(x){ifelse(x<0.05,"!!","")}))
row.names(wynik2)<- names(dane)[7:106]
names(wynik2) <- c("n","srednia","p-value","")
wynik2
#mutacje ktore wplywaja pozytywnie na stopien dolegliwosci
#przy poziomie istotnosci 0.05
#średnia severity jest znaczaco wyzsza
row.names(wynik2)[wynik2[,3]<0.05]
wynik2[,3][wynik2[,3]<0.05]

#oszacowanie FDR z wykorzystaniem pakietu qvalue
#przy progu istotnosci 0.05
qobj <- qvalue(wynik2[,3])
max(qobj$qvalues[qobj$pvalues<=0.05])->fdr
fdr

#model liniowy zaleznosci choroby od ilosciowych skladnikow
model <- lm(severity ~ height + weight + sys_press + dias_press, data=dane[,2:6])
cor(dane$sys_press,dane$dias_press)
plot(dane$sys_press,dane$dias_press)
cor(dane$weight,dane$height)
plot(dane$weight,dane$height)
#step - regresja krokowa, startegia backward
model <- step(model,direction="backward")
summary(model)
#qqnorm(model$residuals)
#qqline(model$residuals)
#hist(model$residuals)
