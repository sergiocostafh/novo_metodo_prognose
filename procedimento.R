# PACOTES
pkgs <- c('ggplot2','minpack.lm','dplyr','cowplot','nlme')
lapply(pkgs,require,character.only=TRUE)

# IMPORTAR DADOS
dados <- read.csv2('dados.csv',header = T)

# PARÂMETROS
nclasses_sitio <- 3

# PAREAMENTO
for(i in 1:nrow(dados)){
  if(i==1){dados$SEQMED[i] <- 1}
  else{if(dados$COD_PAR[i]==dados$COD_PAR[i-1]){dados$SEQMED[i] <- dados$SEQMED[i-1]+1}
    else{dados$SEQMED[i] <- 1}}
}
maxseq <- max(dados$SEQMED)
for(i in 1:(maxseq-1)){
  b1 <- dados%>%filter(SEQMED==i)
  b1 <- b1[,c(1,7:19)]
  for(j in 1:(maxseq-1)){
    b2 <- dados%>%filter(SEQMED==j+1)
    b2 <- b2[,1:19]
    colnames(b1)[-1] <- paste(colnames(b1)[-1],'0',sep='')
    if(i==1 & j==1){dados_par <- inner_join(b2,b1,"COD_PAR")}
    else{b <- bind_rows(dados_par,inner_join(b2,b1))}
  }
}
##########################
# CLASSIFICAÇÃO DE SÍTIO #
##########################
hdom_model <- nlsLM(ALTURADOMINANTE~ALTURADOMINANTE0*((1-exp(-b*(IDADE)))/
                (1-exp(-b*(IDADE0))))^c,
              dados_par,
              start=list(b=.5,c=1.5),control=nls.control(maxiter=1000))

hdom_a <- nlsLM(ALTURADOMINANTE~hdom_a*((1-exp(-coef(hdom_model)[1]*(IDADE)))/
                  (1-exp(-coef(hdom_model)[1]*(Inf))))^coef(hdom_model)[2],
                dados_par,
                  start=list(hdom_a=30),control=nls.control(maxiter=1000))

hdompredict <- function(dados){hd <- coef(hdom_a)[1]*(1-exp(-coef(hdom_model)[1]*dados$IDADE))^coef(hdom_model)[2];return(hd)}

dados_par <- dados_par%>%mutate(FCI = round(ALTURADOMINANTE0/hdompredict(data.frame(IDADE=dados_par$IDADE0)),3),
                        SITIO = round(FCI*hdompredict(data.frame(IDADE=15)),1))
dados <- dados%>%mutate(FCI = round(ALTURADOMINANTE/hdompredict(dados),3),
               SITIO = round(FCI*hdompredict(data.frame(IDADE=15)),1))
range <- round(quantile(dados$SITIO,probs = c(.01,.99)))
limites <- seq(range[1],range[2],diff(range)/nclasses_sitio)
sitios <- c();for(i in 1:nclasses_sitio){sitios[i] <- round((limites[i]+limites[i+1])/2,1)};limites[1] <- -Inf;limites[length(limites)] <- Inf
dados$CLASSE_SITIO <- cut(dados$SITIO,limites,ordered_result = TRUE,labels = sitios)
limites <- seq(range[1],range[2],diff(range)/nclasses_sitio) # RETORNA AO VALOR ORIGINAL DOS LIMITES

curvas <- data.frame(IDADE=rep(seq(6,20,.1),nclasses_sitio*2+1))
curvas <- curvas%>%mutate(SITIO = rep(c(sitios,limites),each=length(seq(6,20,.1))),
                FCI = SITIO/hdompredict(data.frame(IDADE=15)),
                ALTURADOMINANTE = hdompredict(curvas)*FCI,
                TIPO = rep(c(rep('MEDIAS',nclasses_sitio),rep('LIMITES',nclasses_sitio+1)),each=length(seq(6,20,.1))))

ggplot(dados)+
  geom_line(aes(x=IDADE,y=ALTURADOMINANTE,group=COD_PAR),alpha=.25)+
  geom_line(aes(x=IDADE,y=ALTURADOMINANTE,group=SITIO),curvas%>%filter(TIPO=='MEDIAS'),linetype='dashed')+
  geom_line(aes(x=IDADE,y=ALTURADOMINANTE,group=SITIO),curvas%>%filter(TIPO=='LIMITES'))+
  scale_x_continuous(limits=c(0,25),expand=c(0,0))+
  scale_y_continuous(limits=c(0,35),expand=c(0,0),breaks=seq(0,30,5))+
  xlab('Idade (anos)')+
  ylab('Hdom (m)')+
  theme_cowplot()
#############################
# MODELAGEM DE VOLUME TOTAL #
#############################
vol_model <- nlsLM(VOLUMETOTAL~VOLUMETOTAL0*
                     (exp(-(b+b1*FCI)*(1/(IDADE^c))))/
                     (exp(-(b+b1*FCI)*(1/(IDADE0^c)))),
              dados_par,
              start=list(b=.5,b1=1,c=1.5),control=nls.control(maxiter=1000))
erro <- round(summary(vol_model)$sigma/mean(dados_par$VOLUMETOTAL)*100,1)

vol_a <- nlsLM(VOLUMETOTAL~vol_a*
                 (exp(-(coef(vol_model)[1]+coef(vol_model)[2]*FCI)*(1/(IDADE^coef(vol_model)[3]))))/
                 (exp(-(coef(vol_model)[1]+coef(vol_model)[2]*FCI)*(1/(Inf^coef(vol_model)[3])))),
               dados_par,
               start=list(vol_a=1000),control=nls.control(maxiter=1000))

# vol_a <- nlsLM(VOLUMETOTAL~(vol_a0+vol_a1*FCI)*(exp(-(coef(vol_model)[1]+coef(vol_model)[2]*FCI)*(1/(IDADE^coef(vol_model)[3]))))/
#                     (exp(-(coef(vol_model)[1]+coef(vol_model)[2]*FCI)*(1/(Inf^coef(vol_model)[3])))),
#                   dados_par,
#                   start=list(vol_a0=50,vol_a1=100),control=nls.control(maxiter=1000))

volpredict <- function(dados){vol <- coef(vol_a)[1]*(exp(-(coef(vol_model)[1]+coef(vol_model)[2]*dados$FCI)*(1/(dados$IDADE^coef(vol_model)[3]))));return(vol)}
curvas$VOLUMETOTAL <- volpredict(curvas)

ggplot(dados)+
  geom_line(aes(x=IDADE,y=VOLUMETOTAL,group=COD_PAR),alpha=.25)+
  geom_line(aes(x=IDADE,y=VOLUMETOTAL,group=SITIO),curvas%>%filter(TIPO=='MEDIAS'),linetype='dashed')+
  geom_line(aes(x=IDADE,y=VOLUMETOTAL,group=SITIO),curvas%>%filter(TIPO=='LIMITES'))+
  scale_x_continuous(limits=c(0,25),expand=c(0,0))+
  scale_y_continuous(limits=c(0,900),expand=c(0,0),breaks=seq(0,900,100))+
  xlab('Idade (anos)')+
  ylab('Volume Total (m³)')+
  theme_cowplot()
##############################
# MODELAGEM DE SOBREVIVÊNCIA #
##############################
n_model <- nlsLM(DENSIDADE~(DENSIDADE0^(1-b)-a*((1-b)/(1+c))*(ALTURADOMINANTE^(1+c)-ALTURADOMINANTE0^(1+c)))^(1/(1-b)),
            data=dados_par,
            start=list(a=0,b=0.1,c=.2),
            control=nls.control(maxiter=1000))
erro <- round(summary(n_model)$sigma/mean(dados_par$DENSIDADE)*100,1)

n0 <- nlsLM(DENSIDADE~(n0^(1-coef(n_model)[2])-coef(n_model)[1]*((1-coef(n_model)[2])/(1+coef(n_model)[3]))*(ALTURADOMINANTE^(1+coef(n_model)[3])-0.3^(1+coef(n_model)[3])))^(1/(1-coef(n_model)[2])),
               data=dados_par,
               start=list(n0=1000),
               control=nls.control(maxiter=1000))
npredict <- function(dados,n0){(n0^(1-coef(n_model)[2])-coef(n_model)[1]*((1-coef(n_model)[2])/(1+coef(n_model)[3]))*(dados$ALTURADOMINANTE^(1+coef(n_model)[3])-0.3^(1+coef(n_model)[3])))^(1/(1-coef(n_model)[2]))}

curvas$DENSIDADE <- npredict(curvas,coef(n0))
hdom_sim <- seq(10,30,.1)

ggplot()+
  geom_line(aes(x=ALTURADOMINANTE,y=DENSIDADE,group=COD_PAR),dados,alpha=.25)+
  geom_line(aes(x=hdom_sim,y=npredict(data.frame(ALTURADOMINANTE=hdom_sim),rep(unname(coef(n0)),length(hdom_sim)))),size=1.5,color='darkred',alpha=.5)+
  geom_line(aes(x=hdom_sim,y=npredict(data.frame(ALTURADOMINANTE=hdom_sim),rep(unname(coef(n0)*1.1),length(hdom_sim)))),linetype='longdash')+
  geom_line(aes(x=hdom_sim,y=npredict(data.frame(ALTURADOMINANTE=hdom_sim),rep(unname(coef(n0)*.9),length(hdom_sim)))),linetype='longdash')+
  geom_line(aes(x=hdom_sim,y=npredict(data.frame(ALTURADOMINANTE=hdom_sim),rep(unname(coef(n0)*1.2),length(hdom_sim)))),linetype='dashed')+
  geom_line(aes(x=hdom_sim,y=npredict(data.frame(ALTURADOMINANTE=hdom_sim),rep(unname(coef(n0)*.8),length(hdom_sim)))),linetype='dashed')+
  geom_line(aes(x=hdom_sim,y=npredict(data.frame(ALTURADOMINANTE=hdom_sim),rep(unname(coef(n0)*1.3),length(hdom_sim)))),linetype='dotted')+
  geom_line(aes(x=hdom_sim,y=npredict(data.frame(ALTURADOMINANTE=hdom_sim),rep(unname(coef(n0)*.7),length(hdom_sim)))),linetype='dotted')+
  scale_x_continuous(limits=c(0,35),expand=c(0,0))+
  scale_y_continuous(limits=c(0,2000),expand=c(0,0),breaks=seq(0,2000,200))+
  xlab('Hdom (m)')+
  ylab('Densidade (árvores/ha)')+
  theme_cowplot()
##################################
# MODELAGEM DA RELAÇÃO VCOM/VTOT #
##################################
dados <- dados%>%mutate(NV = DENSIDADE/VOLUMETOTAL)
vcom_function <- function(a,b,c,x){
  y <- 1-(a*(1-exp(-b*x^c)))
  ifelse(y<0,0,ifelse(y>1,1,y))}

vcom_model <- nls(VOLUMECOMERCIAL/VOLUMETOTAL~vcom_function(a,b,c,NV),dados,
              start=list(a=.8,b=.02,c=1),control=nls.control(maxiter = 1000))

erro <- round(summary(vcom_model)$sigma/mean(dados$VOLUMECOMERCIAL/dados$VOLUMETOTAL)*100,1)

ggplot()+
  geom_line(aes(NV,VOLUMECOMERCIAL/VOLUMETOTAL,group=COD_PAR),dados,alpha=.5)+
  geom_line(aes(x=seq(min(dados$NV),max(dados$NV)),y=predict(vcom_model,data.frame(NV=seq(min(dados$NV),max(dados$NV))))),color='darkred',size=1.5,alpha=.5)+
  xlab('N/VT (N. fustes/m³)')+
  ylab('Vol.Com/Vol.Tot.')+
  scale_y_continuous(limits=c(.8,1),expand = c(0,0))+
  scale_x_continuous(limits=c(0,10),expand = c(0,0))+
  theme_cowplot()
##########################
# MODELAGEM DO VOLUME S1 #
##########################
dados$S1P <- dados$S1/dados$VOLUMECOMERCIAL

s1_model <- nls(S1P~c/(b)^c*NV^(c-1)*exp(-(NV/(b))^c)*a,
              dados,
              start=list(a=1.148,b=4.337,c=1.721),
              control=nls.lm.control(maxiter = 1000))

erro <- round(summary(s1_model)$sigma/mean(dados$S1P)*100,1)

ggplot()+
  geom_line(aes(NV,S1P,group=COD_PAR),dados,alpha=.5)+
  geom_line(aes(x=seq(min(dados$NV),max(dados$NV),.1),y=predict(s1_model,newdata = data.frame(NV=seq(min(dados$NV),max(dados$NV),.1)))),size=2,alpha=.5,color='darkred')+
  xlab('N/VT (N. fustes/m³/ha)')+
  ylab('Vol.S1/Vol.Com.')+
  scale_color_viridis_c()+
  scale_y_continuous(limits=c(0,.3),expand = c(0,0))+
  scale_x_continuous(limits=c(0,10),expand = c(0,0))+
  theme_cowplot()
##########################
# MODELAGEM DO VOLUME S1 #
##########################
dados$S2P <- dados$S2/dados$VOLUMECOMERCIAL

reexp_alt <- function(a,b,x){
  ifelse(a+b^x>0, a+b^x,0)}

s2_model <- nlsLM(S2P~reexp_alt(a,b,NV),data=dados,
              start=list(a=-0.03,b=.6),control=nls.control(maxiter=1000))

erro <- round(summary(s2_model)$sigma/mean(dados$S2P)*100,1)

ggplot()+
  geom_line(aes(NV,S2P,group=COD_PAR),dados,alpha=.5)+
  geom_line(aes(x=seq(min(dados$NV),max(dados$NV)),y=predict(ms2,data.frame(NV=seq(min(dados$NV),max(dados$NV))))),
            color='darkred',size=1.5,alpha=.5)+
  xlab('N/VT (N. fustes/m³/ha)')+
  ylab('Vol.S2/Vol.Com.')+
  scale_y_continuous(limits=c(0,.6),expand = c(0,0))+
  scale_x_continuous(limits=c(0,10),expand = c(0,0))+
  theme_cowplot()
#########################
# VOLUME POR SORTIMENTO #
#########################
curvas <- curvas%>%
  mutate(NV = DENSIDADE/VOLUMETOTAL)
curvas <- curvas%>%
  mutate(VOLUMECOMERCIAL = predict(vcom_model,curvas)*VOLUMETOTAL,
  S1P = predict(s1_model,curvas),
  S2P = predict(s2_model,curvas),
  PRP = 1-(S1P+S2P))
