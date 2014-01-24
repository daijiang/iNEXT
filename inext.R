install.packages('devtools')
library(devtools)
install_github('iNEXT','JohnsonHsieh')
library(iNEXT)
data(ant)
y1 <- iNEXT(ant$h50m, q=0, datatype="incidence")
y2 <- iNEXT(ant$h50m, q=1, datatype="incidence")
y3 <- iNEXT(ant$h50m, q=2, datatype="incidence")
plot(y1, ylim=c(1, 300), main="Ant at h50m")
lines(y2, col=2)
lines(y3, col=3)
#==========================================================================
library("plyr")
library(reshape2)
library(stringr)

data=read.csv("allsites.csv", stringsAsFactors=F)
str(data)
length(unique(data$Species)) # 512 species
length(unique(data$Site)) # 244 sites
length(unique(data$SampleEvent))
data$site=paste(data$Site, data$SampleEvent,  sep="_") # 262 unique things...

spcode=read.csv("spcode.csv", stringsAsFactors=F)
str(spcode)
spcode$Code=as.character(spcode$Code)

# change code in species from numner to latin names....
for (i in 1: dim(spcode)[1]){
  data$Species[data$Species==spcode$Code[i]]=spcode$Species[i] 
}
unique(data$Species)
write.csv(data, "VegAllSites.NS.csv", row.names=F)

# change data to wide table
data.wide=dcast(data=data[,c(6,4)], site~Species)
row.names(data.wide)=data.wide[,1]
data.wide=data.wide[,-1]
max(data.wide[,]) # 345?
sort(colwise(max)(data.wide))

# how many quadrats per site
x=scan("clipboard")
x1=scan("clipboard")
x2=scan("clipboard")
quadsite=data.frame(site=paste(x1, x, sep="-"), quadnumber=x2)

row.names(data.wide) %in% quadsite$site
a=match(quadsite$site, row.names(data.wide))
data.wide$quadnumber[a[!is.na(a)]]=  quadsite$quadnumber[!is.na(a)]

data.frame(arrange(data.frame(site=quadsite$site[!is.na(a)] , quad=quadsite$quadnumber[!is.na(a)]), site),
q  =data.wide["quadnumber"], si=row.names(data.wide)) # check whether quad number is matched correctly...

data.wide=data.wide[,c(513,1:512)] # the data is ready...
#---------------------------------------------------------------------------------------------------------

# test1$summary  T   U S.obs S.hat  C.hat
# $order
# $interpolation t       qD qD.95.LCL qD.95.UCL        SC SC.95.LCL SC.95.UCL
# $extrapolation

# site.summary=vector("list")
site.detail=vector("list",length=dim(data.wide)[1])
names(site.detail)=row.names(data.wide)
site.detail1=site.detail

for  (i in 1: dim(data.wide)[1]){
  b=t(data.wide[i,])
  test0=iNEXT(b, q=0, datatype="incidence", nboot=100, endpoint=403,knots=60) # maximum quadrat number
  test1=iNEXT(b, q=1, datatype="incidence", nboot=100, endpoint=403,knots=60)
  test2=iNEXT(b, q=2, datatype="incidence", nboot=100, endpoint=403,knots=60)
  test0$interpolation$q=0;test0$extrapolation$q=0
  test1$interpolation$q=1;test1$extrapolation$q=1
  test2$interpolation$q=2;test2$extrapolation$q=2
  site.detail[[i]]=rbind(test0$interpolation, test1$interpolation, test2$interpolation)
  site.detail1[[i]]=rbind(test0$extrapolation, test1$extrapolation, test2$extrapolation)
  site.summary[[i]]=data.frame(site=row.names(data.wide)[i], test0$summary)
}


site.summary=ldply(site.summary)
names(site.summary)[2]="quadrat.size"
names(site.summary)[3]="total.incidence"
names(site.summary)[5]="Chao2"
names(site.summary)[6]="coverage"
site.summary=arrange(site.summary, coverage) # max is 1!
write.csv(site.summary, "siteSummary.csv", row.names=F)

siteDetail=ldply(site.detail)
hh=site.detail1[-which(names(site.detail1)=="1041-46")]
siteDetail1=ldply(hh)
siteDetail2=rbind(siteDetail, siteDetail1)
siteDetail2=arrange(siteDetail2, .id, q,t)
names(siteDetail2)=c("site", "quadrat.size","Hill", "Hill.95.lower", "Hill.95.upper",
                     "coverage","coverage.95.lower", "coverage.95.upper", "q.Hill")
siteDetail2=siteDetail2[,c(1,2,9,3:8)]
siteDetail2$site=str_replace_all(siteDetail2$site, "-", "_")
write.csv(siteDetail2, "rarefy_extrap.csv", row.names=F)

# sp rich ~ log(area), lm intercept, slope...
sprich=subset(siteDetail2, q.Hill=="0")
regsemilog=function(df) {
  model=lm(df$Hill~log(df$quadrat.size))
  coeff=coef(summary(model))
  regression=data.frame(site=unique(df$site), intercept=coeff[1], slope=coeff[2],
                        inter.se=coeff[3],  slope.se=coeff[4], inter.p=coeff[7],
                        slope.p=coeff[8])
  return(regression)
}
semilog.lm=ddply(sprich, .(site), regsemilog)
write.csv(semilog.lm, "semilog_lm.csv", row.names=F)


# diversity of all sites at 0.95 coverage.
max(siteDetail2$coverage)
tt=subset(siteDetail2, quadrat.size==403)
tt=arrange(tt, coverage)
min(tt$coverage) # the min coverage when quadrat size = 403 is 0.9513865

c95=ddply(siteDetail2, .(site), summarize, nn=min(which(coverage>=0.95)))
ttt=merge(c95, siteDetail2, "site")
tttt=ddply(ttt, .(site), summarize, nn=unique(nn))
c95tt=vector("list")
for(i in 1: dim(tttt)[1]){
  aa=ttt[ttt$site==tttt$site[i],]
  bb=aa[c(unique(aa$nn), unique(aa$nn)+60, unique(aa$nn)+120),]
  c95tt[[i]]=bb
}
c95ttt=ldply(c95tt)
spdiveristy.c95=dcast(c95ttt[,c(1,4,5)], site~q.Hill)
names(spdiveristy.c95)=c("site_sampleEvent","sp-rich","shannon", "simpson")
# spdiveristy.c95$site_sampleEvent=str_replace_all(spdiveristy.c95$site_sampleEvent, "-", "_")
write.csv(spdiveristy.c95, "spDiv_95coverage.csv", row.names=F)


