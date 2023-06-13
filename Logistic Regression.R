###=1. PACKAGES=###
library(car)
library(foreign)
library(shapefiles)
library(nlme)
library(maptools)
library(grid)
library(gstat)
library(latticeExtra)
library(splines)
library(spdep)
library(classInt)
library(RColorBrewer)
library(CARBayes)
library(CARBayesST)
library(ctv)
library(sp)
library(rgdal)
library(Matrix)
library(DCluster)
library(CAMAN)
library(rgeos)
library(raster)
library(spatialreg)
library(tmap)
library(gplots)
library(INLA)
library(MoranST)

###=2. INPUT DATA=###
data=read.csv(file.choose(), sep=",", header=TRUE)
data

attach(data)
summary(data)

###=3. SET MAPPING=###
#Mapping#Melihat Peta Jabat
jabar.map=readOGR(dsn="D:/Kampus/Semester 5/Epidem/JABAR.shp")
plot(jabar.map,main="PETA PROVINSI JAWA BARAT",col=c("#FF9966","#996666","#CC6666","#FFFF99","#99CCFF","#66CC99","#9933FF","#339999","#33FF66","#993399","#FF6699","#990000","#CC0033","#FFCCCC","#CCFF66","#00FFCC","#9999CC","#CC33CC","#CCFF66","#FFFF00","#003333","#333399","#336600","#00CC66","#009933","#CC9933","#FF00CC"))
text(jabar.map,jabar.map$KABKOT,cex=0.5)
no=c(1:27)
Coordinate=coordinates(jabar.map)

plot(jabar.map,main="PETA PROVINSI JAWA BARAT",col="white",axes=TRUE)
text(Coordinate-0.001,label=no,cex=1)

getinfo.shape("D:/Kampus/Semester 5/Epidem/JABAR.shp") 
coords<-coordinates(jabar.map)

#Set Data Mapping
northarrow<-list("SpatialPolygonsRescale",layout.north.arrow(),
                 offset<-c(106,108),scale=0.1,which=2)
scalebar<-list("SpatialPolygonRescale",layout.scale.bar(),
               offset<-c(-7,-6),scale=0.1,
               fill<-c("transparent","black"),which=2)
text1<-list("sp.text",c(106,108),"0",which=2)
text2<-list("sp.text",c(-7,-6),"0.1 KM",which=2)

###=4. TAKSIRAN SMR=###
SMR19=data$Y[1:27]/data$Ei[1:27]
SMR20=data$Y[28:54]/data$Ei[28:54]
boxplot(SMR19,SMR20,main="Boxplot SMR",names=c("SMR19","SMR20"))
SMR=c(SMR19,SMR20)
SMR1=data.frame(SMR19,SMR20)
#Deskriptif Risiko Relatif
summary(SMR)

###=5. MATRIKS PEMBOBOTAN SPASIAL=###

jabar.queen=poly2nb(jabar.map,queen=TRUE)
WQ=nb2listw(jabar.queen) 
Wqueen=as(as_dgRMatrix_listw(WQ), "CsparseMatrix") 
Wqueen=as.matrix(Wqueen) 
Wqueen 

###=6. MENGUJI AUTOKORELASI SPASIAL=###
y=data$Y
X1=data$X1
X2=data$X2
X3=data$X3
X4=data$X4
Ei=data$Ei
formula=y~X1+X2+X3+X4
model=glm(formula=y~X1+X2+X3+X4,family="poisson")
summary(model)

spasial=data[,4]
Y=as.matrix(spasial)
W=Wqueen 
MoranSTMCResult<-MoranST.MC(Y,W,nsim=100)
MoranSTMCResult

#moran residual
residual=resid(model)
Moran19=moran.test(residual[1:27],listw=WQ)
Moran20=moran.test(residual[28:54],listw=WQ)

Moran=c(Moran19$estimate[1],Moran20$estimate[1])
pv.Moran=c(Moran19$p.value,Moran20$p.value)
Moran.Indeks=data.frame(Moran,pv.Moran)
Moran.Indeks

###=7. PEMILIHAN MODEL=##
#Menaksir Model#
Model1=y~X1+X2+X3+X4
ModelRun1=inla(Model1,family="poisson",data=data,E=Ei,control.predictor=list(compute=TRUE,link=1),control.compute=list(dic=TRUE,cpo=TRUE))
summary(ModelRun1)

IDs0=data$IDs0
Model2=y~X1+X2+X3+X4+f(IDs0,model="bym",graph=Wqueen)
ModelRun2=inla(Model2,family="poisson",data=data,E=Ei,control.predictor=list(compute=TRUE,link=1),control.compute=list(dic=TRUE,cpo=TRUE))
summary(ModelRun2)

IDt0=data$IDt0
Model3=y~X1+X2+X3+X4+f(IDt0,model="rw2")
ModelRun3=inla(Model3,family="poisson",data=data,E=Ei,control.predictor=list(compute=TRUE,link=1),control.compute=list(dic=TRUE,cpo=TRUE))
summary(ModelRun3)

Model4=y~X1+X2+X3+X4+f(IDs0,model="bym",graph=Wqueen)+f(IDt0,model="rw2")
ModelRun4=inla(Model4,family="poisson",data=data,E=Ei,control.predictor=list(compute=TRUE,link=1),control.compute=list(dic=TRUE,cpo=TRUE))
summary(ModelRun4)

data$IDt1=data$IDt0
Model5=y~X1+X2+X3+X4+f(IDs0,model="bym",graph=Wqueen)+f(IDt0,model="rw2")+f(IDt1,model="rw2",group=IDs0,control.group=list(model="iid"))
ModelRun5=inla(Model5,family="poisson",data=data,E=Ei,control.predictor=list(compute=TRUE,link=1),control.compute=list(dic=TRUE,cpo=TRUE))
summary(ModelRun5)

#Menghitung Ukuran Kecocokan Model#
Y1=ModelRun1$summary.fitted.values$mean*Ei
R.sqr1=sum((Y1-mean(y))^2)/sum((y-mean(y))^2)
Y2=ModelRun2$summary.fitted.values$mean*Ei
R.sqr2=sum((Y2-mean(y))^2)/sum((y-mean(y))^2)
Y3=ModelRun3$summary.fitted.values$mean*Ei
R.sqr3=sum((Y3-mean(y))^2)/sum((y-mean(y))^2)
Y4=ModelRun4$summary.fitted.values$mean*Ei
R.sqr4=sum((Y4-mean(y))^2)/sum((y-mean(y))^2)
Y5=ModelRun5$summary.fitted.values$mean*Ei
R.sqr5=sum((Y5-mean(y))^2)/sum((y-mean(y))^2)

DIC=c(ModelRun1$dic$dic,ModelRun2$dic$dic,ModelRun3$dic$dic,ModelRun4$dic$dic,ModelRun5$dic$dic)
R.sqr=c(R.sqr1,R.sqr2,R.sqr3,R.sqr4,R.sqr5)
data.frame(DIC,R.sqr)

###=8. SUMMARY MODEL TERBAIK=###
#Model Terbaik#
Model_5=y~X1+X2+X3+X4+f(IDs0,model="bym",graph=Wqueen)+f(IDt0,model="rw2")+f(IDt1,model="rw2",group=IDs0,control.group=list(model="iid"))
ModelRun_5=inla(Model_5,family="poisson",data=data,E=Ei,control.predictor=list(compute=TRUE,link=1),control.compute=list(dic=TRUE,cpo=TRUE))
summary(ModelRun_1)


###=9. MENGHITUNG RR MENGGUNAKAN SPATIO-TEMPORAL MODEL=###
#Model Fix#
Modelfix=y~X1+X2+X3+X4+f(IDs0,model="bym",graph=Wqueen)+f(IDt0,model="rw2")+f(IDt1,model="rw2",group=IDs0,control.group=list(model="iid"))
ModelRun=inla(Modelfix,family="poisson",data=data,E=Ei,control.predictor=list(compute=TRUE,link=1),control.compute=list(dic=TRUE,cpo=TRUE))
#Menghitung Nilai Risiko Relatif Spatio-Temporal CAR Model#
ST19=ModelRun$summary.fitted.values$mean[1:27]
ST20=ModelRun$summary.fitted.values$mean[28:54]

boxplot(ST19,ST20, main="Boxplot Spatio-Temporal Model",names=c("ST19","ST20"))
ST=c(ST19,ST20)
ST1=data.frame(ST19,ST20)
#Deskriptif Risiko Relatif
summary(ST1)

#Membuat Boxplot Risiko Relatif Spatio-Temporal Model#
boxplot(SMR,ST,main="Boxplot SMR dan Spatio Temporal Model",names=c("SMR","ST"))
#Detail
var(ST19)
var(SMR19)
var(ST20)
var(SMR20)
var(SMR)
var(ST)
summary(ST)
summary(SMR)

#Membuat Boxplot Perbandingan Penaksir SMR dan Spatio-Temporal Model#
par(mfrow=c(1,1))
boxplot(SMR19,ST19,main="Boxplot SMR dan ST 2019",names=c("SMR19","ST19"))
boxplot(SMR20,ST20,main="Boxplot SMR dan ST 2020",names=c("SMR20","ST20"))


###=10. MEMBUAT BOXPLOT SE PENAKSIR SMR DAN SPATIO-TEMPORAL MODEL=###
Var.SMR19=SMR19/Ei[1:27]
SE.SMR19=sqrt(Var.SMR19)
Var.SMR20=SMR20/Ei[28:54]
SE.SMR20=sqrt(Var.SMR20)


#Menghitung Standar Error Penaksir Spatio-Temporal Model#
SE.ST19<-ModelRun$summary.fitted.values[,5][1:27]
SE.ST20<-ModelRun$summary.fitted.values[,9][28:54]

#Membuat Boxplot Perbandingan Standar Error Penaksir SMR dan Spatio-Temporal Model#
par(mfrow=c(1,1))
boxplot(SE.SMR19,SE.ST19,main="Boxplot SE SMR dan ST 2019",names=c("SE.SMR19","SE.ST19"))
boxplot(SE.SMR20,SE.ST20,main="Boxplot SE SMR dan ST 2020",names=c("SE.SMR20","SE.ST20"))


###=11. MAPPING SEBARAN DBD=###
no=c(1:27)
dataRR=data.frame(no,ST19,ST20)

datamap <- jabar.map[,-c(1:17)]
data.combined <- cbind(datamap, dataRR)

RR=c(ST19[1:27],ST20[1:27])
plotvar=RR
nlcr=5
plotclr=brewer.pal(nlcr,"Reds")
class=classIntervals(plotvar,nlcr,style="quantile")
colcode=findColours(class,plotclr)
plot<- spplot(data.combined,c("ST19","ST20"),names.attr=c("ST19","ST20"),colorkey=list(space="bottom"),
              main="Taksiran Risiko Relatif Spatio-Temporal Model",col.regions=plotclr,at=round(class$brks,digits=3),
              sp.layout=list(northarrow),as.table=TRUE)
label<-layer(sp.text(coordinates(data.combined),txt=c(1:27),
                     txt.cex=0.1, cex=0.7, pos=0.8))
plot+label