#############################################################################
#
# PRACTICA 1
#
# Expresión diferencial de genes de ratón
# Microarray de Affymetrix (Affymetrix Murine Genome U74A version 2 MG_U74Av2
# Origen de los datos: GEO GSE5583 (http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5583)
# Publicación: Mol Cell Biol 2006 Nov;26(21):7913-28.  16940178 (http://www.ncbi.nlm.nih.gov/pubmed/16940178)
#
# Muestras: 3 Wild Type x 3 Histone deacetylase 1 (HDAC1)
#
# R código original (credits): Ahmed Moustafa
#
#
##############################################################################

# Instalar RCurl

if (!requiereNamespace("BiocManager"))
install.packages ("BiocManager")
BiocManager::install ("Rcurl")


install.packages ("RCurl")

# Si esto falla, que seguro lo hace tratar de instalarlo usando el menú, Paquetes, Servidor Spain A Coruña, RCurl

# Cargamos el paquete y los datos

library("RCurl")

"url"=getURL("http://bit.ly/GSE5583_data",followlocation=TRUE)

data=as.matrix(read.table(text =url, row.names = 1,header= T ))

# Chequeamos las dimensiones de los datos, y vemos las primeras y las últimas filas

head(data)
tail(data)
dim(data)
# Hacemos un primer histograma para explorar los datos
hist(data, col="pink", main="GSE5583 -Histogram")

# Transformamos los datos con un logaritmo se pueden cuando representar de imágenes.
# ¿Qué pasa si hacemos una transformación logarítima de los datos? ¿Para qué sirve?
#Se utilizan para representar los datos de manera más atractiva,ayuda a normalizar los datos, hacer que las varianzas sean más constantes, 
#se utilizan para las pruebas de contraste de hipótesis. 
#Para representación de imagenes, Pero luego a la hora de análisis de datos hay que hacerlos con los datos originales sin transformar. 

data2=log2(data)
hist(data2,col="pink", main="GSE5583-Histogram (log2)")


# Hacemos un boxplot con los datos transformados. ¿Qué significan los parámetros que hemos empleado?
# ¿Qué es un boxplot?
#Para conseguir los cuartiles los datos
boxplot(data2, col=c("blue", "blue", "blue", "yellow", "yellow", "yellow"),main="GSE5583 - Boxplots",las=2")



# Hacemos un hierarchical clustering de las muestras basándonos en un coeficiente de correlación
# de los valores de expresión. ¿Es correcta la separación?
#Si es correcta la separación separa las wild type de las mutantes. 

hc = hclust(as.dist(1-cor(data2)))
plot(hc,main="GSE5583 - Hierarchical Clustering")

#######################################
# Análisis de Expresión Diferencial 
#######################################

# Primero separamos las dos condiciones. ¿Qué tipo de datos has generado?
#Matrices que separa, las células wild type de las mutantes.

head(data)
wt<-data[,1:3]
ko<- data[,4:6]

head(wt)
head(ko)

class(wt)
class(ko)
# Calcula las medias de las muestras para cada condición. Usa apply
#Se trabaja con las medias de las 3 réplicas tendrían qeu ser equivalentes con valores muy carcanos

wt.mean<- apply(wt,1,mean)
ko.mean<- apply(ko,1,mean)

head(wt.mean)
head(ko.mean)

# ¿Cuál es la media más alta?
#La media más baja no es 0 es 1,9

limit=max(wt.mean,ko.mean)
limit
limit_repress= min(wt.mean,ko.mean)
limit_repress

# Ahora hacemos un scatter plot (gráfico de dispersión)

plot(ko.mean ~ wt.mean, xlab = "WT", ylab = "KO", main= "GSE5583 - ScatterPlot", col="yellow", xlim=c (0,limit), ylim=c (0,limit))
 
# Añadir una línea diagonal con abline

abline(0,1,col="red")

abline (1,2, col="orange")
abline (h=2, col="green")
abline (v=3, col="violet")


# ¿Eres capaz de añadirle un grid?
grid()

# Calculamos la diferencia entre las medias de las condiciones

diff.mean <- wt.mean - ko.mean
head(diff.mean)

# Hacemos un histograma de las diferencias de medias

hist(diff.mean,col="pink")

# Calculamos la significancia estadística con un t-test.
# Primero crea una lista vacía para guardar los p-values
# Segundo crea una lista vacía para guardar las estadísticas del test.
# OJO que aquí usamos los datos SIN TRANSFORMAR. ¿Por qué?
# ¿Cuántas valores tiene cada condición?

prueba1= wt[1,]
prueba2= ko[1,]
t.test (prueba1, prueba2)
head(wt)
pval=t.test(prueba1, prueba2)$p.value
pval

pvalue=NULL
tstat=NULL

nrow(data)
dim(data)

for (i in 1:nrow (data)){
x = wt[i,]
y = ko[i,]
t=t.test(x,y)
pvalue[i] = t$p.value
tstat[i] = t$statistic
}
head(pvalue)	


# Ahora comprobamos que hemos hecho TODOS los cálculos

length(pvalue)

# Hacemos un histograma de los p-values.
# ¿Qué pasa si le ponemos con una transformación de -log10?

hist(pvalue)
hist(-log10(pvalue))

# Hacemos un volcano plot. Aquí podemos meter la diferencia de medias y la significancia estadística

plot(diff.mean,-log10(pvalue), main= "GSE5583 - Volcano")
plot(diff.mean,pvalue, main= "GSE5583 - Volcano")

# Queremos establecer que el mínimo para considerar una diferencia significativa, es con una diferencia de 2 y un p-value de 0.01
# ¿Puedes representarlo en el gráfico?

diff.mean_cutoff=2
pvalue_cutoff=0.01

abline(v=diff.mean_cutoff, col="yellow", lwd=3)
abline(h=-log10(pvalue_cutoff), col="orange", lwd=3)


# Ahora buscamos los genes que satisfagan estos criterios
# Primero hacemos el filtro para la diferencia de medias (fold)

filter_by_diff.mean=abs(diff.mean)>- diff.mean_cutoff
dim(data[filter_by_diff.mean,])

# Ahora el filtro de p-value

filter_by_pvalue = pvalue <= pvalue_cutoff
dim(data[filter_by_pvalue, ])

# Ahora las combinamos. ¿Cuántos genes cumplen los dos criterios?

filtered_combined=filter_by_diff.mean $ filter_by_pvalue
filtered= data([filtered_combined,]

dim(filtered)
head(filtered)


# Ahora generamos otro volcano plot con los genes seleccionados marcados en rojo


plot(diff.mean, -log10(pvalue), main="GSE5583 - Volcano filtered")
points(diff.mean[filtered_combined], -log10(pvalue[filtered_combined]),col="red")



# Ahora vamos a marcar los que estarían sobreexpresados (rojo) y reprimidos (azul). ¿Por qué parece que están al revés?
plot(diff.mean, -log10(pvalue),main="GSE5583 - Volcano (filtered)")
points(diff.mean[filtered_combined & diff.mean < 0],-log10(pvalue[filtered_combined & diff.mean < 0]),col="red")
points(diff.mean[filtered_combined & diff.mean > 0],-log10(pvalue[filtered_combined & diff.mean > 0]),col="blue")


# Ahora vamos a generar un mapa. Para ello primero tenemos que hacer un cluster de las columnas y los genes 

heatmap(filtered)

rowv = as.dendrogram(hclust(as.dist(1-cor(t(filtered)))))
colv=as.dendrogram(hclust(as.dist(1-cor(filtered))))
heatmap(filtered, Rowv=rowv, Colv=colv, cexCol=0.7, labRow=FALSE)


# ¿Qué es cada parámetro que hemos usado dentro de la función heatmap?
# ¿Eres capaz de cambiar los colores del heatmap? Pista: usar el argumento col y hcl.colors


# Ahora vamos a crear un heatmap más chulo. Para ello necesitamos dos paquetes: gplots y RcolorBrewer

install.packages("gplots")
install.packages("RcolorBrever")

library(gplots)
library(RColorBrewer)


# Hacemos nuestro heatmap
	
heatmap.2(filtered, Rowv=rowv, Colv=colv, cexCol=0.7,col = rev(redblue(256)), scale = "row")


# Guardamos los genes diferencialmente expresados y filtrados en un fichero

#Rowv=determines if and how the row dendrogram should be reordered. 
#Colv=determines if and how the column dendrogram should be reordered. 
#scale=indicates if the values should be centered and scaled in either the row direction or the column direction

# Lo guardamos en un archivo PDF

pdf ("GSE5583_DE_Heatmap_Carlota_Montesinos_Guareño_pdf")
heatmap.2 (filtered, Rowv=rowv, Colv=colv, cexCol=0.7,
	col = rev(redblue(256)), scale = "row",labRow=FALSE)
dev.off()

# Guardamos los genes diferencialmente expresados y filtrados en un fichero

write.table (filtered, "GSE5583_DE_Carlota_Montesinos_Guareño_.txt", sep = "\t",quote = FALSE)
	
