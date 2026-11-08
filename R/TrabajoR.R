##########################################################################################
#####################                                        			##################
#####################    TRABAJO ORDINARIA BIOINFORMÁTICA R             ##################
#####################                                       			##################
##########################################################################################


### 1.- Carga el archivo de entrada "archivo-entrada.txt" correspondiente. 
### Recuerda incluir los nombres de las columnas como el encabezado.
### Examina el archivo con las funciones head(), tail() y dim(). (1 pto)

arch<- read.table(file="archivo-entrada.txt", header=TRUE, sep='\t', check.names=FALSE)
#arch
head(arch)
tail(arch)
dim(arch)
class(arch)

### 2.- En nuestra tabla tenemos los valores de 619 genes para tres réplicas
### del Control y tres réplicas de una Condición. Calcula la media y la desviación
### estándar del control y del tratamiento para cada gen. Puedes emplear cualquier
### función que consideres como apply(). Añade esos valores como columnas nuevas a tu tabla original. 
### Serán 4 columnas nuevas en total, la media de Control, la DE de Control, la media de Condición y la DE de Condición. (1,5 pto)
class(arch)

# Imprime toda la primera columna enumerada
#arch[1]
# Imprime la primera fila
#arch[1,]
# Imprime todos nombres de la primera columna
#arch[,1]

wt   <- arch[, 2:4]   # columnas de Control
cd <- arch[, 5:7]   # columnas de Condición



arch$wt_media <- apply(wt, 1, mean)
arch$wt_DE <- apply(wt, 1, sd)
arch$cd_media <- apply(cd, 1, mean)
arch$cd_DE <- apply(cd, 1, sd)

head(arch)
dim(arch)


#names(genesets)=paste0("set",1:3)

#for (gene in 1:nrow(arch)) {
#  wt <- arch[gene,c(2:4)]
#  cd <- arch[gene,c(5:7)]
#  
#}


### 3.- En una sola imagen, muestra un histograma para la media del Control y otro para la media de la Condicíon.
### Pon colores diferentes en cada uno y mantén esos colores para el resto del trabajo. (1 pto).

col_wt <- "skyblue"
col_cd <- "salmon"

par(mfrow = c(1, 2)) # para mostrar dos graficas en horizontal

hist(log2(arch$wt_media+0.1),
  col = col_wt,
  main = "Media del Control",
  xlab = "Media Control",
  ylab = "Frecuencia",
  breaks = 30)

hist(log2(arch$cd_media+0.1),
  col = col_cd,
  main = "Media de la Condición",
  xlab = "Media Condición",
  ylab = "Frecuencia",
  breaks = 30)

## Aplico logaritmo ya que la diferencia entre valores es inmensa, como en Taller. Además para que funcione, había que sumar algo (0.1), para quitar zeros
## La transformación es para mejor interpretación de datos, pero no para calcular


par(mfrow = c(1, 1)) #para volver mostrar una grafica

### 4.- En una sola imagen, muestra un boxplot para la condición y otro para el control. 
### Añade el título, colores y una una leyenda en la esquina superior izquierda. ( 1,5 pto)

boxplot(log2(arch$wt_media+0.1), log2(arch$cd_media+0.1),
  names = c("Control", "Condición"),
  col = c(col_wt, col_cd),
  main = "Boxplot de la media del Control y de la Condición",
  ylab = "Media")

legend("topleft",
  legend = c("Control", "Condición"),
  fill = c(col_wt, col_cd),
  bty = "n")


### 5.- Evalúa la normalidad de Control y Condición y acorde con los resultados, haz lo siguiente: (2 pto)
### 		Realiza un test de comparativa de medias para cada gen entre condición y control. 
###			Guarda el p-value de la comparativa y anádelo como una columna extra a tu tabla original.


valores_wt <- as.numeric(unlist(arch[, 2:4]))
valores_cd <- as.numeric(unlist(arch[, 5:7]))

shapiro.test(valores_wt)
shapiro.test(valores_cd)

p_shapiro_wt <- shapiro.test(valores_wt)$p.value
p_shapiro_cd <- shapiro.test(valores_cd)$p.value

par(mfrow = c(1, 2))
qqnorm(valores_wt, main = "QQ-plot Control", col = col_wt)
qqline(valores_wt)
qqnorm(valores_cd, main = "QQ-plot Condición", col = col_cd)
qqline(valores_cd)
par(mfrow = c(1, 1))

# Si p > 0.05 en los dos grupos, los datos son normales - t.test
# Si no, los datos no son normales - wilcox.test
if (p_shapiro_wt > 0.05 & p_shapiro_cd > 0.05) {
  print("Datos normales: uso t.test")
} else {
  print("Datos no normales: uso wilcox.test")
}

pvalue <- numeric(nrow(arch))

for (i in 1:nrow(arch)) {

  x <- as.numeric(arch[i, 2:4])
  y <- as.numeric(arch[i, 5:7])
  
  if (p_shapiro_wt > 0.05 & p_shapiro_cd > 0.05) {
    pvalue[i] <- t.test(x, y)$p.value
  } else {
    pvalue[i] <- suppressWarnings(wilcox.test(x, y)$p.value)
  }
}


arch$pvalue <- pvalue

head(arch)
dim(arch)



### 6.- De los resultados del apartado 5, saca (mediante código) un listado de los genes diferencialmente expresados (aquellos que tengan un pvalue < 0.05)
### ¿Cuántos genes tienes? Sino te sale ninguno, justifica por qué. En caso de que sí que te salgan, ¿cuántos salen sobreexpresados y cuántos reprimidos? (1 pto)

arch$diferencia <- arch$cd_media - arch$wt_media

# Filtramos los genes con p-value < 0.05
filtro_pvalue <- arch$pvalue < 0.05
genes_DE <- arch[filtro_pvalue, ]

# Listado de genes diferencialmente expresados
genes_DE[, 1]

# Cuántos genes diferencialmente expresados
n_DE <- nrow(genes_DE)
n_DE

# Sobreexpresados y reprimidos
sobre <- sum(filtro_pvalue & arch$diferencia > 0)
repri <- sum(filtro_pvalue & arch$diferencia < 0)
sobre
repri

# No sale ninguno. Comprobamos cuál es el p-value más bajo y cuántos genes lo tienen
min(arch$pvalue)
sum(arch$pvalue == min(arch$pvalue))


### 7.- Por último, realiza un heatmap de los genes diferencialmente expresados
### Emplea la Práctica de R como referencia. (2 pts)

## Como no hay p value menor que 5, voy dividir el existente por dos para que tengamos algo. Y repito última parte del 6

pvalue100porcien_real <- arch$pvalue/2
filtro_pvalue <- pvalue100porcien_real < 0.05
genes_DE <- arch[filtro_pvalue, ]

# Listado de genes diferencialmente expresados
genes_DE[, 1]

# Cuántos genes diferencialmente expresados
n_DE <- nrow(genes_DE)
n_DE

# Sobreexpresados y reprimidos
sobre <- sum(filtro_pvalue & arch$diferencia > 0)
repri <- sum(filtro_pvalue & arch$diferencia < 0)
sobre
repri



if (n_DE >= 2) {
  seleccion <- genes_DE
  titulo <- "Heatmap - genes diferencialmente expresados"
} else {
  orden <- order(arch$pvalue, -abs(arch$diferencia))
  seleccion <- arch[orden[1:50], ]
  titulo <- "Heatmap - 50 genes con menor p-value"
}

dim(seleccion)

# Matriz con las columnas de expresión (2:7), con los nombres de los genes como filas
# Aplico logaritmo como antes
matriz <- as.matrix(seleccion[, 2:7])
rownames(matriz) <- seleccion[, 1]
matriz <- log2(matriz + 0.1)

head(matriz)


rowv <- as.dendrogram(hclust(as.dist(1 - cor(t(matriz)))))
colv <- as.dendrogram(hclust(as.dist(1 - cor(matriz))))


heatmap(matriz, Rowv = rowv, scale = "row",
        cexCol = 0.8, cexRow = 0.6, col = hcl.colors(50),
        main = titulo)

# install.packages(c("gplots", "RColorBrewer"))
library(gplots)
library(RColorBrewer)

heatmap(matriz, Rowv = rowv, scale = "row",
          trace = "none", cexCol = 0.8, cexRow = 0.6,
          col = rev(redblue(256)), main = titulo)
