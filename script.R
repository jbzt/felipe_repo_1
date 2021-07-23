# cargo los paquetees que voy a usar
library(ggpubr)
library(dplyr)
library(ggplot2)

###################################
# preparar la data para trabajar  #
###################################

# cargar la data (esto lo hice solo la primera vez, despues cree otro archivo que es mas rapido de importar)
# covid_felipe_total = read.delim("data/covid_felipe_total.txt", header = TRUE, sep = "\t")
#saveRDS(object = covid_felipe_total, file = 'data_procesada')


# cargo la data de forma rapida
covid_felipe_total = readRDS('data/data_procesada')

# creo una copia de la data por si la quiero editar
tabla = covid_felipe_total

#emprolijo un poco la data
tabla$positivity = as.factor(tabla$positivity)
tabla$ct = as.numeric(as.character(tabla$ct))
tabla$ct

# creo una lista con los genes que quiero estudiar
genes.interes = c('IFNGR1' ,"IFNGR2" , "IFNG", "JAK1", "JAK2", "NFKB1", "MX1", "STAT6", "STAT1", "JUN", "CBL", "RAPGEF1", "RUNX3", "CEBPB", "MAP2K6", "PRKACA")

# para saber si un gen esta en la tabla
"NFKB1" %in% colnames(tabla)
"GEN_QUE_NO_ESTA" %in% colnames(tabla)

# como saber si todos los genes que quiero estudiar estan en la tabla

genes.interes %in% colnames(tabla) # me da una respuesta por cada gen



###################################################################
#  ejemplos graficos de expresion genica en base a la condicion   #
###################################################################

# para graficar un gen con puntos
ploteo = ggplot(data = tabla, aes(x=positivity, y=IFNGR2 )) +
  geom_jitter(aes(shape=positivity, color=positivity), size=3)+
  xlab("condición") +
  ylab(paste("IFNGR2 expression \n log2 (norm counts +1)")) +
  theme(legend.position = "bottom") +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title =element_text(size = 25),
        legend.position = 'none') +
  stat_compare_means()
ploteo



# para graficar un gen con puntos y agregue barra de promedio
ploteo = ggplot(data = tabla, aes(x=positivity, y=log2(CD9+1))) +
  geom_jitter(aes(shape=positivity, color=positivity), size=3,width = 0.2)+
  xlab("condición") +
  ylab(paste("CD9 expression \n log2 (norm counts +1)")) +
  theme(legend.position = "bottom") +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title =element_text(size = 25),
        legend.position = 'none') +
  stat_compare_means() + 
  stat_summary(fun=mean,
               geom="point",
               shape= '_',
               size=14)
ploteo



# para hacer un violin plot (las partes que tienen # estan "dsactivadas", si se lo sacas, van a aparecer)
ploteo = ggplot(data = tabla, aes(x=positivity, y=log2(CD9+1))) +
  geom_violin(trim=T, aes(color=positivity,fill=positivity)) +
  #Para pintar lineas y rellenar adentro el de arriba)
  geom_boxplot(geom ="errorbar", width=0.25, outlier.shape = NA) +
  geom_jitter(width=0.10, size=0.85) +
  #geom_jitter(aes(shape=positivity, color=positivity), size=3,width = 0.2)+
  xlab("condición") +
  ylab(paste("CD9 expression \n log2 (norm counts +1)")) +
  #Para agregar media en el box
  stat_summary(fun="mean", geom="point", shape=22, size=4) +
  theme(legend.position = "bottom") +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title =element_text(size = 25),
        legend.position = 'none') +
  stat_compare_means() 
# stat_summary(fun=mean,
#              geom="point",
#              shape= '_',
#              size=14)
ploteo


###############################################
#   graficos finales de expresion en loop     #
###############################################

#Creamos una carpeta para meter los graficos
dir.create("graficos_expresion")

# para hacer un loop que grafico y guarde todos los graficos al mismo tiempo (el 16 indica la cantidad de gnes que queiro ver)
for (i in 1:16){
  gen = genes.interes[i]
  ploteo = ggplot(data = tabla, aes(x=positivity, y=log2(tabla[,gen]+1))) +
    geom_violin(trim=FALSE, aes(color=positivity,fill=positivity)) + 
    ylab(paste(gen, "expression \n log2 (norm counts +1)")) +
    xlab("Condición") +
    geom_jitter(width=0.10, size=0.85, alpha=0.4) +
    stat_summary(fun="mean", geom="point", shape=22, size=4) +
    theme(legend.position = "bottom") +
    theme_bw() +
    theme(axis.text = element_text(size = 15),
          axis.title = element_text(size = 15),
          plot.title =element_text(size = 25),
          legend.position = 'none') +
    stat_compare_means(label.x = 1.23) 
  
  #Guardando ando
  ggsave(filename = paste(gen,'ploteo_expresion.png', sep = '_'),
         plot = ploteo,
         dpi = 600,
         device="png",
         path = 'graficos_expresion')
  
  print(ploteo)
}


#########################################################
# ejemplos graficos de correlacion de los genes con ct  #
#########################################################

# grafico de correlacion basico
ploteo.cor = ggplot(data = tabla, aes(x = tabla$ct, y = tabla$STAT1 )) + geom_point()
ploteo.cor

# grafico de correlacion basico pero con la data de expresion logaritmica
ploteo.cor = ggplot(data = tabla, aes(x = tabla$ct, y = log2(tabla$STAT1 + 1) )) + geom_point()
ploteo.cor


# poner bien los nombres de los ejes
# agregar que en el grafico aparezca estadistica (pearson o spearman)
# hacer que los puntos tengan un color si vienen de covid y otro si vienen de sano





###############################################
#   graficos finales de correlacion en loop     #
###############################################


#Creamos una carpeta para meter los graficos
dir.create("graficos_correlacion")
