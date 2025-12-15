####### funcion reflectancia 2.0.0 ##########
###funcion automatizacion reflectancia 
########################
## VERIFICACIÓN E INSTALACIÓN AUTOMÁTICA DE LIBRERÍAS
#########################

paquetes <- c(
  "dplyr",
  "readr",
  "tidyverse",
  "data.table",
  "glue",
  "ggplot2",
  "readxl",
  "scales",
  "zoo",
  "stats",    # stats viene con R base, no necesita instalación
  "mgcv"      # mgcv también viene con R base, pero lo incluimos por completitud
)

# Instalar los que faltan (excluyendo los de base que no necesitan instalación)
paquetes_a_instalar <- paquetes[!paquetes %in% c("stats", "mgcv")]  # estos suelen venir preinstalados
paquetes_faltantes <- paquetes_a_instalar[!paquetes_a_instalar %in% installed.packages()[,"Package"]]

if (length(paquetes_faltantes) > 0) {
  cat("Instalando paquetes faltantes:", paste(paquetes_faltantes, collapse = ", "), "\n")
  install.packages(paquetes_faltantes, dependencies = TRUE)
} else {
  cat("Todos los paquetes necesarios ya están instalados.\n")
}

# Cargar todas las librerías
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyverse)
  library(data.table)
  library(glue)
  library(ggplot2)
  library(readxl)
  library(scales)
  library(zoo)
  library(stats)
  library(mgcv)
})

cat("Librerías cargadas correctamente.\n\n")

##funcion 1: ordena los datos del txt
reflec_order<-function(txt){
  dato<- txt[-c(1:15),] ###quitar filas que no nos sirven
  dato<- as.data.frame(dato) ## pasar a dataframe
  dato<-separate(dato, col = 1, into = c("banda","reflectancia"), sep = "\t") ##separar en columnas de banda y reflectancia
  return(dato)
}


##funcion 2: aplica la funcion de arriba sobre todos los txt en una carpeta 
## y los deja en una lista

upload_reflec<- function(folder_direc){
  folder_direc ###dirección de la carpeta
  archivos <- list.files(path = folder_direc, pattern = "\\.txt", full.names = TRUE) #direccion de los archivos
  nombres_archivos <- tools::file_path_sans_ext(basename(archivos)) ##los nombres de cada uno
  carpeta_arch <- lapply(seq_along(archivos), function(i) {
    r <- read.csv(archivos[i])  ## cargar los archivos
    return(r)
  })
  carpeta_arch<- lapply(carpeta_arch, reflec_order) ##ordenar los archivos en un df usando la funcion 1
  names(carpeta_arch)<-nombres_archivos ##le damos el nombre de cada archivo txt
  return(carpeta_arch)
}


##funcion 3: ordena las reflectancias
lista_df<-function(lista, muestreo = 1){
  df_reflec<-bind_rows(lista, .id = "origen") ##pasar todos los df de una lista a uno solo
  df_reflec<-pivot_wider(df_reflec,names_from = banda, values_from = reflectancia) ##las bandas quedan como columnas
  df_reflec$origen <- gsub("Reflection_", "", df_reflec$origen) ##dejamos solo el numero como dato de profundidad
  df_reflec$origen<-as.numeric(df_reflec$origen) ## profundidad se pasa a numerico
  colnames(df_reflec)[1]<-"depth" ##se cambia origen por depth
  ##reemplazamos la profundidad por los valores del intervalo de muestreo con su unidad de medida
  df_reflec$depth <- seq(1,nrow(df_reflec),by=1)*muestreo
  return(df_reflec)
}

#####funcion 4 100% automatizada integrando todas las anteriores
##### guarda una matriz conlos datos de reflectancia por ancho de banda
manejar_reflectancia<-function(folder, muestreo = 0.2, name = "reflectancia.csv"){
  archivos_lista<-upload_reflec(folder) ##funcion 2
  df<-lista_df(archivos_lista, muestreo = muestreo) ##funcion 3
  write.table(df,name, sep = ";", row.names = FALSE) ##se guardan los datos como un excel
  return(df)
}

##### calculo de indices automatizado

##funcion 5 para ordenar los datos para calculo de reflectancias
calculus_order <- function(df){
  df <- df %>%      ###para pasar la matriz a un solo df con 3 columnas (depth, band, reflectance)
    pivot_longer(
      cols = -depth, ##no transformar depth
      names_to = "band", ##las columnas a un vector de nombre band
      values_to = "reflectance" ##los valores de reflectancia a un vector de nombre reflectance
    )
  df$band<-as.numeric(df$band) ##pasamos a numerico
  df$reflectance<-as.numeric(df$reflectance) ##pasamos a numerico
  return(df)
}

########## para calcular el RABD ################

#### RADB660-670 is calculated based on the formula presented in: 
###  authors: Bert Rein ? Frank Sirocko 
### title: "In-situ reflectance spectroscopy - analysing techniques for high-resolution pigment logging in sediment cores"
### reference: Int J Earth Sci (Geol Rundsch) (2003) 92:143 (DOI 10.1007/s00531-002-0308-5)
#### formula: RABD660-670=((x*R590+y*R730)/x+y)/Rmin(660-670)

### with: 

### RABD660-670 = Relative absorption band depth at band 600-670 nm 
### R590 = Reflectance at 590 nm
### R730 = Reflectance at 730 nm
### Rmin(660-670) = reflectance minimum of the 660-670 nm range
### x = Number of spectral bands between Rmin(660-670) and Reflectance at 730 nm
### y = Number of spectral bands between Rmin(660-670) and Reflectance at 590 nm

### for more information also see: Butz, C., Grosjean, M., Fischer, D., Wunderle, S., Tylmann, W. and Rein, B.
### 2015. Hyperspectral imaging spectroscopy: a promising method for the biogeochemical analysis of lake sediments. 
### Journal of Applied Remote Sensing, 9(1), pp.096031-096031.
##funcion 6
calculate_rabd<-function(df){
  df<-calculus_order(df) ##function 5
  n660to670<-subset(df, band>660 & band<670) ##separarel rango de bandas entre 660 y 670
  min660_670 <- n660to670 %>%  #####buscar la reflectancia minima entre los rango de 660-670
    mutate(minimo = pmin(n660to670$reflectance)) %>%
    group_by(depth) %>% 
    summarise(min_reflec = min(minimo, na.rm = TRUE))
  y <- df %>% ###buscar el numero de bandas espectrales entre 590 a 660
    filter(band > 590.091 & band < 660.019) %>%
    summarise(cantidad = n_distinct(band))
  
  x <- df %>% ###buscar el numero de bandas espectrales entre 670 a 730
    filter(band > 670.136 & band < 730.146) %>%
    summarise(cantidad = n_distinct(band))
  
  reflec_590<-df %>% ##reflectancia en la banda 590
    subset(band==590.091) %>%
    group_by(depth)
  
  reflec_730<-df %>% ##reflectancia en la banda 730
    subset(band==730.146) %>%
    group_by(depth)
  
  ##calculo del RABD
  RABD_660_670<-((((x$cantidad*reflec_590$reflectance) + 
                     (y$cantidad*reflec_730$reflectance))/(x$cantidad+y$cantidad))/min660_670$min_reflec)
  ##agregamos profundidad
  RABD_660_670<-data.frame(RABD_660_670=RABD_660_670,depth = min660_670$depth)
  return(RABD_660_670)
}

#### calcular ratio 660/670 ######

ratio_660_670<-function(df){
  df<-calculus_order(df)
  reflec_660<-subset(df,df$band==660.019)
  reflec_670<-subset(df,df$band==670.136)
  ratio <-reflec_660$reflectance/reflec_670$reflectance
  return(ratio)
}

#### calcular ratio 590/640 ####

ratio_590_640<-function(df){
  df<-calculus_order(df)
  reflec_590<-subset(df,df$band==590.091)
  reflec_640<-subset(df,df$band==640.309)
  ratio <-reflec_590$reflectance/reflec_640$reflectance
  return(ratio)
}

###### calcular ratio 570/630 ########

ratio_570_630<-function(df){
  df<-calculus_order(df)
  reflec_570<-subset(df,df$band==570.113)
  reflec_630<-subset(df,df$band==630.235)
  ratio <-reflec_570$reflectance/reflec_630$reflectance
  return(ratio)
}


#### calcula todos los indices anteriores de manera automatizada #####
indices_reflec<-function(df){
  rabd<-calculate_rabd(df)
  rat_660_670<-ratio_660_670(df)
  rat_590_640<-ratio_590_640(df)
  rat_570_630<- ratio_570_630(df)
  datos = data.frame(rabd,ratio_660_670 =rat_660_670,
                     ratio_590_640 = rat_590_640,
                     ratio_570_630 = rat_570_630)
  return(datos)
}


##### funcion automatica que maneja multiples carpetas de reflectancia a la vez #####

auto_reflectancia <- function(carpeta_principal, muestreo = 0.2) {
  
  # Obtener la fecha actual en formato YYYY-MM-DD
  fecha_hoy <- format(Sys.Date(), "%Y-%m-%d")
  
  # Definir nombres de las carpetas de salida con la fecha
  carpeta_reflectancias <- paste0("reflectancias_", fecha_hoy)
  carpeta_indices       <- paste0("indices_", fecha_hoy)
  
  # Crear las carpetas de salida si no existen
  dir.create(carpeta_reflectancias, showWarnings = FALSE, recursive = TRUE)
  dir.create(carpeta_indices,       showWarnings = FALSE, recursive = TRUE)
  
  # Listar las subcarpetas directas en la carpeta principal
  subcarpetas <- list.dirs(carpeta_principal, full.names = TRUE, recursive = FALSE)
  
  # Filtrar solo directorios
  subcarpetas <- subcarpetas[sapply(subcarpetas, function(d) file.info(d)$isdir)]
  
  total <- length(subcarpetas)
  
  if (total == 0) {
    cat("No se encontraron subcarpetas en la carpeta principal.\n")
    return(invisible(NULL))
  }
  
  # Barra de progreso
  pb <- txtProgressBar(min = 0, max = total, style = 3, width = 60, char = "=")
  cat("Iniciando procesamiento de", total, "subcarpetas...\n")
  cat("Resultados se guardarán en:\n")
  cat("  Reflectancias →", carpeta_reflectancias, "\n")
  cat("  Índices       →", carpeta_indices, "\n\n")
  
  # Procesar cada subcarpeta
  for (i in seq_along(subcarpetas)) {
    subcarpeta <- subcarpetas[i]
    nombre_subcarpeta <- basename(subcarpeta)
    
    # Ruta para guardar el archivo de reflectancia
    archivo_reflectancia <- file.path(carpeta_reflectancias, 
                                      paste0("reflectancia_", nombre_subcarpeta, ".csv"))
    
    # Aplicar manejar_reflectancia (guarda automáticamente el CSV y devuelve el df)
    df_reflectancia <- manejar_reflectancia(subcarpeta, muestreo = muestreo, name = archivo_reflectancia)
    
    # Calcular los índices
    df_indices <- indices_reflec(df_reflectancia)
    
    # Ruta para guardar el archivo de índices
    archivo_indices <- file.path(carpeta_indices, 
                                 paste0("indices_", nombre_subcarpeta, ".csv"))
    
    # Guardar el dataframe de índices como CSV con separador ";"
    write.table(df_indices, archivo_indices, sep = ";", row.names = FALSE, col.names = TRUE)
    
    # Actualizar barra de progreso
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  
  cat("\n¡Procesamiento completado!\n")
  cat("Archivos de reflectancia guardados en:", carpeta_reflectancias, "\n")
  cat("Archivos de índices guardados en:      ", carpeta_indices, "\n")
}



