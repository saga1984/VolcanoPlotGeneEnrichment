#' @title Consenso de DEGs (limma/DESeq2/edgeR), conteos y volcanos con etiquetas
#'
#' @description
#' Lee resultados de **limma**, **DESeq2** y **edgeR**, armoniza columnas
#' (log2FC, Pval, Padj), clasifica genes (Up/Down/Not Significant) con un umbral
#' de \code{|log2FC| >= corte} y \code{Padj < 0.05}, construye un **consenso 2 de 3**,
#' genera tablas de universo y significativos, y produce volcanos (con y sin
#' etiquetas de genes de interés) guardando salidas en disco.
#'
#' @details
#' - Los archivos de entrada deben contener columnas
#' - El "consenso 2 de 3, 2oo3" marca un gen como común cuando, coinciden 2 de las 3 herramientas usadas.
#' - Genera archivos en formato CSV y guarda.
#'
#' @param archivo_limma Ruta al archivo de resultados de \strong{limma} (CSV/TSV; si es TSV usar \code{sep = "\t"}).
#' @param archivo_DESeq2 Ruta al archivo de resultados de \strong{DESeq2} (CSV/TSV; si es TSV usar \code{sep = "\t"}).
#' @param archivo_edgeR Ruta al archivo de resultados de \strong{edgeR} (CSV/TSV; si es TSV usar \code{sep = "\t"}).
#' @param vector_yfg Vector o data.frame con \code{Gene_ID} de genes a etiquetar en el volcan (Your Favorite Genes).
#' @param trinotate_file Ruta al archivo Trinotate para mapear \code{Gene_ID}, símbolo (Blast).
#' @param tratamiento Etiqueta de tratamiento/condición; se usa en nombres de archivos, carpetas y títulos.
#' @param salidas Carpeta base donde se crearán subdirectorios y salidas.
#' @param colores Vector de 3 colores para Down/Not/Up en el volcán. Default: \code{c("#00AFBB", "gray", "#bb0c00")}.
#' @param formatos Formato(s) de imagen a exportar (ej., \code{"jpeg"}, \code{"png"}). Default: \code{"jpeg"}.
#' @param resolucion Resolución (dpi) para las figuras exportadas. Default: \code{300}.
#' @param titulo_grafica (Opcional) Título manual para la gráfica (no utilizado actualmente).
#' @param titulo_legenda (Opcional) Título manual para la leyenda (no utilizado actualmente).
#' @param corte Umbral absoluto de \code{log2FC} para clasificar Up/Down. Default: \code{1}.
#'
#' @return `NULL`. Escribe tablas y figuras en disco y deja un objeto
#'   \code{volcano_label_<tratamiento>} en el \code{.GlobalEnv}.
#'
#' @examples
#' \dontrun{
#' degs_por_condicion(
#'   archivo_limma   = "limma.tsv",
#'   archivo_DESeq2  = "deseq2.tsv",
#'   archivo_edgeR   = "edger.tsv",
#'   vector_yfg      = c("TRINITY_DN1_c0_g1","TRINITY_DN2_c0_g1"),
#'   trinotate_file  = "Trinotate.csv",
#'   tratamiento     = "NaCl_3h",
#'   salidas         = "results",
#'   formatos        = c("png","jpeg"),
#'   corte           = 1
#' )
#' }
#'
#' @seealso
#' \itemize{
#' \item \code{ggrepel::geom_label_repel}, \code{ggplot2::ggplot}
#' \item \code{dplyr::case_when}, \code{dplyr::count}
#' }
#'
#' @import data.table
#' @import tidyverse
#' @importFrom ggrepel geom_label_repel
#' @importFrom stringr str_count
#' @importFrom ggplot2 ggplot aes geom_point geom_vline geom_hline
#' @importFrom ggplot2 scale_color_manual labs scale_x_continuous ggtitle theme_classic theme element_rect
#'
#' @keywords RNA-seq DEGs volcano meta-analysis consensus limma DESeq2 edgeR
#' @export


# definir la funcion
degs_por_condicion <- function(archivo_limma,
                               archivo_DESeq2,
                               archivo_edgeR,
                               vector_yfg,
                               trinotate_file,
                               tratamiento,
                               salidas,
                               colores = c("#00AFBB", "gray", "#bb0c00"),
                               formatos = "jpeg",
                               resolucion = 300,
                               titulo_grafica = "",
                               titulo_legenda = "",
                               corte = 1){
  
  
  # cargar librerias
  library(tidyverse)
  library(data.table)
  library(gridExtra)
  library(ggrepel)
  
  # importar archivos a R
  archivos <- mget(ls()[grepl("archivo_", ls())])
  
  # list para guardar resultados
  datos <- vector(mode = "list", length = length(archivos))
  
  # nombrar data frames
  names(datos) <- names(archivos)
  
  # poblar lista de resultados  
  for(archivo in 1:length(datos)){
    
    # poblar lista d data frames
    datos[[archivo]] <- read.csv(archivos[[archivo]], sep = "\t", header = TRUE, na.strings = NA)
    
    # quitar columnas de datos crudos
    datos[[archivo]] <- datos[[archivo]][, !grepl("_", colnames(datos[[archivo]]))]
    
    # agregar columnas de Gene_ID
    datos[[archivo]]["Gene_ID"] <- rownames(datos[[archivo]])
    
    
    # if statment para cambiar nombres de columnas de cada data frame
    if(names(datos[archivo]) == "archivo_DESeq2"){
      
      # dejar solo columnas de interes
      datos[[archivo]] <- datos[[archivo]][, c("log2FoldChange", "pvalue", "padj", "Gene_ID")]
      # modificar nombres de columnas
      colnames(datos[[archivo]]) <- c("log2FC", "Pval", "Padj", "Gene_ID")
      
    } else if(names(datos[archivo]) == "archivo_edgeR"){
      
      # dejar solo columnas de interes
      datos[[archivo]] <- datos[[archivo]][, c("logFC", "PValue", "FDR", "Gene_ID")]
      # modificar nombres de columnas
      colnames(datos[[archivo]]) <- c("log2FC",  "Pval", "Padj", "Gene_ID") 
      
    } else if(names(datos[archivo]) == "archivo_limma"){
      
      # dejar solo columnas de interes
      datos[[archivo]] <- datos[[archivo]][, c("logFC", "P.Value", "adj.P.Val", "Gene_ID")]
      # modificar nombres de columnas
      colnames(datos[[archivo]]) <- c("log2FC", "Pval", "Padj", "Gene_ID") 
    }
    
    # categorizar por "Downregulated", "Not significant", "Upregulated"
    datos[[archivo]] <- datos[[archivo]] %>% mutate(
      DEGs = case_when(
        Padj > 0.05 & log2FC > -1 * corte & log2FC < corte ~ "Not Significant",
        Padj < 0.05 & log2FC > corte ~ "Upregulated",
        Padj < 0.05 & log2FC < -1 * corte ~ "Downregulated",
        TRUE ~ "Not Significant"
      )
    ) 
    
  }
  
  # Function to sort each data frame by row names
  sort_df_by_rownames <- function(df) {
    df[order(rownames(df)), , drop = FALSE]  # Sort by row names
  }
  
  # Sort each data frame in the list
  sorted_datos <- lapply(datos, sort_df_by_rownames)
  
  # Combine sorted data frames using cbind
  big_df_degs_cbind <- do.call(cbind, sorted_datos)
  
  ########################## obtener universo de datos ###########################
  
  # dividir un data frame en 3 data frames: 1 para pvalue, 1 para pvalue ajustado y 1 para log2fc
  df_universe_log2FC <- big_df_degs_cbind[, grepl("log2FC", colnames(big_df_degs_cbind))]
  df_universe_pval <- big_df_degs_cbind[, grepl("Pval", colnames(big_df_degs_cbind))]
  df_universe_padj <- big_df_degs_cbind[, grepl("Padj", colnames(big_df_degs_cbind))]
  
  # obtener promedio de valores
  df_universe_log2FC[, "log2FC_prom"] <- rowMeans(df_universe_log2FC, na.rm = TRUE)
  df_universe_pval[, "Pval_prom"] <- rowMeans(df_universe_pval, na.rm = TRUE)
  df_universe_padj[, "Padj_prom"] <- rowMeans(df_universe_padj, na.rm = TRUE)
  
  # agregar columna de Gene_id a cada data frame para hacer un merge
  df_universe_log2FC[, "Gene_ID"] <- rownames(df_universe_log2FC) 
  df_universe_pval[, "Gene_ID"] <- rownames(df_universe_pval)
  df_universe_padj[, "Gene_ID"] <- rownames(df_universe_padj)
  
  # hacer merge para obtener data frame final
  df_universe_clean_tmp <- merge(df_universe_log2FC[, c("log2FC_prom", "Gene_ID")], 
                                 df_universe_pval[, c("Pval_prom", "Gene_ID")], 
                                 by = "Gene_ID")
  
  # hacer merge para obtener data frame final
  df_universe_clean <- merge(df_universe_clean_tmp[, c("log2FC_prom", "Pval_prom", "Gene_ID")], 
                             df_universe_padj[, c("Padj_prom", "Gene_ID")], 
                             by = "Gene_ID")
  
  # modificar nombre de columnas
  colnames(df_universe_clean) <- c("Gene_ID", "log2FC", "Pval", "Padj")
  
  # volver a crear columna de DEGs
  df_universe_clean[, "DEGs"] <- case_when(
    df_universe_clean$Padj < 0.05 & df_universe_clean$log2FC > 1 * corte ~ "Upregulated",
    df_universe_clean$Padj < 0.05 & df_universe_clean$log2FC < -1 * corte ~ "Downregulated",
    df_universe_clean$Padj > 0.05 | 
      (df_universe_clean$log2FC < 1 * corte & df_universe_clean$log2FC > -1 * corte)~ "Not Significant")
  
  # QC should be 0
  # nrow(df_universe_clean[is.na(df_universe_clean$DEGs),])
  
  # crear variable de nuevo_dir
  nuevo_dir <- paste(salidas, "/", tratamiento, 
                     sep = "")
  
  # crear directorio de archivos de degs nuevos
  if(dir.exists(nuevo_dir)){
    print(paste("El directorio ", nuevo_dir, ", ya existe", sep = ""))
  } else {
    dir.create(nuevo_dir, recursive = TRUE)
  }
  
  ####### obtener tablas de strings (ids comunes a 2 de 3 herramientas) #######
  
  # filtrar gene clasificacion DEG
  df_degs <- big_df_degs_cbind[, grepl("DEGs", colnames(big_df_degs_cbind))]
  
  # dejar solo nombres de herramientas
  colnames(df_degs) <- gsub("\\..*$", "", colnames(df_degs))
  colnames(df_degs) <- gsub("^.*_", "", colnames(df_degs))
  
  # unique(df_degs[1,])
  
  # obtener lista de Gene_id comunes a 2 de 3 herramientas
  for(i in 1:nrow(df_degs)){
    # contar numero de Not Significants 
    df_degs[i, "Conteo_2d3"] <- sum(str_count(df_degs[i,], "Not Significant"))
  }
  
  # remover filas donde el conteo de no significativos es mayor que 1
  for(i in 1:nrow(df_degs)){
    # crear columna de status
    df_degs[i, "Common_2oo3"] <- ifelse(df_degs$Conteo_2d3[i] > 1, "","Common2oo3")
  }
  
  # remover columna de conteos dos de tres
  df_degs[,"Conteo_2d3"] <- NULL
  
  # agregar columna de Gene_IDs
  df_degs[, "Gene_ID"] <- row.names(df_degs) 
  
  # agregar conteo 
  df_universe_clean <-  merge(df_degs[, c("Gene_ID", "Common_2oo3")],
                              df_universe_clean,
                              by = "Gene_ID")
  
  # degs df
  degs_sig <- df_universe_clean[df_universe_clean$Common_2oo3 == "Common2oo3" &
                                  (df_universe_clean$DEGs == "Upregulated" |
                                     df_universe_clean$DEGs == "Downregulated"),]
  
  ############################# guardar archivos ##################################
  
  # guardar/crear archivo en windows
  write.csv(x = df_universe_clean,
            file = paste(nuevo_dir, "/", "degs_universe_", tratamiento, ".csv", sep = ""),
            row.names = FALSE)
  
  # guardar/crear archivo en windows
  write.csv(x = degs_sig,
            file = paste(nuevo_dir, "/", "degs_sig_", tratamiento, ".csv", sep = ""),
            row.names = FALSE)
  
  # change Treatment name
  tratamiento <-  case_when(
    tratamiento == "DESH" ~ "Dehydration",
    tratamiento == "REH" ~ "Rehydration",
    tratamiento == "ABA" ~ "Abscisic acid",
    tratamiento == "NACL" ~ "NaCl",
    tratamiento == "SORB" ~ "Sorbitol",
    TRUE ~ as.character(tratamiento))
  
  ########################## graficar volcano plot comunes ###########################
  
  # crear nuevo objeto con nomnbre corto
  file <- df_universe_clean
  
  # definir maximos y minimos en x
  maximo_x <- max(file$log2FC)
  minimo_x <- min(file$log2FC)
  
  # definir maximos en y
  maximo_y <- min(file$Padj)
  
  # contar up
  up_c <- count(file[file$DEGs == "Upregulated" & file$Common_2oo3 != "",])
  
  # contar down
  down_c <- count(file[file$DEGs == "Downregulated" & file$Common_2oo3 != "",])
  
  # contar up
  up_t <- count(file[file$DEGs == "Upregulated",])
  
  # contar down
  down_t <- count(file[file$DEGs == "Downregulated",])
  
  # contar_NS
  not_significant_c <- count(file[file$DEGs == "Not Significant" & file$Common_2oo3 != "",])
  not_significant_t <- count(file[file$DEGs == "Not Significant",])
  
  # crear grafica
  volcano <- ggplot(data = file, 
                    aes(x = log2FC, 
                        y = -log10(Padj), 
                        col = DEGs)) +
    geom_vline(xintercept = c(-1*corte, corte), 
               col = "gray", 
               linetype = 'dashed') +
    geom_hline(yintercept = -log10(0.05), 
               col = "gray", 
               linetype = 'dashed') +
    geom_point(size = 2) +
    scale_color_manual(values = colores,                                
                       labels = c(paste("Downregulated", " (", down_c, ")", sep = ""),
                                  "Not significant", 
                                  paste("Upregulated", " (", up_c, ")", sep = ""))) +
    labs(color = 'Diferentially Expressed Genes',
         x = expression("log"[2]*"FC"), y = expression("-log"[10]*"P-value")) +
    scale_x_continuous(breaks = seq(round(minimo_x) - 2, round(maximo_x) + 2, 2)) +
    ggtitle(paste(tratamiento , "- Treated Protonemal DEGs", sep = " ")) +
    theme_classic() +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1.2))
  
  # Create a new column "delabel" to de, that will contain the name of the top 30 differentially expressed genes (NA in case they are not)
  # df$delabel <- ifelse(df$gene_symbol %in% head(df[order(df$padj), "gene_symbol"], 30), df$gene_symbol, NA)
  
  # visualizar en R
  print(volcano)
  
  
  print(paste("Total: Upreglated = ", up_t))
  print(paste("Total: Downreglated = ", down_t))
  print(paste("Total: Not Significant = ", not_significant_t))
  print(paste("Common 2 out of 3: Upreglated = ", up_c))
  print(paste("Common 2 out of 3: Downreglated = ", down_c))
  print(paste("Total Not Significant = ", not_significant_c))
  
  # save output
  out11 <- "Total Upreglated"
  out12 <- up_t
  
  out21 <- "Total Downreglated"
  out22 <- down_t
  
  out31 <- "Total Not Significant"
  out32 <- not_significant_t
  
  out41 <- "Common 2 out of 3 Upreglated" 
  out42 <-  up_c
  
  out51 <- "Common 2 out of 3 Downreglated"
  out52 <- down_c
  
  out61 <- "Common 2 out of Not Significant" 
  out62 <- not_significant_c
  
  # create df
  df_count <- data.frame(matrix(nrow = 6, ncol = 2))
  
  # define colnames
  colnames(df_count) <- c("DEGs", "Count")
  
  # populate df
  df_count$DEGs <- c(out11, out21, out31, out41, out51, out61)
  df_count$Count<- c(out12, out22, out32, out42, out52, out62)
  df_count$Count <- as.numeric(unlist(df_count$Count))
  
  # guardar/crear archivo en windows
  write.csv(x = df_count,
            file = paste(nuevo_dir, "/", "degs_count_total_vs_common_", tratamiento, ".csv", sep = ""),
            row.names = FALSE)
  
  # crear subdirectorio de graficas de volcan
  dir.create(paste(nuevo_dir, "/Graficas_volcano", sep = ""))
  
  # guardar imagen en formatos pre-establecidos
  for(i in formatos) {
    
    # crear y guardar los heatmaps
    match.fun(i)(paste(nuevo_dir, 
                       "/Graficas_volcano/",
                       "Volcano_Plot_", 
                       paste(tratamiento, collapse = "_"), 
                       ".", i, sep = ""),
                 res = resolucion,
                 width = 3000,
                 height = 2000)
    
    # crear y guardar el heatmpat euclidean
    print(
      
      volcano
      
    )
    
    dev.off()
    
  }
  
  # agregar columna de genes, importar trinotate
  trinotate <- read.csv(trinotate_file)
  
  # agregar columna de genes, importar trinotate
  trinotate <- fread(trinotate_file, select = c(1,3,7))
  
  # modificar nombre de columnas
  colnames(trinotate) <- c("Gene_ID", "Blastx", "Blastp")
  
  # eliminar todo lo que va despues del _ incluyendo el mismo _
  trinotate$Blastx <- gsub("\\_.+$", "", trinotate$Blastx)
  trinotate$Blastp <- gsub("\\_.+$", "",trinotate$Blastp)
  
  # crear una nueva columna de blast (considerando tanto blastx como blastp)
  trinotate$Blast <- ifelse(trinotate$Blastp != ".", 
                            trinotate$Blastp, 
                            trinotate$Blastx)
  
  # nombre de columna de vector_yfg
  vector_yfg <- as.data.frame(vector_yfg)
  
  # modificar (asignar) nombre de columnas
  colnames(vector_yfg) <- "Gene_ID"
  
  # unir vector de Gene_id con trinotate
  unido <- merge(trinotate, vector_yfg, by = "Gene_ID")
  
  # filtrar solo columnas importantes
  unido <- unido[,c("Gene_ID", "Blast")]
  
  # sustituir "." por NA
  unido$Blast <- gsub("\\.", NA, unido$Blast)
  
  # remover filas con NA
  unido <- unido[!is.na(unido$Blast),]
  
  # transformar: contar cuantas veces aparece cada gen de blast 
  conteos_genes <- unido %>% group_by(Gene_ID, Blast) %>% dplyr::summarise(conteo = n()) 
  
  # # separar metodo 1: por Gene_id, dplyr
  # conteos_genes_separados <- conteos_genes %>% group_split(Gene_id)
  
  # separar metodo 2: vector de ids unicos
  id_unicos <- unique(conteos_genes$Gene_ID)
  
  # for loop
  for(i in 1:length(id_unicos)){
    assign(paste("split", id_unicos[i], sep = "_"), 
           subset(conteos_genes, Gene_ID == id_unicos[i]))
  }
  
  # crear lista 
  genes_separados <- vector(mode = "list", length = length(id_unicos))
  
  # lista de todos lo objetos
  objetos <- ls()
  
  # filtrar objetos de insteres que poblaran la lista
  objetos_lista <- objetos[grepl("split", objetos)]
  
  # conseguir data frames de Gene_ids para poblar lista de gene_is separados
  genes_separados <- mget(objetos_lista)
  
  # crear vector que guadara la informacion de valores maximos
  maximo <- vector()
  
  # ordenar data frames dentro de lista
  for(tabla in 1:length(genes_separados)){
    # poblar vector de valores maximos por cada data frame
    maximo[tabla] <- max(genes_separados[[tabla]]$conteo)
    # de cada df filtrar solo la fila cuyo valor de conteo sea igual al valor maximo encontrado arriba
    genes_separados[[tabla]] <- genes_separados[[tabla]][genes_separados[[tabla]]$conteo == maximo[tabla],]
  }
  
  # convertir lista en data frame
  genes_separados_df <- do.call(rbind, genes_separados)
  
  # unir file con genes separados
  file <- merge(file, 
                genes_separados_df, 
                by = "Gene_ID", 
                all.x = TRUE)
  
  # remover columna de conteo
  file$conteo <- NULL
  
  # ordenar por Blast
  file <- file[order(file$Blast),]
  
  # convertir NAs de Blast en string vacio
  file[is.na(file$Blast), "Blast"] <- ""
  
  ############################################################################
  
  # crear grafica
  volcano_label <- ggplot(data = file,
                          aes(x = log2FC,
                              y = -log10(Padj),
                              col = DEGs,
                              label = Blast)) +
    geom_vline(xintercept = c(-1*corte, corte),
               col = "gray",
               linetype = 'dashed') +
    geom_hline(yintercept = -log10(0.05),
               col = "gray",
               linetype = 'dashed') +
    geom_point(size = 2) +
    scale_color_manual(values = colores,
                       labels = c(paste("Downregulated", " (", down_c, ")", sep = ""),
                                  "Not significant",
                                  paste("Upregulated", " (", up_c, ")", sep = ""))) +
    labs(color = 'Diferentially Expressed Genes',
         x = expression("log"[2]*"FC"), y = expression("-log"[10]*"P-value")) +
    scale_x_continuous(breaks = seq(round(minimo_x) - 2, round(maximo_x) + 2, 2)) +
    ggtitle(paste(tratamiento , "- Treated Protonemal DEGs", sep = " ")) +
    theme_classic() +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1.2)) +
    # geom_text_repel(max.overlaps = Inf, point.padding = 0.2) +
    geom_label_repel(max.overlaps = Inf, col = "black")
  
  # Create a new column "delabel" to de, that will contain the name of the top 30 differentially expressed genes (NA in case they are not)
  # df$delabel <- ifelse(df$gene_symbol %in% head(df[order(df$padj), "gene_symbol"], 30), df$gene_symbol, NA)
  
  # visualizar en R
  print(volcano_label)
  
  # guardar imagen en formatos pre-establecidos
  for(i in formatos) {
    
    # crear y guardar los heatmaps
    match.fun(i)(paste(nuevo_dir, 
                       "/Graficas_volcano/",
                       "Volcano_Labels_Plot_", 
                       paste(tratamiento, collapse = "_"), 
                       ".", i, sep = ""),
                 res = resolucion,
                 width = 3000,
                 height = 2000)
    
    # crear y guardar el heatmpat euclidean
    print(
      
      volcano_label
      
    )
    
    dev.off()
    
  }
  
  # guardar/crear archivo en windows
  write.csv(x = file,
            file = paste(nuevo_dir, "/", "degs_universe_labels_", tratamiento, ".csv", sep = ""),
            row.names = FALSE)
  
  assign(paste("volcano_label_", tratamiento, sep = ""), volcano_label, envir = .GlobalEnv)
  
}