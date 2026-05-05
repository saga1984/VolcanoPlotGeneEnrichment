#' @title Análisis ORA y GSEA (GO/KEGG) con clusterProfiler a partir de DEGs y universo
#'
#' @description
#' Ejecuta un flujo reproducible para enriquecimiento funcional con **clusterProfiler**
#' (GO y KEGG) usando listas de DEGs y universo anotado con TAIR a partir de
#' salidas tipo Trinotate. Genera tablas intermedias/finales y múltiples
#' visualizaciones (dotplot, barplot y enrichment maps) para ORA y GSEA, así
#' como comparaciones por conjuntos (Induced/Repressed) y por clúster.
#'
#' @details
#' - Requiere anotaciones TAIR válidas (prefijo `AT[1-5]G...`) en el universo.
#' - Supone que `degs`, `universo` y `trinotate` son **CSV** con columnas esperadas
#'   (p. ej., `Gene_ID`, `log2FC`, `Pval`, `Padj`, `DEGs`).
#' - Crea subdirectorios `ORA/` y `GSEA/` junto al archivo `degs` para guardar
#'   resultados (`*.csv`) e imágenes (formatos en `formatos`, p. ej. `"jpeg"`).
#' - Usa **org.At.tair.db** (Arabidopsis) y mapea TAIR→ENTREZ para KEGG.
#' - No devuelve objetos en R; escribe resultados en disco.
#'
#' @param degs Ruta al CSV de genes diferencialmente expresados (DEGs).
#' @param universo Ruta al CSV del universo de genes (incluye `log2FC`, `Padj`).
#' @param trinotate Ruta al CSV de anotaciones de Trinotate.
#' @param tratamiento Etiqueta del tratamiento o condición (usada solo en mensajes).
#' @param min_set_go Tamaño mínimo de conjunto (GO) para ORA/GSEA. Default: 5.
#' @param max_set_go Tamaño máximo de conjunto (GO) para ORA/GSEA. Default: 800.
#' @param min_set_kegg Tamaño mínimo de conjunto (KEGG) para ORA/GSEA. Default: 5.
#' @param max_set_kegg Tamaño máximo de conjunto (KEGG) para ORA/GSEA. Default: 800.
#' @param categorias_kegg Número de categorías a mostrar en gráficos KEGG. Default: 22.
#' @param categorias_go Número de categorías a mostrar en gráficos GO. Default: 24.
#' @param categorias_go_cluster Categorías en compareCluster (GO). Default: 22.
#' @param categorias_kegg_cluster Categorías en compareCluster (KEGG). Default: 20.
#' @param alto_go_cluster Alto (px) para gráficos compareCluster GO. Default: 4000.
#' @param ancho_go_cluster Ancho (px) para gráficos compareCluster GO. Default: 3000.
#' @param alto_kegg_cluster Alto (px) para gráficos compareCluster KEGG. Default: 4000.
#' @param ancho_kegg_cluster Ancho (px) para gráficos compareCluster KEGG. Default: 3600.
#' @param alto_go_emap Alto (px) para enrichment map GO. Default: 3000.
#' @param ancho_go_emap Ancho (px) para enrichment map GO. Default: 3500.
#' @param alto_go Alto (px) para gráficos GO (bar/dot). Default: 3400.
#' @param ancho_go Ancho (px) para gráficos GO (bar/dot). Default: 2600.
#' @param alto_go_gsea Alto (px) para gráficos GSEA-GO. Default: 4000.
#' @param ancho_go_gsea Ancho (px) para gráficos GSEA-GO. Default: 3000.
#' @param alto_kegg Alto (px) para gráficos KEGG (bar/dot/emap). Default: 3400.
#' @param ancho_kegg Ancho (px) para gráficos KEGG (bar/dot/emap). Default: 3400.
#' @param sin_duplicados Remover duplicados en anotaciones Trinotate (TRUE/FALSE). Default: FALSE.
#' @param gene_symbol Columna preferida de símbolo génico (p. ej., `"TAIR"`). No cambia el mapeo interno. Default: "TAIR".
#' @param formatos Formato(s) de salida para imágenes (p. ej., `"jpeg"`, `"png"`). Default: "jpeg".
#' @param resolucion Resolución (dpi) para imágenes. Default: 300.
#'
#' @return `NULL`. Los resultados se guardan en disco en subdirectorios `ORA/` y `GSEA/`
#'   (tablas `*.csv` y figuras en los formatos solicitados).
#'
#' @section Salidas principales:
#' - **ORA/ table_ORA_GO*.csv**, **table_ORA_KEGG*.csv** (incluye _final y _Up/_Down).
#' - **GSEA/ table_GSEA_GO*.csv**, **table_GSEA_KEGG*.csv** (incluye _final).
#' - Figuras: barplot/dotplot/emap para GO y KEGG, por conjunto y por clúster.
#'
#' @examples
#' \dontrun{
#' wraper_cluster_profiler(
#'   degs       = "results/DEGs.csv",
#'   universo   = "results/universe.csv",
#'   trinotate  = "annotation/trinotate.csv",
#'   tratamiento = "NaCl_3h",
#'   formatos   = "png",
#'   resolucion = 300
#' )
#' }
#'
#' @seealso
#' \itemize{
#' \item \code{clusterProfiler::enrichGO}, \code{clusterProfiler::enrichKEGG},
#' \code{clusterProfiler::gseGO}, \code{clusterProfiler::gseKEGG},
#' \code{clusterProfiler::compareCluster}
#' \item \code{enrichplot::dotplot}, \code{enrichplot::emapplot}
#' \item \code{org.At.tair.db}, \code{pathview::pathview}
#' }
#'
#' @import data.table
#' @import RColorBrewer
#' @import ggrepel
#' @import pheatmap
#' @import clusterProfiler
#' @import DOSE
#' @import enrichplot
#' @import ggupset
#' @import pathview
#' @import org.At.tair.db
#' @importFrom stringr str_extract
#' @importFrom tidyr separate_rows
#' @importFrom dplyr mutate
#' @importFrom ggplot2 ggplot aes geom_point labs facet_grid scale_x_continuous theme element_text
#'
#' @keywords enrichment ORA GSEA GO KEGG Arabidopsis TAIR
#' @export

# definir funcion
wraper_cluster_profiler <- function(degs,
                                    universo,
                                    trinotate,
                                    tratamiento,
                                    min_set_go = 5,
                                    max_set_go = 800,
                                    min_set_kegg = 5,
                                    max_set_kegg = 800,
                                    categorias_kegg = 22,
                                    categorias_go = 24,
                                    categorias_go_cluster = 22,
                                    categorias_kegg_cluster = 20,
                                    alto_go_cluster = 4000,
                                    ancho_go_cluster = 3000,
                                    alto_kegg_cluster = 4000,
                                    ancho_kegg_cluster = 3600,
                                    alto_go_emap_cluster = 4000,
                                    ancho_go_emap_cluster = 4500,
                                    alto_go = 3400,
                                    ancho_go = 2600,
                                    alto_go_gsea = 4000,
                                    ancho_go_gsea = 3000,
                                    alto_go_emap = 3000,
                                    ancho_go_emap = 3500,
                                    alto_kegg = 3400,
                                    ancho_kegg = 3400,
                                    sin_duplicados = FALSE,
                                    gene_symbol = "TAIR",
                                    formatos = "jpeg",
                                    resolucion = 300){

  # cargar librerias
  library(data.table)
  library(readxl)
  library(writexl)
  library(RColorBrewer)
  library(ggrepel)
  library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
  library(RColorBrewer) # for a colourful plot
  library(pheatmap)
  library(clusterProfiler) # for PEA analysis
  library(DOSE)
  library(enrichplot) # for visualisations
  library(ggupset) # for visualisations
  library(pathview)
  library("org.At.tair.db", character.only = TRUE) # background genes
  
  ################## crear columna de log2FC en trinotate ########################
  
  # volver a home
  setwd("~")
  
  # agregar columna de genes, importar trinotate
  archivo_degs <- read.csv(degs)
  
  # agregar columna de genes, importar trinotate
  archivo_universe <- read.csv(universo)
  
  # agregar columna de genes, importar trinotate
  archivo_trinotate <- read.csv(trinotate)
  
  # filtrar columnas importantes
  archivo_trinotate <- fread(trinotate, select = c(1, 3, 7, 12))
  
  # modificar nombre de columnas
  colnames(archivo_trinotate) <- c("Gene_ID", "Blastx", "Blastp", "Tair")
  
  # eliminar todo lo que no esta entre un = y un ;
  archivo_trinotate$Blastx_clean1 <- str_extract(archivo_trinotate$Blastx, "(?<=\\=)[^;]+")
  archivo_trinotate$Blastp_clean1 <- str_extract(archivo_trinotate$Blastp, "(?<=\\=)[^;]+")
  archivo_trinotate$Blastx_clean2 <- gsub("\\_.*$", "", archivo_trinotate$Blastx)
  archivo_trinotate$Blastp_clean2 <- gsub("\\_.*$", "", archivo_trinotate$Blastp)
  
  # crear una nueva columna de blast (considerando tanto blastx como blastp)
  archivo_trinotate$Blast1 <- ifelse(is.na(archivo_trinotate$Blastp_clean1), 
                                     archivo_trinotate$Blastx_clean1,
                                     archivo_trinotate$Blastp_clean1)
  
  # crear una nueva columna de blast (considerando tanto blastx como blastp)
  archivo_trinotate$Blast2 <- ifelse(archivo_trinotate$Blastp_clean2 != ".", 
                                     archivo_trinotate$Blastp_clean2, 
                                     archivo_trinotate$Blastx_clean2)
  
  
  # filtrar solo columnas importantes
  archivo_trinotate <- archivo_trinotate[,c("Gene_ID", "Tair", "Blast1", "Blast2")]
  
  # limpiar columna Tair y columna Blast2
  archivo_trinotate$Tair <- gsub("[^:]*:", "", archivo_trinotate$Tair)
  archivo_trinotate$Tair <- gsub("^[0-9]*", "", archivo_trinotate$Tair)
  archivo_trinotate$Tair <- gsub("\\.", NA, archivo_trinotate$Tair)
  archivo_trinotate$Tair[archivo_trinotate$Tair == ""] <- NA
  archivo_trinotate$Blast2 <- gsub("\\.", NA, archivo_trinotate$Blast2)
  
  # Eliminar filas vacias
  archivo_trinotate <- na.omit(archivo_trinotate)
  
  # remover duplicados
  archivo_trinotate <- archivo_trinotate[!duplicated(archivo_trinotate),]
  
  # filtrar archivo trinotate (anotado) por Gene IDs existentes en archivo degs
  archivo_trinotate_universe <- archivo_trinotate[archivo_trinotate$Gene_ID %in% archivo_universe$Gene_ID,]
  
  #################### archivo final de universo #################################
  
  # unir archivos por Gene ID AKA Gene Symbol
  universe_df <- merge(archivo_trinotate, archivo_universe, by = "Gene_ID")
  
  # Remover duplicados por default (se mantienen duplicados)
  if(sin_duplicados){
    archivo_trinotate <- archivo_trinotate[!duplicated(archivo_trinotate),]
  } else {
    archivo_trinotate <- archivo_trinotate
  }
  
  # ordenar en orden decreciente (necesario para cluster profiler)
  universe_df_ordered_logfc <- universe_df[order(universe_df$log2FC, 
                                                 decreasing = TRUE),]
  
  # remove everything that has not TAIR ID
  universe_df$Tair[!grepl("^AT[1-5]G", universe_df$Tair)] <- NA
  universe_df <- universe_df[!is.na(universe_df$Tair),]
  
  # cambiar nombre de columna Tair
  colnames(universe_df)[4] <- "Gene_Symbol"
  colnames(universe_df)[2] <- "TAIR"
  
  ################################################################################
  ################################# inputs ORA ###################################
  ################################################################################
  
  ################################# universo #####################################
  
  # ORDENAR POR LOG2fc PARA ora
  universe_ora_ordered <- universe_df[order(universe_df$log2FC, decreasing = TRUE),]
  
  # crear vector (from object in row 180)
  universe_ora_vec <- universe_ora_ordered$log2FC
  
  # nombrar vector
  names(universe_ora_vec) <- universe_ora_ordered$TAIR
  
  cat("\nvector de DEGs universo\n")
  # inducidos
  cat("Primeros valores", head(universe_ora_vec), sep = "\n", "\n")
  # reprimidos
  cat("Ultimos valores", tail(universe_ora_vec), sep = "\n", "\n")
  
  ########################## genes significativos ################################
  
  # filtrar por valores con valor de TAIR
  degs_sig_ora <- universe_ora_ordered[universe_ora_ordered$Gene_ID %in% archivo_degs$Gene_ID,]
  
  # definir valores
  degs_sig_ora_vec <- degs_sig_ora$log2FC
  
  # definir nombre de vector
  names(degs_sig_ora_vec) <- degs_sig_ora$TAIR
  
  # visualizar vector de nombres
  cat("\nvector de DEGs significativos\n")
  cat("Primeros valores", head(degs_sig_ora_vec), sep = "\n", "\n")
  cat("Ultimos valores", tail(degs_sig_ora_vec), sep = "\n", "\n")
  
  ################################################################################
  
  # directorio de resultados ora
  ora_dir <- paste(dirname(degs), "/ORA/", sep = "")
  
  # guardar el objeto como tabla
  dir.create(ora_dir)
  
  ################################################################################
  ############################### inputs GSEA ####################################
  ################################################################################
  
  # crear columna de ranked
  universe_df[, "Ranking"] <- universe_df$log2FC * -log10(universe_df$Padj)
  
  # ordenar por "Ranking"
  universe_gsea_ordered <- universe_df[order(universe_df$Ranking, decreasing = TRUE),]
  
  # crear vector de universo
  universe_gsea_vec <- universe_gsea_ordered$Ranking
  
  # nombrar vector de universo
  names(universe_gsea_vec) <- universe_gsea_ordered$TAIR
  
  # mensaje a consola de vector de GSEA
  cat("\nVector de Ranking de GSEA\n")
  # verificar vector de universo
  cat("Primeros valores", head(universe_gsea_vec), sep = "\n", "\n")
  # verificar vector de universo
  cat("Ultimos valores" ,tail(universe_gsea_vec), sep = "\n", "\n")
  
  ################################################################################
  
  # directorio de resultados gsea
  gsea_dir <- paste(dirname(degs), "/GSEA/", sep = "")
  
  # guardar el objeto como tabla
  dir.create(gsea_dir)
  
  ################################################################################
  ############################### enrichGO #######################################
  ################################################################################
  
  cat("\n\t\t============= 'COMENZANDO ANALISIS POR ORA' =============\n")
  
  # run cluster profiler
  oraGO <- enrichGO(gene = names(degs_sig_ora_vec),
                    universe = names(universe_ora_vec),
                    OrgDb = org.At.tair.db,
                    keyType = "TAIR",
                    readable = TRUE,
                    ont = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.10)
  
  # verificar keyTypes
  keytypes(org.At.tair.db)
  
  head(keys(org.At.tair.db, keytype = "TAIR"))
  
  # convertir a tabla
  table_GO_enrich <- as.data.frame(oraGO)
  
  # guardar el objeto como tabla  "nuevo bloque"
  write_xlsx(x = table_GO_enrich,
             paste(ora_dir, "table_ORA_GO.xlsx",sep = ""),
             col_names = TRUE, 
             format_headers = TRUE)
  
  # eliminar columnas innecesarias "nuevo bloque"
  table_GO_enrich <- table_GO_enrich[, c("ID",
                                         "Description",
                                         "ONTOLOGY",
                                         "pvalue",
                                         "p.adjust",
                                         "qvalue",
                                         "geneID",
                                         "Count",
                                         "GeneRatio",
                                         "BgRatio")]
  
  # separar columnas
  table_GO_enrich_long <- table_GO_enrich %>%
    separate_rows(geneID, sep = "/")
  
  # TAIR to Blast
  gene2symbol <- oraGO@gene2Symbol
  
  # agregar columna de TAIR
  Tair <- names(gene2symbol)
  
  # agregar columna de gene symbol (Blast)
  Symbol <- gene2symbol
  
  # poblar data frame
  gene2symbol_df <- data.frame(
    TAIR = Tair,
    geneID = Symbol)
  
  # eliminar NAs del data frame
  gene2symbol_df <- na.omit(gene2symbol_df)
  
  # unir data frames
  table_GO_enrich_final_tmp <- merge(table_GO_enrich_long,
                                     gene2symbol_df,
                                     by = "geneID")
  
  # agregar ID Trinity
  table_GO_enrich_final <- merge(table_GO_enrich_final_tmp,
                                 degs_sig_ora[, c("TAIR", 
                                                  "Gene_ID", 
                                                  "log2FC",
                                                  "Pval",       
                                                  "Padj",
                                                  "DEGs")],
                                 by = "TAIR")
  
  # reordenar columnas "nuevo bloque"
  table_GO_enrich_final <- table_GO_enrich_final[, c("ID",
                                                     "TAIR",
                                                     "Gene_ID",
                                                     "geneID",
                                                     "Description",
                                                     "ONTOLOGY",
                                                     "Count",
                                                     "GeneRatio",
                                                     "BgRatio",
                                                     "pvalue",
                                                     "p.adjust",
                                                     "qvalue",
                                                     "log2FC",
                                                     "Pval",       
                                                     "Padj",
                                                     "DEGs")]
  
  # ordenar por description
  table_GO_enrich_final <- table_GO_enrich_final[order(table_GO_enrich_final[, "Description"]),]
  
  # guardar el objeto como tabla "nuevo bloque"
  write_xlsx(x = table_GO_enrich_final,
             paste(ora_dir, "table_ORA_GO_final.xlsx",sep = ""), 
             col_names = TRUE, 
             format_headers = TRUE)
  
  ###################### bar
  
  # guardar imagen
  for(i in formatos) {
    
    # Usar tryCatch para atrapar errores y controlar su salida
    # que el error no detenga la funcion
    tryCatch({
      
      # crear y guardar los heatmaps
      match.fun(i)(paste(ora_dir, "ORA_GO_Bar_Ontology.", i, sep = ""),
                   res = resolucion,
                   width = ancho_go,
                   height = alto_go)
      
      # crear y guardar el heatmpat
      print(
        
        barplot(oraGO,
                showCategory = categorias_go) + 
          facet_grid(~ONTOLOGY) +
          labs(title = "ORA: GO Enrichment")
        
      )
      
      dev.off()
      
    }, error = function(e) {
      # imprimir el mensaje de error
      message("ocurrio un error durante la graficacion de ORA GO", e$message)
    })
    
  }
  
  ##################### dot
  
  # guardar imagen
  for(i in formatos) {
    
    # Usar tryCatch para atrapar errores y controlar su salida
    # que el error no detenga la funcion
    tryCatch({
      
      # crear y guardar los heatmaps
      match.fun(i)(paste(ora_dir, "ORA_GO_Dot_Ontology.", i, sep = ""),
                   res = resolucion,
                   width = ancho_go,
                   height = alto_go)
      
      # crear y guardar el heatmpat
      print(
        
        dotplot(oraGO,
                showCategory = categorias_go) +
          facet_grid(~ONTOLOGY) +
          labs(title = "ORA: GO Enrichment")
        
      )
      
      dev.off()
      
    }, error = function(e) {
      # imprimir el mensaje de error
      message("ocurrio un error durante la graficacion de GO-BP-Down: ", e$message)
    })
    
  }
  
  #################### separar por inducidos y reprimidos ########################
  
  # filtrar positivos == inducidos; negativos == reprimidos
  Activated_ora_go <- names(degs_sig_ora_vec[degs_sig_ora_vec > 0])
  Suppressed_ora_go <- names(degs_sig_ora_vec[degs_sig_ora_vec < 0])
  
  ###################### crear lista de activados ################################
  oraGO_up <- enrichGO(gene = Activated_ora_go,
                       universe = names(universe_ora_vec),
                       OrgDb = org.At.tair.db,
                       keyType = "TAIR",
                       readable = TRUE,
                       ont = "ALL",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.10)
  
  # convertir a data frame y agregar columna de activado
  activated_ora_go <- as.data.frame(oraGO_up) %>%
    mutate(GeneSet = "Activated")
  
  # guardar el objeto como tabla "nuevo bloque"
  write_xlsx(x = activated_ora_go, 
             paste(ora_dir, "table_ORA_GO_Up.xlsx",sep = ""),
             col_names = TRUE, 
             format_headers = TRUE)
  
  # TAIR to Blast
  activated_gene2symbol <- oraGO_up@gene2Symbol
  
  # agregar columna de TAIR
  activated_Tair <- names(activated_gene2symbol)
  
  # agregar columna de gene symbol (Blast)
  activated_Symbol <- activated_gene2symbol
  
  # poblar data frame
  gene2symbol_activated <- data.frame(
    TAIR = activated_Tair,
    geneID = activated_Symbol)
  
  # eliminar NAs del data frame
  gene2symbol_activated <- na.omit(gene2symbol_activated)
  
  # eliminar columnas innecesarias "nuevo bloque"
  activated_ORA_GO <- activated_ora_go[, c("ID",
                                           "Description",
                                           "ONTOLOGY",
                                           "pvalue",
                                           "p.adjust",
                                           "qvalue",
                                           "geneID",
                                           "GeneSet",
                                           "Count",
                                           "GeneRatio",
                                           "BgRatio")]
  
  # separar columnas
  activated_ORA_GO_long <- activated_ORA_GO %>%
    separate_rows(geneID, sep = "/")
  
  # unir data frames
  activated_ORA_GO_final_tmp <- merge(activated_ORA_GO_long,
                                      gene2symbol_activated,
                                      by = "geneID")
  
  # agregar ID Trinity
  activated_ORA_GO_final <- merge(activated_ORA_GO_final_tmp,
                                  degs_sig_ora[, c("TAIR", 
                                                   "Gene_ID", 
                                                   "log2FC",
                                                   "Pval",       
                                                   "Padj",
                                                   "DEGs")],
                                  by = "TAIR")
  
  # reordenar columnas "nuevo bloque"
  activated_ORA_GO_final <- activated_ORA_GO_final[, c("ID",
                                                       "TAIR",
                                                       "Gene_ID",
                                                       "geneID",
                                                       "Description",
                                                       "ONTOLOGY",
                                                       "Count",
                                                       "GeneRatio",
                                                       "BgRatio",
                                                       "GeneSet",
                                                       "pvalue",
                                                       "p.adjust",
                                                       "qvalue",
                                                       "log2FC",
                                                       "Pval",       
                                                       "Padj",
                                                       "DEGs")]
  
  # ordenar por description
  activated_ORA_GO_final <- activated_ORA_GO_final[order(activated_ORA_GO_final[, "Description"]),]
  
  # guardar el objeto como tabla "nuevo bloque"
  write_xlsx(x = activated_ORA_GO_final, 
             paste(ora_dir, "table_ORA_GO_Up_final.xlsx",sep = ""),
             col_names = TRUE, 
             format_headers = TRUE)
  
  ##################### crear lista de reprimidos ################################
  oraGO_down <- enrichGO(gene = Suppressed_ora_go,
                         universe = names(universe_ora_vec),
                         OrgDb = org.At.tair.db,
                         keyType = "TAIR",
                         readable = TRUE,
                         ont = "ALL",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.10)
  
  # convertir a data frame y agregar columna de reprimido
  suppressed_ora_go <- as.data.frame(oraGO_down) %>%
    mutate(GeneSet = "Suppressed")
  
  # guardar el objeto como tabla "nuevo bloque"
  write_xlsx(x = suppressed_ora_go,
             paste(ora_dir, "table_ORA_GO_Down.xlsx",sep = ""),
             col_names = TRUE, 
             format_headers = TRUE)
  
  # TAIR to Blast
  suppressed_gene2symbol <- oraGO_down@gene2Symbol
  
  # agregar columna de TAIR
  suppressed_Tair <- names(suppressed_gene2symbol)
  
  # agregar columna de gene symbol (Blast)
  suppressed_Symbol <- suppressed_gene2symbol
  
  # poblar data frame
  gene2symbol_suppressed <- data.frame(
    TAIR = suppressed_Tair,
    geneID = suppressed_Symbol)
  
  # eliminar NAs del data frame
  gene2symbol_suppressed <- na.omit(gene2symbol_suppressed)
  
  # eliminar columnas innecesarias "nuevo bloque"
  suppressed_ORA_GO <- suppressed_ora_go[, c("ID",
                                             "Description",
                                             "ONTOLOGY",
                                             "pvalue",
                                             "p.adjust",
                                             "qvalue",
                                             "geneID",
                                             "GeneSet",
                                             "Count", 
                                             "GeneRatio",
                                             "BgRatio")]
  
  # separar columnas
  suppressed_ORA_GO_long <- suppressed_ORA_GO %>%
    separate_rows(geneID, sep = "/")
  
  # unir data frames
  suppressed_ORA_GO_final_tmp <- merge(suppressed_ORA_GO_long,
                                       gene2symbol_suppressed,
                                       by = "geneID")
  
  # agregar ID Trinity
  suppressed_ORA_GO_final <- merge(suppressed_ORA_GO_final_tmp,
                                   degs_sig_ora[, c("TAIR", 
                                                    "Gene_ID", 
                                                    "log2FC",
                                                    "Pval",       
                                                    "Padj",
                                                    "DEGs")],
                                   by = "TAIR")
  
  # reordenar columnas "nuevo bloque"
  suppressed_ORA_GO_final <- suppressed_ORA_GO_final[, c("ID",
                                                         "TAIR",
                                                         "Gene_ID",
                                                         "geneID",
                                                         "Description",
                                                         "ONTOLOGY",
                                                         "Count", 
                                                         "GeneRatio",	
                                                         "BgRatio",
                                                         "GeneSet",
                                                         "pvalue",
                                                         "p.adjust",
                                                         "qvalue",
                                                         "log2FC",
                                                         "Pval",       
                                                         "Padj",
                                                         "DEGs")]
  
  # ordenar por description
  suppressed_ORA_GO_final <- suppressed_ORA_GO_final[order(suppressed_ORA_GO_final[, "Description"]),]
  
  # guardar el objeto como tabla "nuevo bloque"
  write_xlsx(x = suppressed_ORA_GO_final,
             paste(ora_dir, "table_ORA_GO_Down_final.xlsx",sep = ""),
             col_names = TRUE, 
             format_headers = TRUE)
  
  ########################### graficar ###########################################
  
  # ordenar por pajustado
  suppressed_ora_go <- suppressed_ora_go[order(suppressed_ora_go$p.adjust),]
  activated_ora_go <- activated_ora_go[order((activated_ora_go$p.adjust)),]
  
  # unir data frames
  resultados_combinados_ora_go <- rbind(activated_ora_go,
                                        suppressed_ora_go)
  
  ################## guardar archivos combinados Up Down #########################
  
  # guardar el objeto como tabla "nuevo bloque"
  write_xlsx(x = resultados_combinados_ora_go,
             paste(ora_dir, "table_ORA_GO_UpDown.xlsx", sep = ""),
             col_names = TRUE, 
             format_headers = TRUE)
  
  # unir archivos finales
  resultados_combinados_ora_go_final <- rbind(activated_ORA_GO_final,
                                              suppressed_ORA_GO_final)
  
  # guardar el objeto como tabla "nuevo bloque"
  write_xlsx(x = resultados_combinados_ora_go_final, 
             paste(ora_dir, "table_ORA_GO_UpDown_final.xlsx", sep = ""),
             col_names = TRUE, 
             format_headers = TRUE)
  
  ########################### greficar continuacion ##############################
  
  # crear col de gene ratio numerica
  for(i in 1:nrow(resultados_combinados_ora_go)){
    resultados_combinados_ora_go[i, "GeneRatioN"] <- eval(
      parse(text = resultados_combinados_ora_go[i, "GeneRatio"]))
  }
  
  # Filtrar mejores resultados para activados y suprimidos por p.adjust
  # up
  resultados_up_ora_go_tmp <- resultados_combinados_ora_go[resultados_combinados_ora_go$GeneSet == "Activated", ]
  resultados_up_ora_go <- resultados_up_ora_go_tmp[(1:(categorias_go/2)), ]
  resultados_up_ora_go <- na.omit(resultados_up_ora_go)
  
  # down
  resultados_down_ora_go_tmp <- resultados_combinados_ora_go[resultados_combinados_ora_go$GeneSet == "Suppressed", ]
  resultados_down_ora_go <- resultados_down_ora_go_tmp[(1:(categorias_go/2)),]
  resultados_down_ora_go <- na.omit(resultados_down_ora_go)
  
  # unir por filas
  resultados_top_ora_go <- rbind(resultados_up_ora_go,
                                 resultados_down_ora_go)
  
  ################################# graficar #####################################
  
  # guardar imagen
  for(i in formatos) {
    
    # Usar tryCatch para atrapar errores y controlar su salida
    # que el error no detenga la funcion
    tryCatch({
      
      # crear y guardar los heatmaps
      match.fun(i)(paste(ora_dir, "ORA_GO_Dot_UpDown.", i, sep = ""),
                   res = resolucion,
                   width = ancho_go,
                   height = alto_go)
      
      # crear y guardar el heatmpat
      print(
        
        ggplot(resultados_top_ora_go, 
               aes(x = GeneRatioN, 
                   y = Description)) +
          geom_point(aes(size = Count, 
                         color = p.adjust)) +
          facet_grid(GeneSet ~ ONTOLOGY,) +
          scale_color_gradient(high = "blue", low = "red") +
          theme_dose() +
          theme(axis.text.y = element_text(size = 12),
                strip.text = element_text(size = 13),
                title = element_text(size = 14)) +
          labs(title = "ORA: GO Enrichment", 
               y = "", 
               x = "Gene Ratio" ) +
          scale_x_continuous(n.breaks = 3)
        
      )
      
      dev.off()
      
    }, error = function(e) {
      # imprimir el mensaje de error
      message("ocurrio un error durante la graficacion de GO-Dot-All: ", e$message)
    })
    
  }
  
  ######################## todos ##############################
  
  # guardar imagen
  for(i in formatos) {
    
    # Usar tryCatch para atrapar errores y controlar su salida
    # que el error no detenga la funcion
    tryCatch({
      
      # crear y guardar los heatmaps
      match.fun(i)(paste(ora_dir, "ORA_GO_Dot_All.", i, sep = ""),
                   res = resolucion,
                   width = 4000,
                   height = 8400)
      
      # crear y guardar el heatmpat
      print(
        
        ggplot(resultados_combinados_ora_go, 
               aes(x = GeneRatioN, 
                   y = Description)) +
          geom_point(aes(size = Count, 
                         color = p.adjust)) +
          facet_grid(GeneSet ~ ONTOLOGY) +
          scale_color_gradient(high = "blue", low = "red") +
          theme_dose() +
          theme(axis.text.y = element_text(size = 14),
                axis.text.x = element_text(size = 14),
                axis.title.x = element_text(size = 16),
                strip.text = element_text(size = 16),
                title = element_text(size = 18),
                legend.title = element_text(size = 16),
                legend.text = element_text(size = 14)) +
          labs(title = "ORA: GO Enrichment", 
               y = "", 
               x = "Gene Ratio" ) +
          scale_x_continuous(n.breaks = 3)
        
      )
      
      dev.off()
      
    }, error = function(e) {
      # imprimir el mensaje de error
      message("ocurrio un error durante la graficacion de GO-Dot-All: ", e$message)
    })
    
  }
  
  ################################################################################
  ########################## EnrichGO: Compare Cluster ###########################
  ################################################################################
  
  #################### separar por inducidos y reprimidos ########################
  
  # crear lista de genes
  gene_lists <- list(
    Induced = names(degs_sig_ora_vec[degs_sig_ora_vec > 0]),
    Repressed = names(degs_sig_ora_vec[degs_sig_ora_vec < 0])
  )
  
  # comparar setd de inducidos y reprimidos
  compare_results_GO <- compareCluster(
    geneCluster = gene_lists,
    fun = "enrichGO",
    OrgDb = org.At.tair.db,
    keyType = "TAIR",
    ont = "ALL",
    universe = names(universe_ora_vec),  # Add the universe argument here
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.10)
  
  
  # verificar keyTypes
  keytypes(org.At.tair.db)
  head(keys(org.At.tair.db, keytype = "TAIR"))
  
  # convertir a tabla
  table_compare_results_go <- as.data.frame(compare_results_GO)
  
  # guardar el objeto como tabla "nuevo bloque"
  write_xlsx(x = table_compare_results_go,
             paste(ora_dir, "table_ORA_ClusterCompare_GO.xlsx",sep = ""),
             col_names = TRUE, 
             format_headers = TRUE)
  
  # eliminar columnas innecesarias "nuevo bloque"
  table_compare_results_go <- table_compare_results_go[, c("ID",
                                                           "Description",
                                                           "ONTOLOGY",
                                                           "Cluster",
                                                           "pvalue",
                                                           "p.adjust",
                                                           "qvalue",
                                                           "geneID",
                                                           "Count",
                                                           "GeneRatio",
                                                           "BgRatio")]
  
  # separar columnas
  table_compare_results_go_long <- table_compare_results_go %>%
    separate_rows(geneID, sep = "/")
  
  # change colname GeneID
  colnames(table_compare_results_go_long)[8] <- "TAIR"
  
  # agregar ID Trinity
  table_compare_results_go_final <- merge(table_compare_results_go_long,
                                          degs_sig_ora[, c("TAIR", 
                                                           "Gene_ID", 
                                                           "log2FC",
                                                           "Pval",       
                                                           "Padj",
                                                           "DEGs")],
                                          by = "TAIR")
  
  # reordenar columnas "nuevo bloque"
  table_compare_results_go_final <- table_compare_results_go_final[, c("ID",
                                                                       "TAIR",
                                                                       "Gene_ID",
                                                                       "Description",
                                                                       "ONTOLOGY",
                                                                       "Cluster",
                                                                       "Count",
                                                                       "GeneRatio",
                                                                       "BgRatio",
                                                                       "pvalue",
                                                                       "p.adjust",
                                                                       "qvalue",
                                                                       "log2FC",
                                                                       "Pval",       
                                                                       "Padj",
                                                                       "DEGs")]
  
  # ordenar por description
  table_compare_results_go_final <- table_compare_results_go_final[order(table_compare_results_go_final[, "Description"]),]
  
  # guardar el objeto como tabla
  write_xlsx(x = table_compare_results_go_final,
             paste(ora_dir, "table_ORA_ClusterCompare_GO_final.xlsx",sep = ""),
             col_names = TRUE, 
             format_headers = TRUE)
  
  ###################### dot
  
  # guardar imagen
  for(i in formatos) {
    
    # Usar tryCatch para atrapar errores y controlar su salida
    # que el error no detenga la funcion
    tryCatch({
      
      # crear y guardar los heatmaps
      match.fun(i)(paste(ora_dir, "ORA_GO_ClusterCompare_Dot_Ontology.", i, sep = ""),
                   res = resolucion,
                   width = ancho_go_cluster,
                   height = alto_go_cluster)
      
      
      # crear y guardar el heatmpat
      print(
        
        dotplot(compare_results_GO,
                showCategory = categorias_go_cluster,
                split = ".sign") +
          facet_grid(~ONTOLOGY) +
          labs(title = "ORA Clusters: GO Enrichment", x = "")
        
      )
      
      dev.off()
      
    }, error = function(e) {
      # imprimir el mensaje de error
      message("ocurrio un error durante la graficacion de Cluster Compare ORA GO", e$message)
    })
    
  }
  
  ################################################################################
  ############################# KEGG (Pathways) ##################################
  ################################################################################
  
  # obtener los Gene_IDs para la funcion enrichKEGG(), (podemos perder algunos genes aquí)
  # pero es un paso innecesario dado a los keyTypes limitados!
  ORA_TAIR_2_ENTREZID <-bitr(names(universe_ora_vec),
                             fromType = "TAIR",
                             toType = "ENTREZID",
                             OrgDb = org.At.tair.db)
  
  ################################################################################
  
  "nuevo bloque"
  # ENTREZ to EC
  ORA_TAIR_2_EC <-bitr(names(universe_ora_vec),
                       fromType = "TAIR",
                       toType = "ENZYME",
                       OrgDb = org.At.tair.db)
  
  # join ENTREZ y ENZYME by TAIR
  ENTREZ_ENZYME_TAIR <- merge(ORA_TAIR_2_ENTREZID,
                              ORA_TAIR_2_EC,
                              by = "TAIR",
                              all.x = TRUE)
  
  # add ENTREZ and ENZYME IDs column
  universe_ora_kegg <- merge(universe_ora_ordered,
                             ENTREZ_ENZYME_TAIR,
                             by = "TAIR",
                             all.x = TRUE)
  
  # reordenar nombres de columnas
  universe_ora_kegg <- universe_ora_kegg[,c("TAIR",
                                            "ENTREZID",
                                            "ENZYME",
                                            "Gene_ID",
                                            "Blast1",
                                            "Gene_Symbol",
                                            "log2FC",
                                            "Pval",
                                            "Padj",
                                            "DEGs")]
  
  # ordenar por ranking
  universe_ora_kegg <- universe_ora_kegg[order(universe_ora_kegg$log2FC, decreasing = TRUE),]
  
  # crear vector de universo
  universe_ora_kegg_vec_entrez <- universe_ora_kegg$log2FC
  universe_ora_kegg_vec_tair <- universe_ora_kegg$log2FC # para pathview()
  
  # nombrar vector de universo
  names(universe_ora_kegg_vec_entrez) <- universe_ora_kegg$ENTREZID
  names(universe_ora_kegg_vec_tair) <- universe_ora_kegg$TAIR # para pathview()
  
  # verificar vector de universo
  cat("\nvector KEGG universo\n")
  cat("ENTREZ:", head(names(universe_ora_kegg_vec_entrez)))
  cat("TAIR:", head(names(universe_ora_kegg_vec_tair)))
  cat("Primeros valores ENTREZ:", head(universe_ora_kegg_vec_entrez), sep = "\n")
  cat("Primeros valores TAIR:", head(universe_ora_kegg_vec_tair), sep = "\n")
  cat("Ultimos valores ENTREZ:", tail(universe_ora_kegg_vec_entrez), sep = "\n")
  cat("Ultimos valores TAIR:", tail(universe_ora_kegg_vec_tair), sep = "\n")
  
  ########### degs significativos kegg ###########
  
  # class(universe_ora_kegg$Gene_ID)
  # class(archivo_degs$Gene_ID)
  
  # length(intersect(universe_ora_kegg$Gene_ID,
  #                  archivo_degs$Gene_ID))
  
  cat("Numero de Genes en universe:", length(unique(universe_ora_kegg$Gene_ID)), "\n")
  cat("Numero de Genes en archivo_degs:", length(unique(archivo_degs$Gene_ID)), "\n")
  cat("Numero de Genes en Intersección:", length(intersect(universe_ora_kegg$Gene_ID,
                                                           archivo_degs$Gene_ID)), "\n")
  
  # filtrar por valores con valor de TAIR
  degs_sig_ora_kegg <- universe_ora_kegg[universe_ora_kegg$Gene_ID %in% archivo_degs$Gene_ID,]
  
  # definir valores
  degs_sig_ora_kegg_vec_entrez <- degs_sig_ora_kegg$log2FC
  degs_sig_ora_kegg_vec_tair <- degs_sig_ora_kegg$log2FC
  
  # definir nombre de vector
  names(degs_sig_ora_kegg_vec_entrez) <- degs_sig_ora_kegg$ENTREZID
  names(degs_sig_ora_kegg_vec_tair) <- degs_sig_ora_kegg$TAIR
  
  # visualizar vector de nombres
  cat("\nvector KEGG estadisticamente significativos\n")
  cat("Primeros valores ENTREZ:", head(degs_sig_ora_kegg_vec_entrez))
  cat("Primeros valores TAIR:", head(degs_sig_ora_kegg_vec_tair))
  cat("Ultimos valores ENTREZ:", tail(degs_sig_ora_kegg_vec_entrez))
  cat("Ultimos valores TAIR:", tail(degs_sig_ora_kegg_vec_tair))
  
  ################################################################################
  ############################## enrichKEGG ######################################
  ################################################################################
  
  # headers from ENTREZ: "818236" "838632" "824262" "838633" "824470" "824470"
  oraKEGG <- enrichKEGG(universe = names(universe_ora_kegg_vec_entrez),
                        gene = names(degs_sig_ora_kegg_vec_entrez),
                        minGSSize = min_set_kegg,
                        maxGSSize = max_set_kegg,
                        organism = "ath",
                        keyType = "ncbi-geneid",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05)
  
  
  # remover especie de filas
  oraKEGG@result$Description <- gsub(
    "- Arabidopsis thaliana \\(thale cress\\)$", "",
    oraKEGG@result$Description)
  
  # convertir a tabla
  table_ORA_KEGG <- as.data.frame(oraKEGG)
  
  # cambiar nombre de columna de IDs de TAIR
  colnames(table_ORA_KEGG)[8] <- "ENTREZID"
  
  # guardar el objeto como tabla "nuevo bloque"
  write_xlsx(x = table_ORA_KEGG,
             paste(ora_dir, "table_ORA_KEGG.xlsx",sep = ""),
             col_names = TRUE, 
             format_headers = TRUE)
  
  # eliminar columnas innecesarias "nuevo bloque"
  table_ORA_KEGG_tmp <- table_ORA_KEGG[, c("ID",
                                           "Description",
                                           "ENTREZID",
                                           "Count", 
                                           "GeneRatio",
                                           "BgRatio",
                                           "pvalue",
                                           "p.adjust",
                                           "qvalue")]
  
  # separar columnas
  table_ORA_KEGG_long <- table_ORA_KEGG_tmp %>% 
    separate_rows(ENTREZID, sep = "/")
  
  # unir datos de entrada con datos de salida en un super data frame
  table_ORA_KEGG_final <- merge(table_ORA_KEGG_long, 
                                degs_sig_ora_kegg,
                                by = "ENTREZID")
  
  # remover duplicados
  table_ORA_KEGG_final <- table_ORA_KEGG_final[!duplicated(table_ORA_KEGG_final),]
  
  # ordenar columnas por Trinity ID
  table_ORA_KEGG_final <- table_ORA_KEGG_final[order(table_ORA_KEGG_final$Gene_ID),]
  
  # cambiar orden de columnas "nuevo bloque"
  table_ORA_KEGG_final <- table_ORA_KEGG_final[,c("Gene_ID",
                                                  "TAIR",
                                                  "ENTREZID",
                                                  "ENZYME",
                                                  "ID",
                                                  "Description", 
                                                  "GeneRatio",
                                                  "BgRatio",
                                                  "pvalue",
                                                  "p.adjust",
                                                  "qvalue",
                                                  "Blast1",
                                                  "Gene_Symbol",
                                                  "log2FC",
                                                  "DEGs",
                                                  "Pval",
                                                  "Padj")]
  
  # guardar el objeto como tabla "nuevo bloque"
  write_xlsx(x = table_ORA_KEGG_final,
             paste(ora_dir, "table_ORA_KEGG_final.xlsx",sep = ""),
             col_names = TRUE, 
             format_headers = TRUE)
  
  ##################### dot plot 
  
  # guardar imagen
  for(i in formatos) {
    
    # Usar tryCatch para atrapar errores y controlar su salida
    # que el error no detenga la funcion
    tryCatch({
      
      # crear y guardar los heatmaps
      match.fun(i)(paste(ora_dir, "ORA_KEGG_Dot.", i, sep = ""),
                   res = resolucion,
                   width = ancho_kegg,
                   height = alto_kegg)
      
      # crear y guardar el heatmpat
      print(
        
        # dot plot of KEGG
        dotplot(oraKEGG,
                showCategory = categorias_kegg) + 
          labs(title = "ORA: Enriched Pathways")
        
      )
      
      dev.off()
      
    }, error = function(e) {
      # imprimir el mensaje de error
      message("ocurrio un error durante la graficacion de ORA-Dot: ", e$message)
    })
    
  }
  
  ##################### bar plot 
  
  # guardar imagen
  for(i in formatos) {
    
    # Usar tryCatch para atrapar errores y controlar su salida
    # que el error no detenga la funcion
    tryCatch({
      
      # crear y guardar los heatmaps
      match.fun(i)(paste(ora_dir, "ORA_KEGG_Bar.", i, sep = ""),
                   res = resolucion,
                   width = ancho_kegg,
                   height = alto_kegg)
      
      # crear y guardar el heatmpat
      print(
        
        # dot plot of KEGG
        barplot(oraKEGG,
                showCategory = categorias_kegg) + 
          labs(title = "ORA: Enriched Pathways")
        
      )
      
      dev.off()
      
    }, error = function(e) {
      # imprimir el mensaje de error
      message("ocurrio un error durante la graficacion de ORA-Dot: ", e$message)
    })
    
  }
  
  ################################################################################
  ########################## EnrichKEGG: Compare Cluster ###########################
  ################################################################################
  
  #################### separar por inducidos y reprimidos ########################
  
  # crear lista de genes
  gene_lists_kegg <- list(
    Induced_kegg = names(degs_sig_ora_kegg_vec_entrez[degs_sig_ora_kegg_vec_entrez > 0]),
    Repressed_kegg = names(degs_sig_ora_kegg_vec_entrez[degs_sig_ora_kegg_vec_entrez < 0])
  )
  
  # comparar setd de inducidos y reprimidos
  compare_results_KEGG <- compareCluster(
    geneCluster = gene_lists_kegg,
    fun = "enrichKEGG",
    universe = names(universe_ora_kegg_vec_entrez),
    minGSSize = min_set_kegg,
    maxGSSize = max_set_kegg,
    organism = "ath", 
    keyType = "ncbi-geneid",
    pvalueCutoff = 0.05)
  
  # remover especie de filas
  compare_results_KEGG@compareClusterResult$Description <- gsub(
    "- Arabidopsis thaliana \\(thale cress\\)$", "",
    compare_results_KEGG@compareClusterResult$Description)
  
  # convertir a tabla
  table_ORA_ClusterCompare_KEGG <- as.data.frame(compare_results_KEGG)
  
  # cambiar nombre de columna de IDs de TAIR
  colnames(table_ORA_ClusterCompare_KEGG)[9] <- "ENTREZID"
  
  # guardar el objeto como tabla "nuevo bloque"
  write_xlsx(x = table_ORA_ClusterCompare_KEGG,
             paste(ora_dir, "table_ORA_ClusterCompare_KEGG.xlsx",sep = ""),
             col_names = TRUE, 
             format_headers = TRUE)
  
  # eliminar columnas innecesarias "nuevo bloque"
  table_ORA_ClusterCompare_KEGG_tmp <- table_ORA_ClusterCompare_KEGG[, c("ID",
                                                                         "Description",
                                                                         "ENTREZID",
                                                                         "Count",
                                                                         "GeneRatio",	
                                                                         "BgRatio",
                                                                         "Cluster",
                                                                         "pvalue",
                                                                         "p.adjust",
                                                                         "qvalue")]
  
  # separar columnas
  table_ORA_ClusterCompare_KEGG_long <- table_ORA_ClusterCompare_KEGG_tmp %>% 
    separate_rows(ENTREZID, sep = "/")
  
  # unir datos de entrada con datos de salida en un super data frame
  table_ORA_ClusterCompare_KEGG_final <- merge(table_ORA_ClusterCompare_KEGG_long,
                                               degs_sig_ora_kegg,
                                               by = "ENTREZID")
  
  # ordenar columnas por Description
  table_ORA_ClusterCompare_KEGG_final <- table_ORA_ClusterCompare_KEGG_final[order(table_ORA_ClusterCompare_KEGG_final$Description),]
  
  # cambiar orden de columnas "nuevo bloque"
  table_ORA_ClusterCompare_KEGG_final <- table_ORA_ClusterCompare_KEGG_final[,c("Gene_ID",
                                                                                "TAIR",
                                                                                "ENTREZID",
                                                                                "ENZYME",
                                                                                "ID",
                                                                                "Description", 
                                                                                "GeneRatio",
                                                                                "BgRatio",
                                                                                "Cluster",
                                                                                "pvalue",
                                                                                "p.adjust",
                                                                                "qvalue",
                                                                                "Blast1",
                                                                                "Gene_Symbol",
                                                                                "log2FC",
                                                                                "DEGs",
                                                                                "Pval",
                                                                                "Padj")]
  
  
  # guardar el objeto como tabla "nuevo bloque"
  write_xlsx(x = table_ORA_ClusterCompare_KEGG_final,
             paste(ora_dir, "table_ORA_ClusterCompare_KEGG_final.xlsx",sep = ""),
             col_names = TRUE, 
             format_headers = TRUE)
  
  ##################### dot
  
  # guardar imagen
  for(i in formatos) {
    
    # Usar tryCatch para atrapar errores y controlar su salida
    # que el error no detenga la funcion
    tryCatch({
      
      # crear y guardar los heatmaps
      match.fun(i)(paste(ora_dir, "ORA_KEGG_ClusterCompare_Dot_Ontology.", i, sep = ""),
                   res = resolucion,
                   width = ancho_kegg_cluster,
                   height = alto_kegg_cluster)
      
      # crear y guardar el heatmpat
      print(
        
        dotplot(compare_results_KEGG, 
                showCategory = categorias_kegg_cluster, 
                split = ".sign") +
          labs(title = "ORA Clusters: Enriched Pathways", x = "")
        
      )
      
      dev.off()
      
    }, error = function(e) {
      # imprimir el mensaje de error
      message("ocurrio un error durante la graficacion de GO-BP-Down: ", e$message)
    })
    
  }
  
  ################################################################################
  ################# Manual: separar por inducidos y reprimidos ###################
  ################################################################################
  
  # inducidos
  Activated_ora_kegg_up <- names(degs_sig_ora_kegg_vec_entrez[degs_sig_ora_kegg_vec_entrez > 0])
  
  # up
  oraKEGG_up <- enrichKEGG(gene = Activated_ora_kegg_up, 
                           universe = names(universe_ora_kegg_vec_entrez),
                           minGSSize = min_set_kegg,
                           maxGSSize = max_set_kegg,
                           organism = "ath", 
                           keyType = "ncbi-geneid",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05)
  
  # remover especie de descripcion
  oraKEGG_up@result$Description <- gsub(
    "- Arabidopsis thaliana \\(thale cress\\)$", "",
    oraKEGG_up@result$Description)
  
  # convertir a tabla
  table_ORA_KEGG_up <- as.data.frame(oraKEGG_up) %>% 
    mutate(GeneSet = "Activated")
  
  # cambiar nombre de columna de IDs de TAIR
  colnames(table_ORA_KEGG_up)[8] <- "ENTREZID"
  
  # guardar el objeto como tabla "nuevo bloque"
  write_xlsx(x = table_ORA_KEGG_up,
             paste(ora_dir, "table_ORA_KEGG_Up.xlsx",sep = ""),
             col_names = TRUE, 
             format_headers = TRUE)
  
  # eliminar columnas innecesarias
  table_ORA_KEGG_up_tmp <- table_ORA_KEGG_up[, c("ID",
                                                 "Description",
                                                 "ENTREZID",
                                                 "Count", 
                                                 "GeneRatio",
                                                 "BgRatio",
                                                 "GeneSet",
                                                 "pvalue",
                                                 "p.adjust",
                                                 "qvalue")]
  
  # separar columnas
  table_ORA_KEGG_up_long <- table_ORA_KEGG_up_tmp %>% 
    separate_rows(ENTREZID, sep = "/")
  
  # unir datos de entrada con datos de salida en un super data frame
  table_ORA_KEGG_up_final <- merge(table_ORA_KEGG_up_long, 
                                   degs_sig_ora_kegg,
                                   by = "ENTREZID")
  
  # remover duplicados
  table_ORA_KEGG_up_final <- table_ORA_KEGG_up_final[!duplicated(table_ORA_KEGG_up_final),]
  
  # ordenar columnas por Trinity ID
  table_ORA_KEGG_up_final <- table_ORA_KEGG_up_final[order(table_ORA_KEGG_up_final$Gene_ID),]
  
  # cambiar orden de columnas "nuevo bloque"
  table_ORA_KEGG_up_final <- table_ORA_KEGG_up_final[,c("Gene_ID",
                                                        "TAIR",
                                                        "ENTREZID",
                                                        "ENZYME",
                                                        "ID",
                                                        "Description",
                                                        "GeneRatio",
                                                        "BgRatio",
                                                        "pvalue",
                                                        "p.adjust",
                                                        "qvalue",
                                                        "GeneSet",
                                                        "Blast1",
                                                        "Gene_Symbol",
                                                        "log2FC",
                                                        "DEGs",
                                                        "Pval",
                                                        "Padj")]
  
  # guardar el objeto como tabla "nuevo bloque"
  write_xlsx(x = table_ORA_KEGG_up_final,
             paste(ora_dir, "table_ORA_KEGG_Up_final.xlsx",sep = ""),
             col_names = TRUE, 
             format_headers = TRUE)
  
  ################################################################################
  
  # reprimidos
  Suppressed_ora_kegg_down <- names(degs_sig_ora_kegg_vec_entrez[degs_sig_ora_kegg_vec_entrez < 0])
  
  # down
  oraKEGG_down <- enrichKEGG(gene = Suppressed_ora_kegg_down,
                             universe = names(universe_ora_kegg_vec_entrez),
                             minGSSize = min_set_kegg,
                             maxGSSize = max_set_kegg,
                             organism = "ath", 
                             keyType = "ncbi-geneid",
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.05)
  
  # remover especie de descripcion
  oraKEGG_down@result$Description <- gsub(
    "- Arabidopsis thaliana \\(thale cress\\)$", "",
    oraKEGG_down@result$Description)
  
  # convertir a tabla
  table_ORA_KEGG_down <- as.data.frame(oraKEGG_down) %>% 
    mutate(GeneSet = "Suppressed")
  
  # cambiar nombre de columna de IDs de TAIR
  colnames(table_ORA_KEGG_down)[8] <- "ENTREZID"
  
  # guardar el objeto como tabla
  write_xlsx(x = table_ORA_KEGG_down,
             paste(ora_dir, "table_ORA_KEGG_Down.xlsx",sep = ""),
             col_names = TRUE, 
             format_headers = TRUE)
  
  # eliminar columnas innecesarias "nuevo bloque"
  table_ORA_KEGG_down_tmp <- table_ORA_KEGG_down[, c("ID",
                                                     "Description",
                                                     "ENTREZID",
                                                     "Count",
                                                     "GeneRatio",
                                                     "BgRatio",
                                                     "GeneSet",
                                                     "pvalue",
                                                     "p.adjust",
                                                     "qvalue")]
  
  # separar columnas
  table_ORA_KEGG_down_long <- table_ORA_KEGG_down_tmp %>% 
    separate_rows(ENTREZID, sep = "/")
  
  # unir datos de entrada con datos de salida en un super data frame
  table_ORA_KEGG_down_final <- merge(table_ORA_KEGG_down_long, 
                                     degs_sig_ora_kegg,
                                     by = "ENTREZID")
  
  # remover duplicados
  table_ORA_KEGG_down_final <- table_ORA_KEGG_down_final[!duplicated(table_ORA_KEGG_down_final),]
  
  # ordenar columnas por Trinity ID
  table_ORA_KEGG_down_final <- table_ORA_KEGG_down_final[order(table_ORA_KEGG_down_final$Gene_ID),]
  
  # cambiar orden de columnas "nuevo bloque"
  table_ORA_KEGG_down_final <- table_ORA_KEGG_down_final[,c("Gene_ID",
                                                            "TAIR",
                                                            "ENTREZID",
                                                            "ENZYME",
                                                            "ID",
                                                            "Description", 
                                                            "GeneRatio",
                                                            "BgRatio",
                                                            "pvalue",
                                                            "p.adjust",
                                                            "qvalue",
                                                            "GeneSet",
                                                            "Blast1",
                                                            "Gene_Symbol",
                                                            "log2FC",
                                                            "DEGs",
                                                            "Pval",
                                                            "Padj")]
  
  
  # guardar el objeto como tabla "nuevo bloque"
  write_xlsx(x = table_ORA_KEGG_down_final,
             paste(ora_dir, "table_ORA_KEGG_Down_final.xlsx",sep = ""),
             col_names = TRUE, 
             format_headers = TRUE)
  
  ############################### graficar #######################################
  
  # ordenar por pajustado
  table_ORA_KEGG_down <- table_ORA_KEGG_down[order(table_ORA_KEGG_down$p.adjust),]
  table_ORA_KEGG_up <- table_ORA_KEGG_up[order((table_ORA_KEGG_up$p.adjust)),]
  
  # unir data frames
  resultados_combinados_ora_kegg <- rbind(table_ORA_KEGG_up, 
                                          table_ORA_KEGG_down)
  
  ##################### guardar archivos combinados ##############################
  
  # guardar el objeto como tabla
  write_xlsx(x = resultados_combinados_ora_kegg,
             paste(ora_dir, "table_ORA_KEGG_UpDown.xlsx", sep = ""),
             col_names = TRUE, 
             format_headers = TRUE)
  
  # unir data frames
  resultados_combinados_ora_kegg_final <- rbind(table_ORA_KEGG_up_final, 
                                                table_ORA_KEGG_down_final)
  
  # guardar el objeto como tabla "nuevo bloque"
  write_xlsx(x = resultados_combinados_ora_kegg_final,
             paste(ora_dir, "table_ORA_KEGG_UpDown_final.xlsx",sep = ""),
             col_names = TRUE, 
             format_headers = TRUE)
  
  ####################### graficar continuacion ##################################
  
  # crear col de gene ratio numerica
  for(i in 1:nrow(resultados_combinados_ora_kegg)){
    resultados_combinados_ora_kegg[i, "GeneRatioN"] <- eval(
      parse(text = resultados_combinados_ora_kegg[i, "GeneRatio"]))
  }
  
  # Filtrar mejores resultados para activados y suprimidos por p.adjust
  # up
  resultados_up_ora_kegg_tmp <- resultados_combinados_ora_kegg[resultados_combinados_ora_kegg$GeneSet == "Activated", ]
  resultados_up_ora_kegg <- resultados_up_ora_kegg_tmp[(1:(categorias_kegg/2)), ]
  resultados_up_ora_kegg <- na.omit(resultados_up_ora_kegg)
  
  # down
  resultados_down_ora_kegg_tmp <- resultados_combinados_ora_kegg[resultados_combinados_ora_kegg$GeneSet == "Suppressed", ]
  resultados_down_ora_kegg <- resultados_down_ora_kegg_tmp[(1:(categorias_kegg/2)),]
  resultados_down_ora_kegg <- na.omit(resultados_down_ora_kegg)
  
  # unir por filas
  resultados_top_ora_kegg <- rbind(resultados_up_ora_kegg,
                                   resultados_down_ora_kegg)
  
  ################################# graficar #####################################
  
  # guardar imagen
  for(i in formatos) {
    
    # Usar tryCatch para atrapar errores y controlar su salida
    # que el error no detenga la funcion
    tryCatch({
      
      # crear y guardar los heatmaps
      match.fun(i)(paste(ora_dir, "ORA_KEGG_Dot_UpDown.", i, sep = ""),
                   res = resolucion,
                   width = ancho_kegg,
                   height = alto_kegg)
      
      # crear y guardar el heatmpat
      print(
        
        ggplot(resultados_top_ora_kegg, 
               aes(x = GeneRatioN, 
                   y = Description)) +
          geom_point(aes(size = Count, 
                         color = p.adjust)) +
          facet_grid(.~ GeneSet) +
          scale_color_gradient(high = "blue", low = "red") +
          theme_dose() +
          theme(axis.text.y = element_text(size = 14),
                axis.text.x = element_text(size = 14),
                axis.title.x = element_text(size = 16),
                strip.text = element_text(size = 16),
                title = element_text(size = 18),
                legend.title = element_text(size = 16),
                legend.text = element_text(size = 14)) +
          labs(title = "ORA: Enriched Pathways",
               y = "",
               x = "Gene Ratio" ) +
          scale_x_continuous(n.breaks = 4)
        
      )
      
      dev.off()
      
    }, error = function(e) {
      # imprimir el mensaje de error
      message("ocurrio un error durante la graficacion de KEGG-Dot-UpDown: ", e$message)
    })
    
  }
  
  ######################## todos ##############################
  
  # guardar imagen
  for(i in formatos) {
    
    # Usar tryCatch para atrapar errores y controlar su salida
    # que el error no detenga la funcion
    tryCatch({
      
      # crear y guardar los heatmaps
      match.fun(i)(paste(ora_dir, "ORA_KEGG_Dot_All.", i, sep = ""),
                   res = resolucion,
                   width = ancho_kegg,
                   height = alto_kegg)
      
      # crear y guardar el heatmpat
      print(
        
        ggplot(resultados_combinados_ora_kegg, 
               aes(x = GeneRatioN, 
                   y = Description)) +
          geom_point(aes(size = Count, 
                         color = p.adjust)) +
          facet_grid(~GeneSet) +
          scale_color_gradient(high = "blue", low = "red") +
          theme_dose() +
          theme(axis.text.y = element_text(size = 14),
                axis.text.x = element_text(size = 14),
                axis.title.x = element_text(size = 16),
                strip.text = element_text(size = 16),
                title = element_text(size = 18),
                legend.title = element_text(size = 16),
                legend.text = element_text(size = 14)) +
          labs(title = "ORA: Enriched Pathways",
               y = "", 
               x = "Gene Ratio" ) +
          scale_x_continuous(n.breaks = 4)
        
      )
      
      dev.off()
      
    }, error = function(e) {
      # imprimir el mensaje de error
      message("ocurrio un error durante la graficacion de KEGG-Dot-All: ", e$message)
    })
    
  }
  
  ################################################################################
  ########################### ORA Encrichment Map ################################
  ################################################################################
  
  ######### GO #########
  ora_res_go <- pairwise_termsim(oraGO)
  
  # guardar imagen
  for(i in formatos) {
    
    # Usar tryCatch para atrapar errores y controlar su salida
    # que el error no detenga la funcion
    tryCatch({
      
      # crear y guardar los heatmaps
      match.fun(i)(paste(ora_dir, "ORA_GO_EMap.", i, sep = ""),
                   res = resolucion,
                   width = ancho_go_emap,
                   height = alto_go_emap)
      
      # crear y guardar el heatmpat
      print(
        
        emapplot(x = ora_res_go, 
                 showCategory = categorias_go) + 
          labs(title = "ORA: GO Enriched")
        
      )
      
      dev.off()
      
    }, error = function(e) {
      # imprimir el mensaje de error
      message("ocurrio un error durante la graficacion de ORA-Map: ", e$message)
    })
    
  }
  
  
  ######### GO UpDown #########
  
  # crear lista con nombres
  combinados_go <- list(Induced = oraGO_up, Repressed = oraGO_down)
  
  # unir listas 
  oraGO_updown <- merge_result(combinados_go)
  
  # create input file
  ora_res_go_updown <- pairwise_termsim(oraGO_updown)
  
  # guardar imagen
  for(i in formatos) {
    
    # Usar tryCatch para atrapar errores y controlar su salida
    # que el error no detenga la funcion
    tryCatch({
      
      # crear y guardar los heatmaps
      match.fun(i)(paste(ora_dir, "ORA_GO_UpDown_EMap.", i, sep = ""),
                   res = resolucion,
                   width = ancho_go_emap,
                   height = alto_go_emap)
      
      # crear y guardar el heatmpat
      print(
        
        emapplot(x = ora_res_go_updown, 
                 showCategory = categorias_go) + 
          labs(title = "ORA: GO Enriched", fill = "Gene Set")
        
      )
      
      dev.off()
      
    }, error = function(e) {
      # imprimir el mensaje de error
      message("ocurrio un error durante la graficacion de ORA-Map: ", e$message)
    })
    
  }
  
  ######### GO Up #########
  ora_res_go_up <- pairwise_termsim(oraGO_up)
  
  # guardar imagen
  for(i in formatos) {
    
    # Usar tryCatch para atrapar errores y controlar su salida
    # que el error no detenga la funcion
    tryCatch({
      
      # crear y guardar los heatmaps
      match.fun(i)(paste(ora_dir, "ORA_GO_Up_EMap.", i, sep = ""),
                   res = resolucion,
                   width = ancho_go_emap,
                   height = alto_go_emap)
      
      # crear y guardar el heatmpat
      print(
        
        emapplot(x = ora_res_go_up, 
                 showCategory = categorias_go) + 
          labs(title = "ORA: GO Enriched Activated")
        
      )
      
      dev.off()
      
    }, error = function(e) {
      # imprimir el mensaje de error
      message("ocurrio un error durante la graficacion de ORA-Map: ", e$message)
    })
    
  }
  
  ######### GO Down #########
  ora_res_go_down <- pairwise_termsim(oraGO_down)
  
  # guardar imagen
  for(i in formatos) {
    
    # Usar tryCatch para atrapar errores y controlar su salida
    # que el error no detenga la funcion
    tryCatch({
      
      # crear y guardar los heatmaps
      match.fun(i)(paste(ora_dir, "ORA_GO_Down_EMap.", i, sep = ""),
                   res = resolucion,
                   width = ancho_go_emap,
                   height = alto_go_emap)
      
      # crear y guardar el heatmpat
      print(
        
        emapplot(x = ora_res_go_down, 
                 showCategory = categorias_go) + 
          labs(title = "ORA: GO Enriched Suppressed")
        
      )
      
      dev.off()
      
    }, error = function(e) {
      # imprimir el mensaje de error
      message("ocurrio un error durante la graficacion de ORA-Map: ", e$message)
    })
    
  }
  
  ######### GO Clusters #########
  ora_res_go_clusters <- pairwise_termsim(compare_results_GO)
  
  # guardar imagen
  for(i in formatos) {
    
    # Usar tryCatch para atrapar errores y controlar su salida
    # que el error no detenga la funcion
    tryCatch({
      
      # crear y guardar los heatmaps
      match.fun(i)(paste(ora_dir, "ORA_GO_Clusters_EMap.", i, sep = ""),
                   res = resolucion,
                   width = 3000,
                   height = 3000)
      
      # crear y guardar el heatmpat
      print(
        
        emapplot(x = ora_res_go_clusters, 
                 showCategory = categorias_go_cluster) + 
          labs(title = "ORA Clusters: GO Enriched")
        
      )
      
      dev.off()
      
    }, error = function(e) {
      # imprimir el mensaje de error
      message("ocurrio un error durante la graficacion de ORA-Map: ", e$message)
    })
    
  }
  
  ############################### KEGG ###########################################
  
  # crear arcivo de entrada
  enrich_res_kegg <- pairwise_termsim(oraKEGG)
  
  # guardar imagen
  for(i in formatos) {
    
    # Usar tryCatch para atrapar errores y controlar su salida
    # que el error no detenga la funcion
    tryCatch({
      
      # crear y guardar los heatmaps
      match.fun(i)(paste(ora_dir, "ORA_KEGG_EMap.", i, sep = ""),
                   res = resolucion,
                   width = ancho_kegg,
                   height = alto_kegg)
      
      # crear y guardar el heatmpat
      print(
        
        emapplot(x = enrich_res_kegg, 
                 showCategory = categorias_kegg) + 
          labs(title = "ORA: KEGG Enriched")
        
      )
      
      dev.off()
      
    }, error = function(e) {
      # imprimir el mensaje de error
      message("ocurrio un error durante la graficacion de KEGG-EMap: ", e$message)
    })
    
  }
  
  ######### KEGG Clusters #########
  ora_res_kegg_clusters <- pairwise_termsim(compare_results_KEGG)
  
  # guardar imagen
  for(i in formatos) {
    
    # Usar tryCatch para atrapar errores y controlar su salida
    # que el error no detenga la funcion
    tryCatch({
      
      # crear y guardar los heatmaps
      match.fun(i)(paste(ora_dir, "ORA_KEGG_Clusters_EMap.", i, sep = ""),
                   res = resolucion,
                   width = ancho_kegg_emap,
                   height = alto_kegg_emap)
      
      # crear y guardar el heatmpat
      print(
        
        emapplot(x = ora_res_kegg_clusters, 
                 showCategory = categorias_kegg_cluster) + 
          labs(title = "ORA Clusters: Enriched Pathways")
        
      )
      
      dev.off()
      
    }, error = function(e) {
      # imprimir el mensaje de error
      message("ocurrio un error durante la graficacion de ORA-KEGG-Map: ", e$message)
    })
    
  }
  
  ######### KEGG UpDown #########
  
  # crear lista con nombres
  combinados_kegg <- list(Induced = oraKEGG_up, Repressed = oraKEGG_down)
  
  # unir listas 
  oraKEGG_updown <- merge_result(combinados_kegg)
  
  # create input file
  ora_res_kegg_updown <- pairwise_termsim(oraKEGG_updown)
  
  # guardar imagen
  for(i in formatos) {
    
    # Usar tryCatch para atrapar errores y controlar su salida
    # que el error no detenga la funcion
    tryCatch({
      
      # crear y guardar los heatmaps
      match.fun(i)(paste(ora_dir, "ORA_KEGG_UpDown_EMap.", i, sep = ""),
                   res = resolucion,
                   width = ancho_go_emap,
                   height = alto_go_emap)
      
      # crear y guardar el heatmpat
      print(
        
        emapplot(x = ora_res_kegg_updown, 
                 showCategory = categorias_go) + 
          labs(title = "ORA: Enriched Pathways", fill = "Gene Set")
        
      )
      
      dev.off()
      
    }, error = function(e) {
      # imprimir el mensaje de error
      message("ocurrio un error durante la graficacion de ORA-KEGG-Map: ", e$message)
    })
    
  }
  
  ######### KEGG Activated #########
  
  # Wrap the function call in tryCatch
  enrich_res_kegg_up <- tryCatch({
    pairwise_termsim(oraKEGG_up)
  }, error = function(e) {
    message("An error occurred while running pairwise_termsim: ", e$message)
    NULL  # Return NULL or a default value to allow the script to continue
  })
  
  # Check if enrich_res_kegg_up is not NULL before proceeding
  if (!is.null(enrich_res_kegg_up)) {
    print("pairwise_termsim executed successfully!")
    
    # guardar imagen
    for(i in formatos) {
      
      # crear y guardar los heatmaps
      match.fun(i)(paste(ora_dir, "ORA_KEGG_Up_EMap.", i, sep = ""),
                   res = resolucion,
                   width = ancho_kegg,
                   height = alto_kegg)
      
      # crear y guardar el heatmpat
      print(
        
        emapplot(x = enrich_res_kegg_up, 
                 showCategory = categorias_kegg) + 
          labs(title = "ORA: KEGG Enriched Activated")
        
      )
      
      dev.off()
      
    }
    
    
  } else {
    print("No enriched terms found, skipping further analysis.")
  }
  
  ######### KEGG Suppressed #########
  
  # Wrap the function call in tryCatch
  enrich_res_kegg_down <- tryCatch({
    pairwise_termsim(oraKEGG_down)
  }, error = function(e) {
    message("An error occurred while running pairwise_termsim: ", e$message)
    NULL  # Return NULL or a default value to allow the script to continue
  })
  
  # Check if enrich_res_kegg_up is not NULL before proceeding
  if (!is.null(enrich_res_kegg_down)) {
    print("pairwise_termsim executed successfully!")
    
    # guardar imagen
    for(i in formatos) {
      
      # crear y guardar los heatmaps
      match.fun(i)(paste(ora_dir, "ORA_KEGG_Down_EMap.", i, sep = ""),
                   res = resolucion,
                   width = ancho_kegg,
                   height = alto_kegg)
      
      # crear y guardar el heatmpat
      print(
        
        emapplot(x = enrich_res_kegg_down, 
                 showCategory = categorias_kegg) + 
          labs(title = "ORA: KEGG Enriched Suppressed")
        
      )
      
      dev.off()
      
    }
    
    
  } else {
    print("No enriched terms found, skipping further analysis.")
  }
  
  ################################################################################
  ###################### ORA KEGG individual pathways ############################
  ################################################################################
  ################################################################################
  
  ########################## rutas personalizadas ################################
  
  # # mostrar las rutas enriquecidas en KEGGgraph
  # cat("\nLas rutas enriquecidas inducidas son:", nrow(table_ORA_ClusterCompare_KEGG), "\n")
  # 
  # if(nrow(table_ORA_ClusterCompare_KEGG) == 0){
  #   cat("saltar analisis")
  # } else {
  # 
  #   for(i in 1:nrow(table_ORA_ClusterCompare_KEGG)){
  #     cat(paste(i,".- ", table_ORA_ClusterCompare_KEGG$ID[i],
  #               "\t",
  #               table_ORA_ClusterCompare_KEGG$Cluster[i],
  #               " = ",
  #               table_ORA_ClusterCompare_KEGG$Description[i], sep = ""), sep = "\n")
  #   }
  # 
  #   # pedir rutas inducidas a usuario
  #   KEGG_paths_cluster = readline(prompt = "Elige las rutas enriquecidas que deseas analizar: " )
  # 
  #   # volver vector
  #   KEGG_paths_cluster_vector <- unlist(strsplit(KEGG_paths_cluster, split = " |,|:|;"))
  # 
  #   # movernos a directorio destino
  #   setwd(dir = ora_dir)
  # 
  #   # obtener pathways
  #   for(ruta in KEGG_paths_cluster_vector){
  #     # Manejo de errores con tryCatch
  #     tryCatch({
  #       # Produce una gráfica de KEGG (PNG)
  #       pathview(gene.data = degs_sig_ora_kegg_vec,
  #                pathway.id= ruta,
  #                species = "ath",
  #                gene.idtype = gene.idtype.list[11])
  #     }, error = function(e) {
  #       cat("Error procesando la ruta:", ruta, ":", conditionMessage(e), "\n")
  #     })
  # 
  #   }
  # 
  #   # volver a home
  #   setwd("~")
  # }
  # 
  # # volver a home
  # setwd("~")
  # getwd() 
  
  ################################ graficar todas ################################
  
  # mostrar las rutas enriquecidas en KEGGgraph
  cat("\nLas rutas enriquecidas inducidas y reprimidas son:", nrow(table_ORA_ClusterCompare_KEGG), "\n")
  
  if(nrow(table_ORA_ClusterCompare_KEGG) == 0){
    cat("saltar analisis")
  } else {
    
    for(i in 1:nrow(table_ORA_ClusterCompare_KEGG)){
      cat(paste(i,".- ", table_ORA_ClusterCompare_KEGG$ID[i],
                "\t",
                table_ORA_ClusterCompare_KEGG$Cluster[i],
                " = ",
                table_ORA_ClusterCompare_KEGG$Description[i], sep = ""), sep = "\n")}
  }
  
  # obtener pathways
  for(ruta in table_ORA_ClusterCompare_KEGG$ID){
    # Manejo de errores con tryCatch
    tryCatch({
      # Produce una gráfica de KEGG (PNG)
      pathview(gene.data = degs_sig_ora_kegg_vec_tair,
               pathway.id= ruta,
               species = "ath",
               gene.idtype = "TAIR")
    }, error = function(e) {
      cat("Error procesando la ruta:", ruta, ":", conditionMessage(e), "\n")
    })
    
  }
  
  # volver a home
  setwd("~")
  
  ################################################################################
  ################################### GSEA #######################################
  ################################################################################
  
  # mostrar vector de universo de GSEA
  cat("\n\t\t============= 'COMENZANDO ANALISIS POR GSEA' =============\n")
  cat("\nvector de universo para GSEA\n")
  head(universe_gsea_vec)
  tail(universe_gsea_vec)
  
  # remover dupicados
  # universe_gsea_vec <- universe_gsea_vec[!duplicated(names(universe_gsea_vec))]
  
  ##################### dot
  
  # objeto de GSEA
  gseaGO <- gseGO(geneList = universe_gsea_vec, 
                  ont ="ALL", 
                  keyType = "TAIR", 
                  minGSSize = min_set_go, 
                  maxGSSize = max_set_go, 
                  pvalueCutoff = 0.05, 
                  pAdjustMethod = "BH",
                  verbose = TRUE, 
                  OrgDb = org.At.tair.db)
  
  # convertir a tabla
  table_GSEA_GO_tmp <- as.data.frame(gseaGO)
  
  # cambiar nombre de columna
  colnames(table_GSEA_GO_tmp)[12] <- "TAIR" 
  
  # guardar el objeto como tabla
  write_xlsx(x = table_GSEA_GO_tmp,
             path = paste(gsea_dir, "table_GSEA_GO.xlsx",sep = ""),
             col_names = TRUE,
             format_headers = TRUE)
  
  # eliminar columnas innecesarias
  table_GSEA_GO_tmp <- table_GSEA_GO_tmp[, c("ID",
                                             "Description",
                                             "ONTOLOGY",
                                             "TAIR",
                                             "setSize",
                                             "rank",
                                             "enrichmentScore",
                                             "NES",
                                             "leading_edge",
                                             "pvalue",
                                             "p.adjust",
                                             "qvalue")]
  
  # separar columnas
  table_GSEA_GO_long <- table_GSEA_GO_tmp %>% 
    separate_rows(TAIR, sep = "/")
  
  # unir datos de entrada con datos de salida en un super data frame
  table_GSEA_GO_final <- merge(table_GSEA_GO_long, 
                               universe_gsea_ordered,
                               by = "TAIR")
  
  # remover duplicados
  table_GSEA_GO_final <- table_GSEA_GO_final[!duplicated(table_GSEA_GO_final),]
  
  # ordenar columnas por Ranking
  table_GSEA_GO_final <- table_GSEA_GO_final[order(table_GSEA_GO_final$Ranking,decreasing = TRUE),]
  
  # cambiar orden de columnas
  table_GSEA_GO_final <- table_GSEA_GO_final[,c("Gene_ID",
                                                "TAIR",
                                                "ID",
                                                "Description",
                                                "ONTOLOGY",
                                                "rank",
                                                "enrichmentScore",
                                                "NES",
                                                "setSize",
                                                "leading_edge",
                                                "pvalue",
                                                "p.adjust",
                                                "qvalue",
                                                "Blast1",
                                                "Gene_Symbol",
                                                "log2FC",
                                                "DEGs",
                                                "Pval",
                                                "Padj",
                                                "Ranking")]
  
  # guardar el objeto como tabla
  write_xlsx(x = table_GSEA_GO_final,
             path = paste(gsea_dir, "table_GSEA_GO_final.xlsx", sep = ""),
             col_names = TRUE, 
             format_headers = TRUE)
  
  ##################### dot plot ONTOLOGY
  
  # guardar imagen
  for(i in formatos) {
    
    # Usar tryCatch para atrapar errores y controlar su salida
    # que el error no detenga la funcion
    tryCatch({
      
      # crear y guardar los heatmaps
      match.fun(i)(paste(gsea_dir, "GSEA_GO_ONTOLOGY1_Dot.", i, sep = ""),
                   res = resolucion,
                   width = ancho_go_gsea,
                   height = alto_go_gsea)
      
      # crear y guardar el heatmpat
      print(
        
        # dot plot of GO ALL
        dotplot(gseaGO,
                showCategory = categorias_go) +
          facet_grid(~ONTOLOGY) 
        
      )
      
      dev.off()
      
    }, error = function(e) {
      # imprimir el mensaje de error
      message("ocurrio un error durante la graficacion de GSEA-Dot: ", e$message)
    })
    
  }
  
  ##################### dot plot
  
  # guardar imagen
  for(i in formatos) {
    
    # Usar tryCatch para atrapar errores y controlar su salida
    # que el error no detenga la funcion
    tryCatch({
      
      # crear y guardar los heatmaps
      match.fun(i)(paste(gsea_dir, "GSEA_GO_ONTOLOGY2_Dot.", i, sep = ""),
                   res = resolucion,
                   width = ancho_go_gsea,
                   height = alto_go_gsea)
      
      # crear y guardar el heatmpat
      print(
        
        # dot plot of GO ALL
        dotplot(gseaGO,
                showCategory = categorias_go,
                split=".sign") +
          facet_grid(~ONTOLOGY) 
        
      )
      
      dev.off()
      
    }, error = function(e) {
      # imprimir el mensaje de error
      message("ocurrio un error durante la graficacion de GSEA-Dot: ", e$message)
    })
    
  }
  
  ################################################################################
  ############################# GSEA KEGG (Pathways) #############################
  ################################################################################
  
  "nuevo bloque"
  # TAIR to EC
  GSEA_TAIR_2_EC <-bitr(names(universe_gsea_vec),
                        fromType = "TAIR",
                        toType = "ENZYME",
                        OrgDb = org.At.tair.db)
  
  # obtener los Gene_IDs para la funcion enrichKEGG(), (podemos perder algunos genes aquí)
  GSEA_TAIR_2_ENTREZIDs <-bitr(names(universe_gsea_vec),
                               fromType = "TAIR",
                               toType = "ENTREZID",
                               OrgDb = org.At.tair.db)
  
  # unir ENTREZ, TAIR y ENZYME
  GSEA_TAIR_ENTREZ_ENZYME <- merge(GSEA_TAIR_2_ENTREZIDs, 
                                   GSEA_TAIR_2_EC, 
                                   by = "TAIR", 
                                   all.x = TRUE)
  
  # add ENTREZ IDs column
  universe_gsea_kegg <- merge(GSEA_TAIR_ENTREZ_ENZYME,
                              universe_gsea_ordered, 
                              by = "TAIR")
  
  # reordenar nombres de columnas
  universe_gsea_kegg <- universe_gsea_kegg[,c("Gene_ID",
                                              "TAIR",
                                              "ENTREZID",
                                              "ENZYME",
                                              "Ranking",
                                              "Blast1",
                                              "Gene_Symbol",
                                              "Common_2oo3",
                                              "log2FC",
                                              "Pval",
                                              "Padj",
                                              "DEGs",
                                              "Treatment")]
  
  # ordenar por ranking
  universe_gsea_kegg <- universe_gsea_kegg[order(universe_gsea_kegg$Ranking, decreasing = TRUE),]
  
  # crear vector de universo
  universe_gsea_kegg_vec_entrez <- universe_gsea_kegg$Ranking
  universe_gsea_kegg_vec_tair <- universe_gsea_kegg$Ranking
  
  # nombrar vector de universo
  names(universe_gsea_kegg_vec_entrez) <- universe_gsea_kegg$ENTREZID
  names(universe_gsea_kegg_vec_tair) <- universe_gsea_kegg$TAIR
  
  cat("\nvector GSEA KEGG\n")
  # verificar vector de universo
  cat("Primeros valores ENTREZID:", head(universe_gsea_kegg_vec_entrez))
  cat("Primeros valores TAIR:", head(universe_gsea_kegg_vec_tair))
  # verificar vector de universo
  cat("Ultimos valores ENTREZID:", tail(universe_gsea_kegg_vec_entrez))
  cat("Ultimos valores TAIR:", tail(universe_gsea_kegg_vec_tair))
  
  ################################################################################
  ################################ gseKEGG #######################################
  ################################################################################
  
  # Arabidopsis thaliana
  gseaKEGG <- gseKEGG(geneList = universe_gsea_kegg_vec_entrez,
                      organism = "ath",
                      pvalueCutoff = 0.05, 
                      minGSSize = min_set_kegg, 
                      maxGSSize = max_set_kegg,
                      pAdjustMethod = "BH",
                      keyType = "ncbi-geneid")
  #?gseKEGG()
  
  # remover especie de filas
  gseaKEGG@result$Description <- gsub(
    "- Arabidopsis thaliana \\(thale cress\\)$",
    "",
    gseaKEGG@result$Description)
  
  # convertir a tabla
  table_GSEA_KEGG <- as.data.frame(gseaKEGG)
  
  # cambiar nombre de columna de IDs de TAIR
  colnames(table_GSEA_KEGG)[11] <- "ENTREZID"
  
  # guardar el objeto como tabla
  write_xlsx(x = table_GSEA_KEGG,
             path = paste(gsea_dir, "/", "table_GSEA_KEGG.xlsx",sep = ""),
             col_names = TRUE,
             format_headers = TRUE)
  
  # separar columnas
  table_GSEA_KEGG_long <- table_GSEA_KEGG %>% 
    separate_rows(ENTREZID, sep = "/")
  
  # unir datos de entrada con datos de salida en un super data frame
  table_GSEA_KEGG_final <- merge(universe_gsea_kegg, 
                                 table_GSEA_KEGG_long, 
                                 by = "ENTREZID", 
                                 all.Y = TRUE)
  
  # remover duplicados
  table_GSEA_KEGG_final <- table_GSEA_KEGG_final[!duplicated(table_GSEA_KEGG_final),]
  
  # ordenar columnas por Trinity ID
  table_GSEA_KEGG_final <- table_GSEA_KEGG_final[order(table_GSEA_KEGG_final$Gene_ID),]
  
  # cambiar orden de columnas
  table_GSEA_KEGG_final <- table_GSEA_KEGG_final[,c("Gene_ID",
                                                    "TAIR",
                                                    "ENTREZID",
                                                    "ENZYME",
                                                    "ID",
                                                    "Description",
                                                    "enrichmentScore",
                                                    "Ranking",
                                                    "NES",
                                                    "pvalue",
                                                    "p.adjust",
                                                    "qvalue",
                                                    "Blast1",
                                                    "Gene_Symbol",
                                                    "Common_2oo3",
                                                    "log2FC",
                                                    "DEGs",
                                                    "Pval",
                                                    "Padj")]
  
  # guardar el objeto como tabla
  write_xlsx(x = table_GSEA_KEGG_final,
             path = paste(gsea_dir, "/", "table_GSEA_KEGG_final.xlsx",sep = ""),
             col_names = TRUE, 
             format_headers = TRUE)
  
  ##################### dot plot 
  
  # guardar imagen
  for(i in formatos) {
    
    # Usar tryCatch para atrapar errores y controlar su salida
    # que el error no detenga la funcion
    tryCatch({
      
      # crear y guardar los heatmaps
      match.fun(i)(paste(gsea_dir, "GSEA_KEGG_Dot.", i, sep = ""),
                   res = resolucion,
                   width = ancho_kegg,
                   height = alto_kegg)
      
      # crear y guardar el heatmpat
      print(
        
        # dot plot of GO ALL
        dotplot(gseaKEGG,
                showCategory = categorias_kegg,
                split=".sign") +
          facet_grid(. ~ .sign) 
        
      )
      
      dev.off()
      
    }, error = function(e) {
      # imprimir el mensaje de error
      message("ocurrio un error durante la graficacion de GSEA-Dot: ", e$message)
    })
    
  }
  
  ################################################################################
  ########################### GSEA Encrichment Map ###############################
  ################################################################################
  
  ######### GO #########
  gsea_res_go <- pairwise_termsim(gseaGO)
  
  # guardar imagen
  for(i in formatos) {
    
    # Usar tryCatch para atrapar errores y controlar su salida
    # que el error no detenga la funcion
    tryCatch({
      
      # crear y guardar los heatmaps
      match.fun(i)(paste(gsea_dir, "GSEA_GO_EMap.", i, sep = ""),
                   res = resolucion,
                   width = ancho_go_emap,
                   height = alto_go_emap)
      
      # crear y guardar el heatmpat
      print(
        
        emapplot(x = gsea_res_go,
                 showCategory = categorias_go,
                 split = "NES") + 
          labs(title = "GSEA: GO Enriched")
        
      )
      
      dev.off()
      
    }, error = function(e) {
      # imprimir el mensaje de error
      message("ocurrio un error durante la graficacion de GSEA-Map: ", e$message)
    })
    
  }
  
  ############################### KEGG ###########################################
  
  enrich_res_kegg <- pairwise_termsim(gseaKEGG)
  
  # guardar imagen
  for(i in formatos) {
    
    # Usar tryCatch para atrapar errores y controlar su salida
    # que el error no detenga la funcion
    tryCatch({
      
      # crear y guardar los heatmaps
      match.fun(i)(paste(gsea_dir, "GSEA_KEGG_EMap.", i, sep = ""),
                   res = resolucion,
                   width = ancho_kegg,
                   height = alto_kegg)
      
      # crear y guardar el heatmpat
      print(
        
        emapplot(x = enrich_res_kegg, 
                 showCategory = categorias_kegg,
                 split = "NES") + 
          labs(title = "GSEA: KEGG Enriched")
        
      )
      
      dev.off()
      
    }, error = function(e) {
      # imprimir el mensaje de error
      message("ocurrio un error durante la graficacion de GSEA-EMap: ", e$message)
    })
    
  }
  
  
  ################################################################################
  ###################### GSEA KEGG individual pathways ############################
  ################################################################################
  
  # mostrar las rutas enriquecidas en KEGGgraph
  cat("\nLas rutas enriquecidas totales son:", nrow(table_GSEA_KEGG),"\n")
  for(i in 1:nrow(table_GSEA_KEGG)){
    cat(paste(i,".- ", table_GSEA_KEGG$ID[i], "\t", "= ",
              table_GSEA_KEGG$Description[i], sep = ""), sep = "\n")
  }
  
  # # pedir rutas a ususario
  # KEGG_paths = readline(prompt = "Elige las rutas que deseas analizar: " )
  # 
  # # volver vector
  # KEGG_paths_vector <- unlist(strsplit(KEGG_paths, split = " |,|:|;"))
  # 
  # # movernos a directorio destino
  # setwd(dir = gsea_dir)
  
  ################################################################################
  # Todos juntos
  for(ruta in table_GSEA_KEGG$ID){
    # Manejo de errores con tryCatch
    tryCatch({
      # Produce una gráfica de KEGG (PNG)
      pathview(gene.data = universe_gsea_kegg_vec_tair,
               pathway.id= ruta,
               species = "ath",
               gene.idtype = "TAIR")
    }, error = function(e) {
      cat("Error procesando la ruta:", ruta, ":", conditionMessage(e), "\n")
    })
    
  }
  
  # volver a home
  setwd("~")
  
}