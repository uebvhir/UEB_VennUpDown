###Plot Venn diagramm from 2 to 5 comparisons
###Function based in previous work of Ferran Briansó and Miriam Mota ("https://raw.githubusercontent.com/rgonzalo/GitBackupScripts/master/VennDiagrams/CreateVennEuler.R") 
###Parameters available to tune:
#topTabs = lista de los data.frame con las toptables
#compNames = lista con los nombres de cada comparación
#label = nombre que quieres que salga en el archivo
#colFeat = Columna por la que busca genes comunes
#colPVal = columna para hacer la selección (normalmente Pvalue)
#pval =  p valor por el que se quiere hacer la selección
#pltR = enseñar la imagen en R,
#pltPdf = hacer el pdf,
#eul = dibujar el euler plot,
#venn = TRUE,
#csv = crear l'arxiu csv amb els gens comuns+
#colors = si volem colors adhoc (sino fa rainbow)
#trans = transparencia dels colors
#cex1 = tamaño de las etiquetas
#cex2 = tamaño de los números
#rotation = rotación del venn en global
#position = vector con las posiciones en grados de cada etiqueta de comparación
#FC = valor de logFC para seleccionar genes

#include = dirección en que quieren seleccionarse los genes (up / down / both)


createVennEuler <- function(topTabs, compNames, label = "selected", colFeat = "X", 
                            colPVal = "P.Value", pval = 0.05, pltR = TRUE, 
                            pltPdf = TRUE, pltPng=FALSE, venn = TRUE, eul = TRUE, csv = TRUE, 
                            colors = rainbow(length(compNames)), trans = 0.5, 
                            cex1 = 1.5, rotation, position, cex2 = 1.2, FC, include = ""){
    
    ## Initializing lists
    list_genes_sel <- list()
    ## Reading input data
    for (i in 1:length(topTabs)) {
        colpval <- which(names(topTabs[[i]]) == colPVal)
        colFC <- which(names(topTabs[[i]]) == "logFC")
        colFeature <- which(names(topTabs[[i]]) == colFeat)
        if (include==""){list_genes_sel[[i]] <- as.character(topTabs[[i]][, colFeat][topTabs[[i]][, colpval] < pval & abs(topTabs[[i]][, colFC]) > FC])}
        if (include=="Up"){list_genes_sel[[i]] <- as.character(topTabs[[i]][, colFeat][topTabs[[i]][, colpval] < pval & topTabs[[i]][, colFC] > FC])}
        if (include=="Down"){list_genes_sel[[i]] <- as.character(topTabs[[i]][, colFeat][topTabs[[i]][, colpval] < pval & topTabs[[i]][, colFC] < -FC])}
    }
    
    ## Creating Venn Diagram
    if (venn) {
        if (include=="") {titulo <- paste0("Venn diagram for comparison: \n", label, "\n(" , colPVal, " < ", pval," & abs(logFC) > ", FC, ") \n")}
        if (include=="Up") {titulo <- paste0("Venn diagram for comparison: \n", label, "\nUp-regulated genes (" , colPVal, " < ", pval," & logFC > ", FC, ") \n")}
        if (include=="Down") {titulo <- paste0("Venn diagram for comparison: \n", label, "\nDown-regulated genes (" , colPVal, " < ", pval," & logFC < -", FC, ") \n")}
        venn.plot <- venn.diagram(list_genes_sel,
                                  category.names = compNames,
                                  fill = colors, 
                                  alpha = trans,
                                  resolution = 800,
                                  cat.cex = cex1,
                                  cex = cex2,
                                  main.cex=cex2,
                                  main = titulo,
                                  filename = NULL,
                                  rotation.degree = rotation,
                                  cat.pos = position)
        if (pltPdf) {
            pdf(paste0("VennDiagram", include, ".",label, ".", colPVal, pval, ".logFC",FC,".pdf"))
            grid.draw(venn.plot)
            dev.off()
        }
        if (pltPng) {
            png(paste0("VennDiagram", include, ".",label, ".", colPVal, pval, ".logFC",FC,".png"))
            grid.draw(venn.plot)
            dev.off()
        }
        if (pltR) {grid.draw(venn.plot)}
    }
    
    ## Creating Euler Diagram
    if (eul) {
        set <- NULL
        for (i in 1:length(compNames)) {
            set <- c(set, rep(compNames[i],length(list_genes_sel[[i]])))
        }
        v <- venneuler(data.frame(elements = c(unlist(list_genes_sel)),
                                  sets = set))
        if (pltPdf) {
            pdf(paste0("EulerDiagram", include, ".", label, ".", colPVal, pval, ".pdf"))
            plot(v, main = paste0("Euler diagram (" , colPVal, " < ", pval,")"))
            dev.off()
        }
        if (pltR) {plot(v, main = paste0("Euler diagram (" , colPVal, " < ", pval,")"))}
    }
    
    ## Obtaining lists of combined elements and computing shared elements
    
    names(list_genes_sel) <- paste0(compNames, "_")
    combs <-  unlist(lapply(1:length(list_genes_sel),
                            function(j) combn(names(list_genes_sel), j, simplify = FALSE)),
                     recursive = FALSE)
    names(combs) <- sapply(combs, function(i) paste0(i, collapse = ""))
    #str(combs)
    
    elements <- lapply(combs, function(i) Setdiff(list_genes_sel[i], list_genes_sel[setdiff(names(list_genes_sel), i)]))
    n.elements <- sapply(elements, length)
    list_res <- list(elements = elements, n.elements = n.elements)
    
    seq.max <- seq_len(max(n.elements))
    mat <- sapply(elements, "[", i = seq.max)
    mat[is.na(mat)] <- ""
    sharedElements <- rbind(t(data.frame(n.elements)),data.frame(mat))
    
    ## Writing table of shared elements as csv
    if (csv) {
        WriteXLS(sharedElements, ExcelFileName = paste0("sharedElements", include, ".", label, ".", colPVal, pval,".logFC",FC,".xls"), 
                 row.names = FALSE)
    }
    
    ## Returning table as data.frame
    return(sharedElements)
}
#####################################################################


###################################################################
###### Internal functions used to extract the lists of shared elements 
###### both x and y are, in all cases, lists of character strings
###################################################################
Intersect <- function(x) {
    if (length(x) == 1) {
        unlist(x)
    } else if (length(x) == 2) {
        intersect(x[[1]], x[[2]])
    } else if (length(x) > 2) {
        intersect(x[[1]], Intersect(x[-1]))
    }
}

Union <- function(x) {
    if (length(x) == 1) {
        unlist(x)
    } else if (length(x) == 2) {
        union(x[[1]], x[[2]])
    } else if (length(x) > 2) {
        union(x[[1]], Union(x[-1]))
    }
}

Setdiff <- function(x, y) {
    xx <- Intersect(x)
    yy <- Union(y)
    setdiff(xx, yy)
}
