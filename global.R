# root.path= "C:/Users/Siamak/" 
# shinyRappToyDataset_SiamakPlay.csv
# 800-H1-H20-RNA-Seq-SingleCell-Retina-OMRF-03-29-19_FPKM_v2_SiamakPlay.csv
library(Seurat)
singleCellDataAll_H <- read.csv(file= "800-H1-H20-RNA-Seq-SingleCell-Retina-OMRF-03-29-19_FPKM_v2_SiamakPlay.csv", head=TRUE, sep= ",", stringsAsFactors= FALSE)  # Which Plate? data has duplicated genes  

singleCellDataGenes_H= singleCellDataAll_H[, c(8, 11:ncol(singleCellDataAll_H))]                # all THY1+
rownames(singleCellDataGenes_H)= singleCellDataAll_H$Ensembl.ID

# Removing duplicates
selectHighestExpressedDuplicatedGenesFun <- function(singleCellDataGenes){
  # ------------- Smart and quick way to find genes that are more expressed from duplicated genes ------------------------ #
  # creating a temporary label that is a combination of gene and the amount of expression
  # *Note: sep="-" is very important to avoid gene name confusion
  singleCellDataGenes$tmp.label= paste(singleCellDataGenes$Gene.Symbol, rowSums( singleCellDataGenes[, c(2:ncol(singleCellDataGenes))] > 0 )+1000, sep="-")
  # order data based on this new label
  single.tmp= singleCellDataGenes[order(singleCellDataGenes$tmp.label), ]
  # pick the last duplicated label (belongs to highest expressed gene)
  single.tmp.new= single.tmp[!duplicated(single.tmp$Gene.Symbol, fromLast= TRUE),]
  
  #----------- Test Routine: to make sure cell with greater number of genes selected, 1/3/2018
  t1= single.tmp[!duplicated(single.tmp$Gene.Symbol, fromLast= TRUE),]
  t2= single.tmp[!duplicated(single.tmp$Gene.Symbol, fromLast= FALSE),]
  print("Two test routines are perfomring:")
  r1= sum( rowSums( t1[, c(2:ncol(singleCellDataGenes))] > 0 ) < rowSums( t2[, c(2:ncol(singleCellDataGenes))] > 0 ) )  # sum should be zero
  r2= sum(duplicated(t1$Gene.Symbol))                                       # final check after removing duplicates
  print(paste("Number of genes with wrong lower expresseion is: ", toString(r1), sep= ""))
  print(paste("Number of duplicated genes is: ", toString(r2), sep=""))
  # ---------------------------------------------------------------------------------------------------------------------- #
  return ( single.tmp.new[,1:(ncol(singleCellDataGenes)-1)] )
}

singleCellDataGenes_H= selectHighestExpressedDuplicatedGenesFun(singleCellDataGenes_H)  # [25,394, 801]

rownames(singleCellDataGenes_H)= singleCellDataGenes_H[,1]

#____________________________________________ Meta data
rgcMasterList <- read.csv(file= "CellTypes_RGC_Master_08Dec2018.csv", head=TRUE, sep= "," , stringsAsFactors = FALSE)      # RGC master list provided by ROB
rgc1000List   <- read.csv(file= "RobTop1001.csv", head=TRUE, sep= "," , stringsAsFactors = FALSE)     # Rob's top 1000 genes, 4/11/2018
rgc1000List.genes= rgc1000List$Symbol

#______________________________________________________ Number of genes per cell
singleCellDataGenes= singleCellDataGenes_H # singleCellDataGenes_ABCDEFGIJKLM

# plot this on Shiny
n.cells= ncol(singleCellDataGenes)-1
n.genes.per.cell= colSums( singleCellDataGenes[ , 2: (n.cells+1)] > 0 ) # 
# hist(n.genes.per.cell, breaks= 100, main="", xlab="Number of genes expressed", ylab= "Frequency", col= c("gray"), cex.axis= 1.5, cex.lab= 1.6 )

# ___________________________________________________ RGC markers
# Selecting cells (e.g., those with at least 2 of the 6 canonical RGC markers are present) and genes
# "Thy1", "Rbpms", "Rpbms2", "Jam2", "G3bp1", 
# Provided by Rob: these are rgc amrkers, if at least two of them are expressed in a cell, more likely rgc 
canonical.rgc.markers.1= c("Thy1", "Rbpms", "Rbpms2", "Jam2", "G3bp1", "Ywhaz")
canonical.rgc.markers.2= c("Nrcam", "Rtn1", "Tecr", "Spock2", "Ica1l", "Ywhah", "Vsnl1", "Ica1", "Chrnb3", "Nefl", "Stmn2", 
                           "Sncg", "Uchl1", "Tubb2a", "Fxyd7", "Lynx1", "Chrna6", "Tubb2b", "Atoh7", "Tubb4a", "Nefm", "Nrn1", "Pou4f2")
reference.rgc.markers.3= c("Txn1", "Ywhah", "Apbb2", "Epb4.1l3", "Thy1", "Gam2", "Nell1", "G3bp1", "Fgf12", "Rbpms", "Arhgap44", "Foxp2", "Tmod2", "Smad2", 
                           "Nefl", "Atl1", "Fstl4", "Bend6", "Cbx6", "Rbpms2")  # RGC master based on Macosco and high RGC retina enrichment and they all have expression above in whole retina > 8
reference.rgc.markers.3= reference.rgc.markers.3[ which( reference.rgc.markers.3 %in% singleCellDataGenes$Gene.Symbol ) ]  # only those genes that exist in our dataset
# pan RGC markers based on Rheaume Nature Communicatins article 2018
reference.rgc.markers.4= c("Tubb3", "Rbpms", "Sox4", "Sox11", "Sox12", "Pou4f1", "Pou4f2", "Pou4f3", "Isl1", "Isl2")
subtypes.rgc.markers.4=  c("Ebf1", "Ebf2", "Ebf3", "Ebf4", "Satb1", "Satb2", "Barhl2", "Neurod2", "Foxp1", "Foxp2", "Cartpt", "Col25a1", "Cdh6", "Fstl4", "Npy", 
                           "Prdm16", "Sdk1", "Sdk2", "Tgfb1", "Jam2", "Cntnap4", "Mmp17", "Pvalb", "Pde1a", "Eomes", "Igf1", "Opn4", "Spp1", "Gna14", "Trhr", "Pcp2", "Dlx2", "Htr2a", "Htr2c")
reference.rgc.markers.4= reference.rgc.markers.4[ which( reference.rgc.markers.4 %in% singleCellDataGenes$Gene.Symbol ) ]  # only those genes that exist in our dataset

canonical.rgc.markers.all= c(canonical.rgc.markers.1, canonical.rgc.markers.2, reference.rgc.markers.3)

n.canonical.rgc.markers.per.cell= colSums( singleCellDataGenes[which(singleCellDataGenes$Gene.Symbol %in% reference.rgc.markers.3), 2: (n.cells+1)] > 0 )

# Barplot of reference RGC markers
can.cell.counts= table(n.canonical.rgc.markers.per.cell)

# ____________________________________ log transform
data.log= log(singleCellDataGenes[ , 2: (n.cells+1)]+1)   # selected genes, selected cells 

# ____________________________________ making Seurat object
sc.object  <- CreateSeuratObject(raw.data= data.log, min.cells= 0, min.genes= 0, project= "Single_cell_data" )  # 2000 genes, 900 genes in Macosko, genes expressed in fewer than 0.05% of cells are exluded


# to find a gene name simply
# names(ave.genes.expressed)[grep(pattern= "Cry", names(ave.genes.expressed), value= FALSE)]

#_______________________________________ adding meta-data and normalization
mito.genes <- grep(pattern = "^mt-", x = rownames(x = sc.object@data), value = TRUE)
percent.mito <- Matrix::colSums(sc.object@raw.data[mito.genes, ])/Matrix::colSums(sc.object@raw.data)

# adding row indent to cells
cells.selected= which( names(x = data.log)  %in% names(x = sc.object@data) )                      # which cells are selected for downstream analysis

row.ident= factor( c( rep( paste("R", 1:40, sep= ""), times= 20*1) ))[cells.selected] # , rep( paste("R", 81:120, sep= ""), times= 10*1) ) )
row.ident.values= c( rep( 1:40, times= 20*1 ) )  #, rep( 81:120, times= 10*1 ) )  
sc.object <- AddMetaData(object= sc.object, metadata= row.ident, col.name= "row.ident")      # adding row identity metaData
sc.object@meta.data$row.ident= row.ident.values[cells.selected]  # or use candidate RGC cells [idx.rgc.A] 

sc.object <- NormalizeData(object= sc.object, normalization.method= "LogNormalize", scale.factor= 10000)  
dim(sc.object@raw.data)  # [14,366 * 794]

# Scaled expression
ave.expressed.genes= log ( x= rowMeans( x= exp(as.matrix(sc.object@data)) -1) + 1 )
ave.expressed.genes["Thy1"]
# dis.dispersed.genes= log ( x= apply(exp(as.matrix(sc.object@data)) -1, 1, var) +1 )
# dis.dispersed.genes["Thy1"]

# Non-scaled expression
ave.expressed.genes.2= log ( x= rowMeans( x= exp(sc.object@raw.data) -1) + 1 )
ave.expressed.genes.2["Thy1"]

# back to original scale from log [cell expr/sum(cell expr)*10.000 +1]
expressed.genes.3= log (    ( ( exp(as.matrix(sc.object@data)) -1)/10000.0 )*matrix( rep(colSums(exp(sc.object@raw.data) -1), each= nrow(sc.object@raw.data)), nrow=nrow(sc.object@raw.data) )  + 1 )
names(expressed.genes.3)= names(ave.expressed.genes)
expressed.genes.3["Thy1"]


# ExpMean= log(x = mean(x = exp(x = x) - 1) + 1)
# LogVMR= log(x = var(x = exp(x = x) - 1) / mean(x = exp(x = x) - 1)), VMR: variance-to-mean ratio

# variable genes must be identified beforehand 
# VariableGenePlot(sc.object, do.text= TRUE, cex.use= 0.5, cex.text.use= 0.5, do.spike= FALSE, pch.use= 19, col.use= "gray",  spike.col.use= "red", plot.both= FALSE, do.contour= FALSE,
#                   contour.lwd= 3, contour.col= "white", contour.lty= 2, x.low.cutoff= -Inf, x.high.cutoff= Inf, y.cutoff= -Inf)
# pass.cutoff.1= which( rownames(sc.object.1@hvg.info) %in% c(reference.rgc.markers.4, subtypes.rgc.markers.4) ) # Rheaume pan and subtype RGCs
# pass.cutoff.2= which( rownames(sc.object.1@hvg.info) %in% c(reference.rgc.markers.3) )                         # Rob's RGC genes
# text( sc.object.1@hvg.info$gene.mean[pass.cutoff.1], sc.object.1@hvg.info$gene.dispersion.scaled[pass.cutoff.1], rownames(sc.object.1@hvg.info)[pass.cutoff.1], col= "red", cex= 0.8)
# text( sc.object.1@hvg.info$gene.mean[pass.cutoff.2], sc.object.1@hvg.info$gene.dispersion.scaled[pass.cutoff.2], rownames(sc.object.1@hvg.info)[pass.cutoff.2], col= "blue", cex= 0.8)

