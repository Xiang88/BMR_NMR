options(stringsAsFactors = FALSE)
library(dplyr)
library(plotly)

setwd("/home/xiang/DNAm/github_dnamage/BMR")
base_name = c(
	"/home/xiang/DNAm/N13.2018-9246RodentsSkinGorbunova",
	"/home/xiang/DNAm/N24.2019-9058MammalsLiverVeraGorbunova"
	) [index_tissue]
base_number = c("13","24") [index_tissue]

info_name = "/home/xiang/DNAm/annotation_info"



cpg_set_for_average = c(
	"intersect( top_n_conserved , union( mm10_island_cpg , hg38_island_cpg ))" ,
	"setdiff( top_n_conserved , union( mm10_island_cpg , hg38_island_cpg ) )" ,

##old
	"setdiff( all_cpg , mappable_cpg )" , 
	"intersect(all_cpg , mappable_cpg )" , 
	"top1000conserved" , 
	"mm10_island_cpg" , 
	"hg38_island_cpg" ) [index_set]

keep_species = c("Mus musculus" , "Rattus norvegicus" , "Spalax galili" , "Nannospalax ehrenbergi")
hl_species = c('Nannospalax ehrenbergi' , 'Spalax galili')

cat("base_name: ",base_name,"\n")
cat("base_number: ",base_number,"\n")
cat("cpg_set_for_average: ",cpg_set_for_average,"\n")
cat("top_n: ",top_n,"\n")

######################################################
######################################################

if(!exists("uenv")) uenv = new.env()

load(file.path(base_name , "NormalizedData" , "all_probes_sesame_normalized.Rdata") , uenv)
uenv$ss = read.csv(file.path(base_name,sprintf("SampleSheetMinimal%s.csv", base_number)))

	
if(!exists("island",envir=uenv))
	uenv$island = read.csv(file.path(info_name,"probes_CGislands_mm10_hg38.csv"))

if(!exists("coord",envir=uenv))
	uenv$coord = read.csv(file.path(info_name,"HorvathMammalChip_SpeciesGenomeCoordsUnique_v1.0.csv"))

if(!exists("mappability",envir=uenv))
	uenv$mappability = read.csv(file.path(info_name,"MammalChip_GenomeCoords_021119.tsv"),sep='\t')

## Based on old coord file
## mappable_cpg = (uenv$mappability %>% filter( !is.na(.data[['Northern_Israeli_blind_mole_rat']]) ))$CpG_ID
## mappable_cpg = substr(mappable_cpg,1,10)

mappable_cpg = (uenv$coord %>% filter(.data$NannospalaxGalili!="" ))$probeID
hg38_island_cpg = (uenv$island %>% filter( hg38_chr_island!="" ))$probeID
mm10_island_cpg = (uenv$island %>% filter( mm10_chr_island!="" ))$probeID

top_n_conserved <-
	(function(){
		# Grab top1000 most conserved CpGs 
		# minus 1 because counting $CpG_ID 
		#order and take top 1000 most common

		alist = uenv$coord
		alist$num_in_species = rowSums(!(alist == "")) - 1 
		trunlist = alist[order(-(alist$num_in_species)),,drop=F] 
		trunlist[1:top_n,,drop=F]
	})()
top_n_conserved = top_n_conserved$probeID



dnam_prepare <- function(ss , beta) {
	## Merge beta values with Sample Sheet
	## Assume the columns in beta are CGid , <Basename 1> , <Basename 2> , ...
	## Assume a column Basename exists in ss

	d = beta
	e = ss
	sid = setdiff( names(d) , "CGid")
	cpgid = d$CGid
	df = data.frame( t(d[,sid] ) , row.names=sid  )
	names(df) = cpgid
	m = merge(e , df , by.x="Basename", by.y="row.names", all.x=TRUE)
	m
}

d = dnam_prepare(beta = uenv$normalized_betas_sesame , ss = uenv$ss)
d = d[ d$SpeciesLatinName %in% keep_species,]
d$grp = ifelse(d$SpeciesLatinName %in% hl_species , "Grp" , " Other")
all_cpg = names(d) [ grep("^cg\\d|ch\\.\\d",names(d)) ]

data_col = NULL
data_col = eval(parse(text=cpg_set_for_average))

d$gross_mean = rowMeans( d[, data_col] )


p = NULL
p = plot_ly(d, y = ~gross_mean, x = ~SpeciesLatinName ,
	color = ~SpeciesLatinName, 
	size=~c(10,30)[factor(grp)], 
	type="box" ,
	width = 1000, height = 800) %>%
    layout( yaxis = list(range=c(0,1)) )

print(p)

