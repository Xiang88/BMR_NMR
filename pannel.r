top_n = 5000

index_tissue=1
index_set=1
source("run.r")
d11 = d
p11 = p %>% layout(yaxis=list(title="Mean Methylation CpG Island"),
	annotations = list(x = 1 , y =.95, text = "A", showarrow=F)
		)

index_tissue=1
index_set=2
source("run.r")
d12 = d
p12 = p %>% layout(yaxis=list(title="Mean Methylation Non-CpG Island"),
	annotations = list(x = 1 , y =.95, text = "B", showarrow=F)
	)

index_tissue=2
index_set=1
source("run.r")
d21 = d
p21 = p %>% layout(yaxis=list(title="Mean Methylation CpG Island"),
	annotations = list(x = 1 , y =.95, text = "C", showarrow=F)
	)

index_tissue=2
index_set=2
source("run.r")
d22 = d
p22 = p %>% layout(yaxis=list(title="Mean Methylation Non-CpG Island"),
	annotations = list(x = 1 , y =.95, text = "D", showarrow=F)
)

p = NULL
p = subplot(p11 ,
	p12 ,
	p21 ,
	p22 ,
	nrows=2 , titleY=T , margin=.05) %>%
    layout(showlegend = FALSE )

print(p)





##############################################
##  Bar plot combining two rows for tissue  ##
##############################################
d11$island = 1
d21$island = 1
d12$island = 0
d22$island = 0
d_total=rbind(d11,d12,d21,d22)
summary( aov(gross_mean~SpeciesLatinName+Tissue+island,data=d_total) )

s_1 = rbind(d11,d21) %>% group_by(Tissue,SpeciesLatinName) %>% summarise(m=mean(gross_mean),sd=sd(gross_mean))
s_2 = rbind(d12,d22) %>% group_by(Tissue,SpeciesLatinName) %>% summarise(m=mean(gross_mean),sd=sd(gross_mean))

p_1 = plot_ly(s_1, x=~SpeciesLatinName , y=~m , type='bar' , color=~Tissue , showlegend=T ,
		error_y=list(array=~sd,color = '#000000')
	) %>%
	layout(yaxis=list(title="Mean Methylation CpG Island",
		range=c(0,1)),barmode='group' ,
		annotations = list(x = 1 , y =.98, text = "A", showarrow=F)
	)

p_2 = plot_ly(s_2, x=~SpeciesLatinName , y=~m , type='bar' , color=~Tissue , showlegend=F ,
		error_y=list(array=~sd,color = '#000000')
	) %>%
	layout(yaxis=list(title="Mean Methylation Non-CpG Island",
		range=c(0,1)) , barmode='group' ,
		annotations = list(x = 1 , y =.98, text = "B", showarrow=F)
	)

p_ = subplot(p_1 , p_2 , titleY=T , margin=.08 ) %>%
	layout( legend = list(y = 0.85) )

p__ = browsable(
	div(style="display: block; width: 700px;",
		div(style="display: block;" , p_),
		div(style="display: block; padding: 0 14% 0 10%; font-family: sans-serif;",
			p(strong("Fig. ") , "Average methylation among three species (Mus musculus, Rattus norvegicus, and Nannospalax ehrenbergi), two tissues (liver and skin), and within CpG island (A) and out side CpG island (B). CpG's are limited to the top 5000 converved sites in the genome of mammals. The results qualitatively are the same no matter which top-number we choose. A 3-way ANOVA gives p=0.013 for species, p=1e-9 for tissue, and p=2e-16 for CpG island.")
			)
	))

print(p__)



