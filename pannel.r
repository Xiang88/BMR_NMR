top_n = 5000

index_tissue=1
index_set=1
source("run.r")
p11 = p %>% layout(yaxis=list(title="Mean Methylation CpG Island"),
	annotations = list(x = 1 , y =.95, text = "A", showarrow=F)
		)

index_tissue=1
index_set=2
source("run.r")
p12 = p %>% layout(yaxis=list(title="Mean Methylation Non-CpG Island"),
	annotations = list(x = 1 , y =.95, text = "B", showarrow=F)
	)

index_tissue=2
index_set=1
source("run.r")
p21 = p %>% layout(yaxis=list(title="Mean Methylation CpG Island"),
	annotations = list(x = 1 , y =.95, text = "C", showarrow=F)
	)

index_tissue=2
index_set=2
source("run.r")
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



