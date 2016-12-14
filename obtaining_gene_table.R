setwd("//SCOMP0856/veen119$/My Documents/Minor 2016/BIF 30806 Advanced Bioinformatics/Project/ballgown")
bis1 <- read.delim(file = "\\\\SCOMP0856\\veen119$\\My Documents\\Minor 2016\\BIF 30806 Advanced Bioinformatics\\Project\\ballgown\\t_data_bis1.ctab", header = TRUE)
bis2 <- read.delim(file = "\\\\SCOMP0856\\veen119$\\My Documents\\Minor 2016\\BIF 30806 Advanced Bioinformatics\\Project\\ballgown\\t_data_bis2.ctab", header = TRUE)
bis3 <- read.delim(file = "\\\\SCOMP0856\\veen119$\\My Documents\\Minor 2016\\BIF 30806 Advanced Bioinformatics\\Project\\ballgown\\t_data_bis3.ctab", header = TRUE)
cont1 <- read.delim(file = "\\\\SCOMP0856\\veen119$\\My Documents\\Minor 2016\\BIF 30806 Advanced Bioinformatics\\Project\\ballgown\\t_data_control1.ctab", header = TRUE)
cont2 <- read.delim(file = "\\\\SCOMP0856\\veen119$\\My Documents\\Minor 2016\\BIF 30806 Advanced Bioinformatics\\Project\\ballgown\\t_data_control2.ctab", header = TRUE)
cont3 <- read.delim(file = "\\\\SCOMP0856\\veen119$\\My Documents\\Minor 2016\\BIF 30806 Advanced Bioinformatics\\Project\\ballgown\\t_data_control3.ctab", header = TRUE)

dim(bis1)
dim(bis2)
dim(bis3)
dim(cont1)
dim(cont2)
dim(cont3)

length(which(bis2$gene_id != cont1$gene_id))

total_data = c(59492,6)
total_data$gene_id <- bis1$gene_id
total_data$gene_name <- bis1$gene_name

total_data$bis1 <- bis1$FPKM
total_data$bis2 <- bis2$FPKM
total_data$bis3 <- bis3$FPKM
total_data$control1 <- cont1$FPKM
total_data$control2 <- cont2$FPKM
total_data$control3 <- cont3$FPKM
total_data <- as.data.frame(total_data)
total_data <- total_data[,3:10]

test_data <- total_data[,3:8]
total_data$gene_id[which(total_data$gene_name != ".")] <- total_data$gene_name


bis1 = read.delim(file = "t_data_bis1.ctab", header = TRUE)
bis2 = read.delim(file = "t_data_bis2.ctab", header = TRUE)
bis3 = read.delim(file = "t_data_bis3.ctab", header = TRUE)
control1 = read.delim(file = "t_data_control1.ctab", header = TRUE)
control2 = read.delim(file = "t_data_control2.ctab", header = TRUE)
control3 = read.delim(file = "t_data_control3.ctab", header = TRUE)

genenames = bis1$gene_name[bis1$gene_name != "."]
bis1FPKM = bis1$FPKM[bis1$gene_name != "."]
bis2FPKM = bis2$FPKM[bis2$gene_name != "."]
bis3FPKM = bis3$FPKM[bis3$gene_name != "."]
control1FPKM = control1$FPKM[control1$gene_name != "."]
control2FPKM = control2$FPKM[control2$gene_name != "."]
control3FPKM = control3$FPKM[control3$gene_name != "."]

exp_data = data.frame(bis1FPKM, bis2FPKM, bis3FPKM, control1FPKM, control2FPKM, control3FPKM, row.names = genenames)
sds = apply(exp_data, 1, sd)
idxs = which(sds == 0)
exp_data = exp_data[-idxs, ]
dim(exp_data)
head(exp_data)

write.table(exp_data, file="gene_expression_table.csv", sep=",", row.names=TRUE)
