###########  Load required packages  #############

require(gstudio)
require(reshape2)
require(ggplot2)
require(xlsx)
require(MASS)
require(grid)

Gdat <- function() {
  slct=1
	while(slct>0 || is.nan(slct) || is.na(slct))
	{
		cat("\n\n\t\t\tYou must complete the following tasks in succession\n\n0) Move to data analysis!\n\n1) Select a folder where you would like to save your data\n
2) Specify the location of your data file\n\n3) Load your data into R and start Forensic Pipe\n\n4) Set the names of your populations\n\n")
		slct <- as.numeric(readline(prompt="Please enter your choice: "))
		if(is.nan(slct) || is.na(slct))
		{
			cat("\a\n\n\t\t#####################################################\n\n\t\t\tYou have entered an invalid choice\n\n\t\t#####################################################\n\n")
		}
		else if(slct==1)
		{
		swd <- choose.dir(caption="Select save directory")
		setwd(swd)
		}

########### Read specified data file##############

		else if(slct==2)
		{
			n1 <- file.choose()
		}

###########     Import data from    ##############
###########     file into gstudio   ##############


		else if(slct==3)
		{
			data <- read_population(n1, type = "genepop")
		}

###########   Change all pop names  ##############


		else if(slct==4)
		{
			d <- unique(data$Population)
			for (i in 1:length(d))
			{
				x <- readline(prompt="The name of population: ")
				data$Population <- gsub(paste("Pop-", i, sep=""), x, data$Population)
			}
		}

		else if(slct==0)
		{
			break
		}
		else
		{
		cat("You have entered an incorrect value")
		}
	}
  return(data)
}

ForenStat <- function(y){
###########Compute allele frequencies#############
###########Number of unique Alleles##############
	Allnum <- A(y)
###########Allele frequencies of all##############
	x <- frequencies(y)
	freqall.table <- (xtabs(Frequency ~ Allele + Locus, data = x))
	x.melt <- melt(x)
###########Allele frequencies of pops#############
	pops <- partition(y, stratum = "Population")
	FLp <- frequencies(y, stratum="Population")
	freqpop.table2 <- lapply(pops, function(x){(xtabs(Frequency ~ Allele + Locus, data = FLp))})
	FLp.melt <- melt(FLp)

	nc <- length(names(y))
	nb <- length(names(pops))

	freqpop <- data.frame()
	freqpop.table <- list()

		for(i in 1:nb)
		{
			freqpop <- FLp[which(FLp$Stratum == (unique(FLp$Stratum))[[i]]),]
			freqpop.table[[i]] <- xtabs(Frequency ~ Allele + Locus, data = freqpop)
		}

	fpt <- matrix(freqpop.table)
	
########## Find duplicated genotypes ############

		Dupli_geno <- y[duplicated(y[3:20]) | duplicated(y[,3:20], fromLast=TRUE),1:20]

	
	
###########       Polymorphic        #############
###########   Information Capacity   #############

	PIC <- fpt
	PICa <- data.frame()
	PICb <- data.frame()
	PIC2 <- list()
	PIC4 <- list()
	PIC_tot_temp <- data.frame()
	PIC_tot <- list()
	PIC_tot2 <- list()

	for (i in 1:nrow(fpt))
	{
		for(j in 1:length(fpt[[1,1]][1,]))
		{
			PICa <- apply(PIC[[i,1]], 2, function(x) sum(x*x))
			PICb <- apply(PIC[[i,1]], 2, function(x) sum(x*x*x*x))
			PIC2[[i]] <- PICa
			PIC4[[i]] <- PICb

			PIC_tot_temp <- 1 - PIC2[[i]][j] - (PIC2[[i]][j]*PIC2[[i]][j]) + PIC4[[i]][j]
			PIC_tot[j] <- PIC_tot_temp
		}
		PIC_tot2[[i]] <-PIC_tot
	}

	PIC_mat <- matrix(PIC_tot2)
	PICDF <- do.call(rbind,PIC_tot2)
	rownames(PICDF) <- names(pops)
	colnames(PICDF) <- names(fpt[[1,1]][1,])
	PICDFt <- t(PICDF)
	
###########  Expected heterozygosity #############

	hexp <- lapply(pops, function(x) return(He(x)))

	d <- unique(names(pops))
	hexp_data <- data.frame()
	hexp_data <- data.frame(Loci = hexp[[1]]$Locus)
	p <- data.frame()
		for(i in 2:(length(d)+1))
		{
			p <- paste('hexp.', d[i-1], sep="")
			hexp_data[,i] <- data.frame(p = hexp[[i-1]]$He)
			names(hexp_data)[i] <- p[]
		}

	shexp_data <- t(hexp_data[order(hexp_data[,1]),])
	
########### Observed heterozygosity  #############

	hobs <- lapply(pops, function(x) return(Ho(x)))

	d <- unique(names(pops))
	hobs_data <- data.frame()
	hobs_data <- data.frame(Loci = hobs[[1]]$Locus)
	p <- data.frame()
		for(i in 2:(length(d)+1))
		{
			p <- paste('hobs.', d[i-1], sep="")
			hobs_data[,i] <- data.frame(p = hobs[[i-1]]$Ho)
			names(hobs_data)[i] <- p[]
		}

	shobs_data <- t(hobs_data[order(hobs_data[,1]),])

#####Frequency of heterozygotes of each loci#####

	Hobsmelt <- melt(hobs_data)
	Hexpmelt <- melt(hexp_data)

###########       Homozygosity       #############
########### Frequency of homozygotes #############
###########   and PI of each loci    #############

	Homomelt = Hobsmelt
	Homomelt$value <- (1-Hobsmelt$value)

	homozy <- data.frame(hobs_data)
	PI <- data.frame(hobs_data)
	rownames(homozy) <- hobs_data[,1]
		for (i in 1:nrow(homozy))
		{
			for (j in 2:length(homozy))
			{
				homozy[[i,j]] <- (1 - (homozy[[i,j]]))
				PI[[i,j]] <- (1/(2*(homozy[i,j])))
			}
		}
	homozy <- t(homozy)
	homozy2 <- homozy[2:nrow(homozy),]
	homozy2 <- homozy2[,order(colnames(homozy2))]
	PI.melt <- melt(PI)

########### Combined Paternity index #############

	CPI <- data.frame(colnames=names(pops), CPI=1)

		for(i in 2:length(PI))
		{
			CPI[i-1,2] <- (prod(PI[i]))
		}

###########    Power of Exclusion    #############

	pexcl_pops <- homozy2
	nc1 <- length(pexcl_pops[1,])
	nb1 <- nrow(pexcl_pops)

		for(i in 1:nb1)
		{
			for(j in 1:nc1)
			{
				hom <- as.numeric(pexcl_pops[i,j])
				he <- 1 -as.numeric(pexcl_pops[i,j])
				pexcl_pops[[i,j]] <- (he*he)*(1-(2*(he*(hom*hom))))
			}
		}

	pexcl_tot_pop <- sprintf("%.12f", apply(pexcl_pops, 1, function(x) (1-prod(1-as.numeric(x)))))
	pexcl_pops.melt <- melt(pexcl_pops)
	PE <- as.numeric(pexcl_pops.melt[,3])

###########   Genotype Frequencies   #############

	temp <- list()
	all <- list()

		for(i in 1:length(names(pops)))
		{
			nc <- 3
			for(j in 3:length(names(y)))
			{
				nr <- nrow(genotype_frequencies(pops[[i]][,j]))
				vec <- data.frame(x=numeric(nr), y=numeric(nr), z=numeric(nr))
				z <- (pops[[i]][,j])
				vec <- genotype_frequencies(z)
				temp[[j]] <- vec
			}
			all[[i]] <- temp
		}

	DF <- do.call(cbind,all)

########### Random Match Probability #############

	nc <- length(names(pops))
	nr <- length(names(y))
	RMP <- matrix(ncol=nc, nrow=nr)
	rownames(RMP) <- names(y)
	colnames(RMP) <- names(pops)

		for(x in 1:nr)
		{
			for(z in 1:nc)
			{
				RMP[x,z] <- (sum((DF[[x,z]][,2]/nrow(pops[[z]]))^2))
			}
		}

	RMP.df <- as.data.frame(RMP)
	RMP.df1 <- RMP.df[3:nr,]
	RMP.df2 <- t(RMP.df1[order(rownames(RMP.df1)),order(colnames(RMP.df1)) ])

###########Combined Match Probability#############

	CMP <- data.frame(apply(RMP.df2, 1, function(x) prod(x)))
	colnames(CMP) <- "CMP"


###########Graphing Match Probability#############
###########         per loci         #############

	GF.rmp.melt <- data.frame()
	GF.rmp.melt <- melt(RMP.df2)
	GF.rmp.melt$Y1 <- cut(GF.rmp.melt$value,breaks = c(-Inf,0:0.05,0.05:0.10,0.10:0.15,0.15:0.20,0.20:0.25,0.25:0.30,0.30:0.35,0.35:0.40,0.40:0.45,0.45:0.50,0.50:0.55,0.55:0.60,0.60:0.65,0.65:0.70,Inf),right = FALSE)

###########  Discrimination Capacity #############


###########  Discrimination Capacity #############
###########         per loci         #############

	DCpl <- RMP.df2

		for(x in 1:length(RMP.df2[2,]))
		{
			for(z in 1:length(RMP.df2[,1]))
			{
				DCpl[[z,x]] <- (1 - RMP.df2[z,x])
			}
		}

	DCpl = t(DCpl)

###########     plot DC per loci     #############

	DCpl.melt <- melt(t(DCpl))

###########        Combined DC       #############
###########         per loci         #############

	CDC <- CMP
		for(i in 1:nrow(CMP))
		{
			CDC[i,1] <- sprintf("%.20f", 1 - CMP[i,1])
		}

	
	
###########       compute hwe         #############
	
	HWE <- hwe(y, mode = c("Chi","Permute"), supress_warnings = TRUE)

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	#####   print outputs #######
	
	
	
	
	write.csv(Dupli_geno, paste0(getwd(),"/dpulications.csv"), row.names=TRUE)
	write.xlsx(Allnum, paste0(getwd(),"/All_num.xlsx"), sheetName = "Alleles", append = FALSE, row.names=FALSE, col.names=TRUE)
	write.csv(freqall.table, paste0(getwd(),"/fpt_all.csv"), row.names=TRUE)
	
	tiff(paste0(getwd(),"/FreqLoci_all.tiff"), height = 6, width = 6, units = 'in', res=300)
	  print(ggplot(FLp.melt, aes(x = Locus, y = Allele)) +
	          ggtitle("Allele frequencies of all loci of all individuals") +
	          geom_tile(aes(fill = value), colour = "black") +
	          scale_fill_gradient(low = "yellow", high = "red") +
	          theme(plot.title = element_text(lineheight=.8, face="bold"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(), axis.text.x=element_text(angle=90, hjust = 1, size = 5), axis.text.y=element_text(size = 4, hjust = 1)))
	dev.off()
	
	###########  Plot Allele frequencies #############
	
	for (i in 3:length(y))
	{
	  tiff(paste0(getwd(),"/plot_", names(y[i]), ".tiff", sep=""), height = 6, width = 6, units = 'in', res=300)
	    flname <- paste("Allele frequency of", names(y[i]), "within all individuals", sep = " ")
	    print(plot(y[[i]]) + ggtitle(paste0(flname)) + 
	          theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10)))
	  dev.off()
	  
	  tiff(paste0(getwd(),"/FreqLoci_pop.tiff"), height = 10, width = 12, units = 'in', res=300)
	  print(ggplot(FLp.melt, aes(x = Locus, y = Allele)) +
	          facet_wrap(~Stratum, ncol=nb) + ggtitle("Allele frequencies of all loci within populations") +
	          geom_tile(aes(fill = value), colour = "black") +
	          scale_fill_gradient(low = "yellow", high = "red") +
	          theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(), axis.text.x=element_text(angle=90, hjust = 1, size = 6), axis.text.y=element_text(hjust = 1, size = 6)))
	  dev.off()
	}
	
	d <- unique(FLp[,2])
	for (i in 1:length(d))
	{
	  tiff(paste0(getwd(),"/freqsplot_", d[[i]], ".tiff", sep=""), height = 6, width = 6, units = 'in', res=300)
	  flname <- paste("Allele frequency of", d[[i]], "between populations", sep = " ")
	  f <- FLp[FLp$Locus %in% c(d[i]),]
	  print(ggplot(f) + ggtitle(paste0(flname)) + geom_frequencies(f) + facet_grid(Stratum ~ .) + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6)))
	  dev.off()
	}
	###########   Tabulate frequencies   #############
	
	rownames(fpt) <- names(pops)
	for(i in 1:nb)
	{
	  fpt2 <- do.call(cbind,fpt[i])
	  fpt2 <- fpt2[,order(colnames(fpt2)) ]
	  write.xlsx(fpt2, paste0(getwd(),"/fpt", names(pops)[[i]], ".xlsx"), sheetName = "fpt", append = TRUE, row.names=TRUE)
	}
  
		
	write.xlsx(PICDFt, paste0(getwd(),"/PICDF.xlsx"), sheetName = "PIC", append = TRUE, row.names=TRUE, col.names=TRUE)
	
	write.xlsx(shexp_data, paste0(getwd(),"/hobs_data.xlsx"), sheetName = "ExpHe", col.names=FALSE, row.names=TRUE)
	write.xlsx(shobs_data, paste0(getwd(),"/hobs_data.xlsx"), sheetName = "ObsHe", append = T ,col.names=FALSE, row.names=TRUE)
	
	tiff(paste0(getwd(),"/Expected_heterozygosity.tiff"), height = 6, width = 6, units = 'in', res=300)
	  print(ggplot(Hexpmelt, aes(x = Loci, y = variable)) +
	        ggtitle("Expected heterozygosity of all loci in each population") +
	        geom_tile(aes(fill = value), colour = "black", height = 1) +
	        scale_fill_gradient(low = "yellow", high = "red") +
	        theme(axis.text.x=element_text(angle=90, hjust = 1),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank())+
	        ylab("Population groups") +
	        xlab("Locus"))
	dev.off()
	
			
	tiff(paste0(getwd(),"/Observed_heterozygosity.tiff"), height = 6, width = 6, units = 'in', res=300)
			print(ggplot(Hobsmelt, aes(x = Loci, y = variable)) +
			        ggtitle("Observed heterozygosity of all loci in each population") +
			        geom_tile(aes(fill = value), colour = "black", height = 1) +
			        scale_fill_gradient(low = "yellow", high = "red") +
			        theme(axis.text.x=element_text(angle=90, hjust = 1),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank())+
			        ylab("Population groups") +
			        xlab("Locus"))
	dev.off()
			
	write.xlsx(homozy2, paste0(getwd(),"/hobs_data.xlsx"), sheetName = "Homozygotes", append = TRUE, col.names= TRUE, row.names=TRUE)
			
	tiff(paste0(getwd(),"Observed_homozygosity.tiff"), height = 6, width = 6, units = 'in', res=300)
	  	print(ggplot(Homomelt, aes(x = Loci, y = variable)) + ggtitle("Observed heterozygotes of all loci in each population") + geom_tile(aes(fill = value), colour = "black", height = 1) + scale_fill_gradient(low = "yellow", high = "red") + theme(axis.text.x=element_text(angle=90, hjust = 1),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank())+ ylab("Population groups") + xlab("Locus"))
	dev.off()
			
	write.xlsx(PI, paste0(getwd(),"/Paternity_Index.xlsx"), sheetName = "paternity index", row.names=FALSE)
			
	tiff(paste0(getwd(),"/Pat_Index.tiff"), height = 6, width = 6, units = 'in', res=300)
	  	print(ggplot(PI.melt, aes(x = Loci, y = variable)) + ggtitle("The paternity index of all loci in each population") + geom_tile(aes(fill = value), colour = "black", height = 1) + scale_fill_gradient(low = "yellow", high = "red") + theme(axis.text.x=element_text(angle=90, hjust = 1),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank()))
	dev.off()
			
	write.xlsx(CPI, paste0(getwd(),"/Paternity_Index.xlsx"), sheetName = "CPI", row.names=FALSE, append = TRUE)
			
	write.xlsx(RMP.df2, paste0(getwd(),"/RMP.xlsx"), col.names=TRUE, sheetName = "RMP", append = TRUE, row.names=TRUE) 
			
	write.xlsx(CMP, paste0(getwd(),"/RMP.xlsx"), col.names=TRUE, sheetName = "CMP", append = TRUE, row.names=TRUE)

	tiff(paste0(getwd(),"/RMP_heatmap.tiff"), height = 6, width = 6, units = 'in', res=300)
	  	 print(ggplot(GF.rmp.melt, aes(x = Var2, y = Var1)) + ggtitle("Random match probability of loci between populations") +	geom_tile(aes(fill = Y1), colour = "white") + scale_fill_brewer(palette = 1) + theme(axis.text.x=element_text(angle=90, hjust = 1),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank()) + ylab("Population groups") + xlab("Locus"))
	dev.off()
		
	write.xlsx(DCpl, paste0(getwd(),"/RMP.xlsx"), col.names=TRUE, sheetName = "DC", append = TRUE, row.names=TRUE)
	
	tiff(paste0(getwd(),"DC_heatmaps.tiff"), height = 6, width = 6, units = 'in', res=300)
	    print(ggplot(DCpl.melt, aes(x = Var2, y = Var1)) + ggtitle("Discrimination capacity of each loci within populations")+ geom_tile(aes(fill = value), colour = "black", height = 1) + scale_fill_gradient(low = "yellow", high = "red") + theme(axis.text.x=element_text(angle=90, hjust = 1),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank()) + ylab("Population groups") + xlab("Locus"))
	dev.off()
			
	write.xlsx(CDC, paste0(getwd(),"/RMP.xlsx"), sheetName = "CDC", col.names=FALSE, append = TRUE, row.names=TRUE)
		  
	write.xlsx(pexcl_pops, paste0(getwd(),"/pexcl_all_pops.xlsx"), sheetName = "pexcl_pops", col.names=TRUE, append = TRUE, row.names=TRUE)
			
	write.xlsx(pexcl_tot_pop, paste0(getwd(),"/pexcl_all_pops.xlsx"), sheetName = "CPE_pops", col.names=TRUE, append = TRUE, row.names=TRUE)

	tiff(paste0(getwd(),"/PE_heatmap.tiff"), height = 6, width = 6, units = 'in', res=300)
		  print(ggplot(pexcl_pops.melt, aes(x = Var2, y = Var1)) +
		          ggtitle("Power of exclusion of loci within populations") +
		          geom_tile(aes(fill = PE), colour = "black", height = 1) +
		          scale_fill_gradient(low = "yellow", high = "red") +
		          theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(), axis.text.x=element_text(angle=90, hjust = 1)) +
		          ylab("Population groups") +
		          xlab("Locus"))
	dev.off()
		
	write.xlsx(HWE, paste0(getwd(),"/hobs_data.xlsx"), sheetName = "HWE", append = TRUE, row.names=FALSE, col.names=TRUE)
	
}



