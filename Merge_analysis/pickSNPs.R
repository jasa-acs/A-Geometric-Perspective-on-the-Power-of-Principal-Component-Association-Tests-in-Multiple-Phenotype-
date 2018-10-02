### http://www.gettinggeneticsdone.com/2011/03/prune-gwas-data-in-rstats.html
#-------------------------------------------------------------------------------------
#   Select a set of SNPs > dist bp (default 100kb) apart
#   Map is a matrix with colnames "chrom" and "position".
#   Matrix MUST BE sorted by chrom,position!!!!!!!!!!!!!!!!!
#   Can have more columns but the colnames "chrom" and "position" can't be changed
#   The function returns a vector of row indices corresponding
#   to the rows you should pick for your subset.
#-------------------------------------------------------------------------------------
pickSNPs<-function(map, dist = 100000) { 	   
	t=as.data.frame(table(map$chrom))				
	vec = map$position
	subs = c(1,rep(NA,nrow(map)-1));  # length(subs) = nrow, but the 1st element is 1 => always select the 1st snp  
	for (k in 1:nrow(t)) { #  t: count of SNPs per chr
		if (k==1) i=1 else i=sum(t[1:(k-1),2])+1; # the 1st SNP on each ch
		subs[i] = i
		stop=sum(t[1:k,2])
		while (i<stop) {
			for (j in (i+1):stop) {
				if ((vec[j]-vec[i]) > dist) {
						#cat(i, vec[i], j, vec[j],vec[j]-vec[i], x[j],'\n'); 
						subs[j]= j; 
						i=j; 
						next; # jump out of loop
				} else if (j==stop) i=stop
			}
		}
	}
	subs[!is.na(subs)]	#  row number of selected SNPs
}
#-------------------------------------------------------------------------------------