library(data.table)
library(matlib)
library(Rcpp)
library(RcppArmadillo)
library(bigmemory)
library(mvtnorm)
library(MASS)
library(magic)

suppressMessages(library('plink2R'))
suppressMessages(library("optparse"))
suppressMessages(library(aSPU2))
source("dist_support.R")

# Some codes are directly copied from TWAS source codes
option_list = list(
make_option("--sumstats", action="store", default=NA, type='character',
help="summary statistics (rds file and must have SNP and Z column headers) [required]"),
make_option("--out", action="store", default=NA, type='character',
help="Path to output files [required]"),
make_option("--weights", action="store", default=NA, type='character',
help="File listing molecular weight (rds files and must have columns WGT,ID,CHR,P0,P1) [required]"),
make_option("--weights_dir", action="store", default=NA, type='character',
help="Path to directory where weight files (WGT column) are stored [required]"),
make_option("--ref_ld", action="store", default=NA, type='character',
help="Reference LD files in binary PLINK format [required]"),
make_option("--pathway_list", action="store", default=NA, type='character',
help="Pathway list we want to analyze[required]"),
make_option("--force_model", action="store", default=NA, type='character',
help="Force specific predictive model to be used, no flag (default) means select most significant cross-val. Options: blup,lasso,top1,enet")
)


opt = parse_args(OptionParser(option_list=option_list))


sumstat = readRDS(opt$sumstats)
wgtlist = read.table(opt$weights,head=T,as.is=T)

fc <- file(opt$pathway_list)
path.final <- strsplit(readLines(fc), ",")

pathway.tmp = NULL
for (job in 1:length(path.final)) {
    pathway.tmp = c(pathway.tmp, path.final[[job]])
}


pathway.tmp = unique(pathway.tmp)
wgtlist = wgtlist[wgtlist[,"ID"] %in% pathway.tmp,]


out.put.final = as.data.frame(matrix(NA,length(path.final),7))
colnames(out.put.final) = c("pathway","# genes","# SNPs","PathSPU(1)","PathSPU(2)","aSPUpath2","time")

for (job in 1:length(path.final) ) {
    tryCatch({
        pathway = path.final[[job]]
        
        wgtlist0 = wgtlist[wgtlist[,"ID"] %in% pathway,]
        
        sumstat.orgin = sumstat
        # only save related information
        snp.inf = list()
        
        chr.inf = unique(wgtlist0$CHR)
        for ( w in chr.inf )  {
            genos = read_plink(paste(opt$ref_ld,w,sep=''),impute="avg")
            snp.inf[[w]] = genos$bim[,2]
        }
        snp.inf = unlist(snp.inf)
        snp.inf = unique(snp.inf)
        
        m = match(snp.inf, sumstat.orgin$SNP)
        ## For each wgt file:
        sumstat.orgin = sumstat.orgin[m,]
        
        
        out.res = as.data.frame(matrix(NA,nrow(wgtlist0),16))
        
        n.chr = unique(wgtlist0[,"CHR"])
        
        out.gene.i = 1
        snp.name = list()
        U.out = NULL
        V.list = list()
        inf.list = list()
        out.i = 1
        weight.out = NULL
        
        for ( w in n.chr ) {
            
            tryCatch({
                
                # Load in summary stats
                sumstat = sumstat.orgin
                # Load in reference data
                # w = 1
                genos = read_plink(paste(opt$ref_ld,w,sep=''),impute="avg")
                
                # Match summary data to input, record NA where summary data is missing
                m = match( genos$bim[,2] , sumstat$SNP )
                sum.missing = is.na(m)
                sumstat = sumstat[m,]
                sumstat$SNP = genos$bim[,2]
                sumstat$A1[ sum.missing ] = genos$bim[sum.missing,5]
                sumstat$A2[ sum.missing ] = genos$bim[sum.missing,6]
                
                # QC / allele-flip the input and output
                qc = allele.qc( sumstat$A1 , sumstat$A2 , genos$bim[,5] , genos$bim[,6] )
                
                # Flip Z-scores for mismatching alleles
                sumstat$Z[ qc$flip ] = -1 * sumstat$Z[ qc$flip ]
                sumstat$A1[ qc$flip ] = genos$bim[qc$flip,5]
                sumstat$A2[ qc$flip ] = genos$bim[qc$flip,6]
                
                # Remove strand ambiguous SNPs (if any)
                if ( sum(!qc$keep) > 0 ) {
                    genos$bim = genos$bim[qc$keep,]
                    genos$bed = genos$bed[,qc$keep]
                    sumstat = sumstat[qc$keep,]
                }
                
                # Load weights
                wgtlist2 = wgtlist0[wgtlist0$CHR==w,]
                wgt.matrix2 = NULL
                snps2 = NULL
                for(i in 1:nrow(wgtlist2)) {
                    wgt.file = paste(opt$weights_dir,"/",wgtlist2$WGT[i],sep='')
                    load(wgt.file)
                    # Remove NAs (these should not be here)
                    wgt.matrix[is.na(wgt.matrix)] = 0
                    
                    # Identify the best model -- use enet
                    mod.best = (which.max(cv.performance[1,]))
                    
                    if ( names(mod.best) == "top1" ) {
                        # cat( "WARNING: top eQTL is the best predictor for this gene, continuing with 2nd-best model\n" )
                        mod.best = names( which.max(cv.performance[1,colnames(cv.performance)!="top1"]) )
                        mod.best = which( colnames(cv.performance) == mod.best )
                    }
                    wgt.matrix = wgt.matrix[,mod.best]
                    
                    denom = sum(abs(wgt.matrix))
                    if(denom==0) {
                        denom  = 0.1
                    }
                    wgt.matrix = wgt.matrix/denom
                    
                    wgt.matrix = as.matrix(wgt.matrix)
                    wgt.matrix2 = rbind(wgt.matrix2,wgt.matrix)
                    wgt.matrix2 = as.matrix(wgt.matrix2)
                    snps2 = rbind(snps2,snps)
                    
                    snp.name[[out.gene.i]] = snps[,2]
                    out.gene.i = out.gene.i + 1
                }
                
                wgt.matrix = wgt.matrix2
                snps = snps2
                
                # Match up the SNPs and weights
                m = match( snps[,2] , genos$bim[,2] )
                m.keep = !is.na(m)
                snps = snps[m.keep,]
                wgt.matrix = wgt.matrix[m.keep,]
                wgt.matrix = as.matrix(wgt.matrix)
                
                cur.genos = scale(genos$bed[,m[m.keep]])
                cur.bim = genos$bim[m[m.keep],]
                # Flip WEIGHTS for mismatching alleles
                qc = allele.qc( snps[,5] , snps[,6] , cur.bim[,5] , cur.bim[,6] )
                wgt.matrix[qc$flip,] = -1 * wgt.matrix[qc$flip,]
                rm(snps)
                rm(snps2)
                rm(wgt.matrix2)
                
                
                cur.FAIL = FALSE
                
                # Match up the SNPs and the summary stats
                m = match(cur.bim[,2] , sumstat$SNP)
                cur.Z = sumstat$Z[m]
                
                # Compute LD matrix
                if ( length(cur.Z) == 0 ) {
                    cat( "WARNING : " , unlist(wgtlist0[w,]) , " had no overlapping SNPs\n")
                    cur.FAIL = TRUE
                } else {
                    cur.LD = t(cur.genos) %*% cur.genos / (nrow(cur.genos)-1)
                    
                    cur.miss = is.na(cur.Z)
                    # Impute missing Z-scores
                    if ( sum(cur.miss) != 0 ) {
                        if ( sum(!cur.miss) == 0 ) {
                            cat( "WARNING : " , unlist(wgtlist0[w,]) , " had no overlapping GWAS Z-scores\n")
                            cur.FAIL = TRUE
                        } else {
                            cur.wgt =  cur.LD[cur.miss,!cur.miss] %*% solve( cur.LD[!cur.miss,!cur.miss] + 0.1 * diag(sum(!cur.miss)) )
                            cur.impz = cur.wgt %*% cur.Z[!cur.miss]
                            cur.r2pred = diag( cur.wgt %*% cur.LD[!cur.miss,!cur.miss] %*% t(cur.wgt) )
                            cur.Z[cur.miss] = cur.impz / sqrt(cur.r2pred)
                        }
                    }
                    
                    if ( !cur.FAIL ) {
                        
                        U = cur.Z
                        U = as.matrix(U)
                        
                        V = cur.LD
                        weight = wgt.matrix
                        name = rownames(V)
                        rownames(U) = name
                        
                        U = as.matrix(U)
                        # remove SNPs corresponding to zero weight
                        weight.tmp = abs(weight)
                        index = rowSums(weight.tmp) > 0
                        U = U[index,]
                        V = V[index,]
                        V = V[,index]
                        weight = weight[index,]
                        weight = as.matrix(weight)
                        U = as.matrix(U)
                        V = as.matrix(V)
                        V.list[[out.i]] = V
                        out.i = out.i + 1
                        weight.out = rbind(weight.out,weight)
                        U.out = rbind(U.out, U)
                        
                    }
                }
                print(w)
            }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
        }
        
        snp.len = length(unlist(snp.name))
        
        inf = as.data.frame(matrix(NA,dim(U.out)[1],4))
        colnames(inf) = c("snpindex","gene.id","chrom.id","gene.name")
        
        rs.name = rownames(U.out)
        inf[,1] = rs.name
        
        inf.tmp= inf
        for (i in 1:nrow(wgtlist0)) {
            tmp = c(i,wgtlist0[i,c("CHR","ID")])
            inf[rs.name %in% snp.name[[i]],2:4] = tmp
        }
        
        V = do.call(adiag, V.list)
        weight = weight.out
        U = U.out
        methy.info = inf
        
        
        time.start = proc.time()
        weight_diag <- diag(as.vector(weight),nrow = length(weight))
        Zstat.w <- weight_diag %*% U
        corSNP.w <- weight_diag %*% V %*% t(weight_diag)
        
        tmp.name = rownames(U)
        rownames(Zstat.w) = tmp.name
        rownames(corSNP.w) = colnames(corSNP.w) = tmp.name
        
        Zstat.w1 =Zstat.w
        corSNP.w1 = corSNP.w
        methy.info1 = methy.info
        
        asy.aSPUpath = path.asy(Zstat.w1,corSNP.w1,methy.info1)
        time.end = proc.time()
        time.elapse1 = time.end[3] - time.start[3]
        
        out =   c(pathway[1],length(unique(methy.info[,4])),dim(U)[1],asy.aSPUpath,time.elapse1)
        
        out.put.final[job,] =  out
        print(job)
        
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
    
}

saveRDS(out.put.final,opt$out)
write.table(out.put.final, "output.txt")



