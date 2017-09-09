### SSU #################

SumSqU<-function(U, CovS,cr){
    if (is.null(dim(CovS))) {# only one-dim:
        Tscore<- sum(U^2 /CovS)
        if (is.na(Tscore) || is.infinite(Tscore) || is.nan(Tscore)) Tscore<-0
        pTg1<-as.numeric(1-pchisq(Tscore, 1))
    }
    else {
        #it's possible U=0 and Cov(U)=0:
        if (all(abs(U)<1e-20)) pTg1<-1 else{
            Tg1<- t(U) %*% U
            ##distr of Tg1 is sum of cr Chisq_1:
            #cr<-eigen(CovS,symmetric = TRUE, only.values=TRUE)$values
            
            ##approximate the distri by alpha Chisq_d + beta:
            alpha1<-sum(cr*cr*cr)/sum(cr*cr)
            beta1<-sum(cr) - (sum(cr*cr)^2)/(sum(cr*cr*cr))
            d1<-(sum(cr*cr)^3)/(sum(cr*cr*cr)^2)
            alpha1<-as.double(alpha1)
            beta1<-as.double(beta1)
            d1<-as.double(d1)
            pTg1<-as.numeric(1-pchisq((Tg1-beta1)/alpha1, d1))
        }
    }
    return(pTg1)
}

##########SumTest########################
Sum<-function(U, CovS){
    #it's possible U=0 and Cov(U)=0:
    if (all(abs(sum(U))<1e-20)) pTsum<-1 else{
        a<-rep(1, length(U))
        Tsum<- sum(U)/(sqrt(as.numeric(t(a) %*% CovS %*% (a))))
        pTsum <- as.numeric( 1-pchisq(Tsum^2, 1) )
    }
    pTsum
}

path.asy <- function(U, V, methy.info) { #cr is the eigenvalue of CovS
    
    chrs <- unique(methy.info[,"chrom.id"])
    cr = NULL
    for (i in 1:length(chrs)) {
        c = chrs[i]
        Covtemp = unname(V[methy.info[methy.info[,"chrom.id"]==c,1],methy.info[methy.info[,"chrom.id"]==c,1]])
        eV <- eigen(Covtemp,symmetric = TRUE, only.values=TRUE)$values
        cr = c(cr,eV)
    }
    
    p.sum = Sum(U,V)
    p.ssu = SumSqU(U,V, cr)
    p.final = 1 - (1- min(p.sum,p.ssu))^2
    c(p.sum,p.ssu,p.final)
}

# This function (allele.qc) is downloaded from TWAS (http://gusevlab.org/projects/fusion/#typical-analysis-and-output)
allele.qc = function(a1,a2,ref1,ref2) {
    ref = ref1
    flip = ref
    flip[ref == "A"] = "T"
    flip[ref == "T"] = "A"
    flip[ref == "G"] = "C"
    flip[ref == "C"] = "G"
    flip1 = flip
    
    ref = ref2
    flip = ref
    flip[ref == "A"] = "T"
    flip[ref == "T"] = "A"
    flip[ref == "G"] = "C"
    flip[ref == "C"] = "G"
    flip2 = flip;
    
    snp = list()
    snp[["keep"]] = !((a1=="A" & a2=="T") | (a1=="T" & a2=="A") | (a1=="C" & a2=="G") | (a1=="G" & a2=="C"))
    snp[["flip"]] = (a1 == ref2 & a2 == ref1) | (a1 == flip2 & a2 == flip1)
    
    return(snp)
}
