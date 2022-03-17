#!/usr/bin/Rscript
countAll=502414
mendDis=read.table("../data/blairMendelian.dat",sep="\t")
rownames(mendDis)=mendDis[,1];mendDis=mendDis[,2,drop=F]
commDis=read.table("../data/blairCommon.dat",sep="\t")
rownames(commDis)=commDis[,1];commDis=commDis[,2,drop=F]
counts=read.table(pipe("wc -l ../data/blair/* | awk '{print $2,$1}' | xargs -I{} basename {} | head -n -1"))
rownames(counts)=counts[,1]; counts=counts[,-1,drop=F]
mendKeep=paste0("m",c(17,20,34,41,47,51,52,53,56,69,70,72,76,77,91))

enrich=function(){
for(x in 1:ncol(data))
        for(y in 1:nrow(data)){
                fish=fisher.test(rbind(c(data[y,x],counts[colnames(data)[x],]-data[y,x]),c(counts[rownames(data)[y],]-data[y,x],countAll-counts[colnames(data)[x],]-counts[rownames(data)[y],]+data[y,x])))
                or[y,x]<<-fish$estimate
                p[y,x]<<-fish$p.value
                if(p[y,x]>=.05)or[y,x]<<-NA
        }
}

data=read.table("../data/blairMendMend.dat",sep="\t")
data=data[mendKeep,mendKeep]
or=p=matrix(0,ncol=ncol(data),nrow=nrow(data),dimnames=list(mendDis[rownames(data),],mendDis[colnames(data),]))
enrich()
write.table(or,"../data/blairMendMend_or.dat",quote=F,sep="\t")
write.table(p,"../data/blairMendMend_p.dat",quote=F,sep="\t")

data=read.table("../data/blairMendComm.dat",sep="\t")
data=data[mendKeep,]
or=p=matrix(0,ncol=ncol(data),nrow=nrow(data),dimnames=list(mendDis[rownames(data),],commDis[colnames(data),]))
enrich()
write.table(or,"../data/blairMendComm_or.dat",quote=F,sep="\t")
write.table(p,"../data/blairMendComm_p.dat",quote=F,sep="\t")
