awk -vdis=$dis '($1==dis){print $2}' ../data/maf01_mac3.dat > ../results/dis/${dis}.markers
awk -vdis=$dis 'NR==FNR{m[$1]=1;next} ($1 in m){gs[$3]=gs[$3]"\t"$2} END{for(g in gs)print g gs[g]}' ../results/dis/${dis}.markers ../data/genes.anno > ../results/dis/${dis}.genes

Rscript step2_SPAtests.R \
        --vcfFile=../data/allchr_revel.vcf.gz \
        --vcfFileIndex=../data/allchr_revel.vcf.gz.tbi \
        --vcfField=DS \
        --IsDropMissingDosages=FALSE \
        --minMAF=0 \
        --minMAC=0.5 \
        --GMMATmodelFile=../data/dis/${dis}/${dis}.rda \
        --varianceRatioFile=../data/dis/${dis}/${dis}.varianceRatio.txt \
        --SAIGEOutputFile=../results/dis/${dis}.SAIGE.gene.txt \
        --LOCO=FALSE    \
        --groupFile=../results/dis/${dis}.genes    \
        --sparseSigmaFile=../data/dis/${dis}/${dis}.varianceRatio.txt_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseSigma.mtx       \
        --dosageZerodCutoff=0   \
        --IsSingleVarinGroupTest=TRUE > ../results/dis/${dis}.log
