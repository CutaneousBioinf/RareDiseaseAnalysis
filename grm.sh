~/plink/plink2 --bfile allchr --exclude <(awk '$1==6&&$4>=26000000&&$4<=34000000||$1>22{print $2}' allchr.bim) --keep ../data/gen.pats --freq --memory 10240 --out ../data/ma

~/plink/plink2 --bfile allchr --extract <(awk 'NR>1&&$5>=0.05&&$5<=0.95{print $2}' ../data/ma.afreq) --indep-pairwise 100 5 0.2 --memory 10240 --out ../data/ma

~/plink/plink2 --bfile allchr --keep ../../data/gen.pats --extract ../data/ma.prune.in --make-bed --memory 10240 --out ../data/ma.prune

Rscript createSparseGRM.R     \
               --plinkFile=../data/ma.prune \
               --nThreads=4  \
               --outputPrefix=../data/sparseGRM        \
               --numRandomMarkerforSparseKin=2000      \
               --relatednessCutoff=0.125
