awk '{print $1"\t"$2"\t"$4"\t"$6"\t"$5}' allchr.bim > variants.txt
java -Xmx8g -jar snpEff.jar GRCh38.99 variants.txt > ../data/anno.vcf

java -jar SnpSift.jar dbnsfp -db ../data/dbNSFP4.1a.txt.gz ../data/anno.vcf > ../data/anno.nsfp.vcf
bgzip ../data/anno.nsfp.vcf && tabix -pvcf ../data/anno.nsfp.vcf.gz
