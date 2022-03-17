~/plink/plink2 --bed ukb.bed --bim chr9.bim --fam ukb.fam --extract <(echo 9:5073770:G:T) --export A --out ../data/jak2.mut
sed 1d ukb.txt | cut -f$(head -1 ukb.txt | tr '\t' '\n' | awk '$1=="21022-0.0"{print NR+1}' | tr '\n' ',' && echo 2) | awk 'NR==FNR{age[$1]=$2;next} FNR>1&&$NF!="NA"&&age[$2]{print $2"\t"$NF"\t"age[$2]}' - ../data/jak2.mut.raw > ../data/jak2.mut
