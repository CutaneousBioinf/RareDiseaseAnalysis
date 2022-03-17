mendDis=$(wc -l ../data/blair/* | awk '/m/&&$1!=0{print $2}')
commDis=$(wc -l ../data/blair/* | awk '/c/&&$1!=0{print $2}')

{ echo "$mendDis" | xargs -I{} basename {} | tr '\n' '\t' | sed 's/\t$/\n/'
for x in $mendDis; do
        basename $x | tr '\n' '\t'
        for y in $mendDis; do
                awk 'NR==FNR{pat[$1]=1;next} ($1 in pat)' $x $y | wc -l | tr '\n' '\t'
        done | sed 's/\t$/\n/'
done; } > ../data/blairMendMend.dat

{ echo "$commDis" | xargs -I{} basename {} | tr '\n' '\t' | sed 's/\t$/\n/'
for x in $mendDis; do
        basename $x | tr '\n' '\t'
        for y in $commDis; do
                awk 'NR==FNR{pat[$1]=1;next} ($1 in pat)' $x $y | wc -l | tr '\n' '\t'
        done | sed 's/\t$/\n/'
done; } > ../data/blairMendComm.dat
