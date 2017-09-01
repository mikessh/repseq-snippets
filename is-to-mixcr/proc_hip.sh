parallel -j 10 --workdir $PWD bash ./mixcr_convert.sh {} ::: `ls ../*.tsv`

#for f in `ls ../*.tsv`
#do
#   echo $f
#   bash mixcr_convert.sh $f
#done