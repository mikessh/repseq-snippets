#!/bin/bash
IN=$1
OUT=${IN##*/}

groovy Immunoseq2Fastq.groovy $IN $OUT.fq

mixcr -Xmx8G align --save-description -f $OUT.fq $OUT.vdjca
mixcr -Xmx8G exportAlignments -descrR1 -nFeature CDR3 -aaFeature CDR3 \
                       -vFamilies -dGenes -jFamilies -positionOf CDR3Begin \
                       -positionOf VEndTrimmed -positionOf DBeginTrimmed -positionOf DEndTrimmed -positionOf JBeginTrimmed -f \
                       $OUT.vdjca $OUT.mixcr

rm $OUT.fq $OUT.vdjca

Rscript finalize.R $OUT.mixcr

# some bug with family names in mixcr
sed -ri 's/\*00//g' $OUT.mixcr