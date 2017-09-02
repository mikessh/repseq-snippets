rm PublicCoincidence.class 
javac PublicCoincidence.java

SS=`ls ../H*.txt | paste -sd "," -`  # all samples

# Other parameters:
# freq threshold for publics, P&log10 odds threshold for association - comma-separated
# file with pooled clonotypes, incl incidence count

java -Xmx100G -cp . PublicCoincidence 0.1,0.001,0.3 ../pool/pool.aa.table.txt "$SS" hip_assoc