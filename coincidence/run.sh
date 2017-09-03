rm PublicCoincidence.class 
javac PublicCoincidence.java

# Fetch symlinked
SS=`readlink -f ../samples/[HK]* | paste -sd "," -`  # all samples
SP=`readlink -f ../hip.pool.txt`

# Other parameters:
# freq threshold for publics, P&log10 odds threshold for association - comma-separated
# file with pooled clonotypes, incl incidence count

java -Xmx100G -cp . PublicCoincidence 0.05,0.001,0.3 "$SP" "$SS" hip_assoc