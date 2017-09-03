rm PublicCoincidence.class 
javac PublicCoincidence.java

# Fetch symlinked
cd ../samples
SS=`ls [HK]*.txt | paste -sd "," -`  # all samples
cd ..

# Other parameters:
# freq threshold for publics, P&log10 odds threshold for association - comma-separated
# file with pooled clonotypes, incl incidence count

java -Xmx100G -cp coincidence/ PublicCoincidence 0.05,0.001,0.3 hip.pool.txt "$SS" hip_assoc