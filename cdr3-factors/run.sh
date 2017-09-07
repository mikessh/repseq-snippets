rm BuildFactors.class
javac BuildFactors.java

# Fetch symlinked
cd ../samples/
SS=`readlink -f [HK]* | paste -sd "," -`  # all samples
cd ..

java -Xmx100G -cp cdr3-factors/ BuildFactors "$SS" factors.txt