wget https://github.com/antigenomics/vdjdb-db/archive/master.zip
unzip master.zip
cd vdjdb-db-master/src/
groovy -cp . BuildDatabase.groovy
cd ../../
cp vdjdb-db-master/database/vdjdb.slim.txt .
rm master.zip
rm -r vdjdb-db-master
