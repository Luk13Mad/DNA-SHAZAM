# DNA-SHAZAM
prototype of adaptation of the Shazam algorithm for DNA

FASTA files of mtDNA can for example be obtained from MITOMAP database (https://www.mitomap.org) 

implementation in rust as command line tool [here](https://github.com/Luk13Mad/dnasr)

-----------------------------------------------------------
to run on a FASTA file with methylated data execute:

python3 Shannon_entropy_new.py -in $input.fasta -m Y -w 10 -k 5 -np 15  
for i in *[0-9]_meth.csv.gz ; do python3 feature_extraction_new.py -in $i -w 24 -s 2; done   
for i in *[0-9]_meth.csv.gz.featvec.csv.gz ; do python3 generate_adresses_new.py -in $i -db Y -off 1 ; done  
for i in adress_list_db*[0-9].csv.gz ; do zcat $i >> merged_db_meth.csv ; done  
