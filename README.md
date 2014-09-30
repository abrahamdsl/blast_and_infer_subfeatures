blast_and_infer_subfeatures
===========================

A (quite messy but usable) Perl program to infer gene and its subfeatures on a (longer, not sure with shorter) suspected corresponding gene sequence in a new, but almost identical genome. Uses NCBI BLAST+ whenever possible.

How to use.

1) Have a MySQL database. Then create a database and insert the SQL dump db_schema.sql.
2) Edit the Perl script with your database settings. ( Lines 54, 312-314 ).
3) Run Perl program according to the parameters needed. 

Ex:

  ./blast_and_infer_genes_and_subfeatures.pl --genome1=../../../Sbicolor_v2.1_255.fa  \
    --genome2=thisnewgenome.fa \
    --gff1=../Sbicolor_v2.1_255_gene_exons.gff3  \
    --gff2=../somegenesfromsomeweirdprocess.gff3
  
4) Wait for so long. Well, currently addressing that. In my run, 30 hours were spent for 5349 genes. I recommend
running the program under nohup.
