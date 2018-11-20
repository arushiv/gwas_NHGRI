All associations v1.0.2 - with added ontology annotations, GWAS Catalog study accession numbers and genotyping technology -
- Data downloaded from https://www.ebi.ac.uk/gwas/docs/file-downloads on 06-12-2018
```
wget https://www.ebi.ac.uk/gwas/api/search/downloads/alternative -O gwas_catalog_v1.0.2-associations_e92_r2018-05-29.tsv
```


Prune: Take GWAS studies with European population and prune using 1000g to r2<0.2

```
snakemake --cluster-config configuration_for_pruning/cluster.yaml  --cluster "sbatch --time {cluster.time} --mem {cluster.mem} --cpus-per-task {cluster.cpus} --job-name {cluster.jobname} -o {cluster.output} -e {cluster.error} --parsable" -j 60 -nps Snakefile --configfile configuration_for_pruning/config.yaml --keep-going
```