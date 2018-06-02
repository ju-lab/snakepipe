# snakepipe
Snakemake pipeline to process FASTQ -> Processed BAM -> basic SNV,SV, and CNV calling 


# usage without PBS/Torque
```
$ snakemake -p done
```

# usage with PBS/Torque
```
$ snakemake -p done --cluster 'python qsnake.py' -j 100
```



