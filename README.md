# snakepipe
Snakemake pipeline to process FASTQ -> Processed BAM -> basic SNV,SV, and CNV calling 


# usage without PBS/Torque
```
$ snakemake -p done
```

# usage with PBS/Torque
모든 job들을 스케쥴러에 바로 등록하기 위해서는 PBS-Torque 프로파일을 사용하여야함. `./profile` 폴더에 PBS-Torque에 submission하는 스크립트들이 이미 설치되어 있음. 

```
$ snakemake -p done --profile profile -j 24 --immediate-submit
```
## Argument 설명
-p done: 최종적으로 모든 job들이 완성이 되면 `touch done`이란 커맨드가 실행되어 `done`이라는 파일이 생성된다. -p는 `done`이라는 파일을 만들기 위해 역으로 모든 rule들을 추적하여 커맨드를 실행시키라는 명령입니다. 
-np done: 위의 명령어와 같은 것이지만 `-n`은 dryrun을 시도해봐서 문제가 없는지 확인하는 명렁입니다. 
-j 24: 동시에 사용할 코어를 최대 24개로 제한합니다. 

## DAG 그리기
```
snakemake --forceall --dag | dot -Tpng > dag.png
```
job들의 dependency를 directed acyclic graph로 visualization을 시켜줍니다. 



