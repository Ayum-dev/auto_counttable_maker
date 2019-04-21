#! bin/bash
set -xeu

<<COMMENTOUT

$ bash MakeCountTable_Ion_fastq.sh Ion_fastq.csv mouse

- fastqかSRRの判別
- trimmomatic
- gtf, transcript file をGENCODEから
- salmon

COMMENTOUT


# 実験テーブル.csv
EX_MATRIX_FILE=$1
RUNINDOCKER=1
THREADS=4
REF_SPIECE=$2

DOCKER=docker
# DOCKER=udocker # udockerも指定できる。

SCRIPT_DIR=$(cd $(dirname $0); pwd)

if [[ $REF_SPIECE = mouse ]]; then
  BASE_REF_TRANSCRIPT=ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M19
  REF_TRANSCRIPT=gencode.vM19.transcripts.fa.gz
  SALMON_INDEX=salmon_index_mouse
#   REF_GTF=gencode.vM19.annotation.gtf.gz
  TX2SYMBOL=gencode.vM19.metadata.MGI.gz

elif [[ $REF_SPIECE = human ]]; then
  BASE_REF_TRANSCRIPT=ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29
  # REF_TRANSCRIPT=gencode.v29.pc_translations.fa.gz
  REF_TRANSCRIPT=gencode.v29.transcripts.fa.gz
  SALMON_INDEX=salmon_index_human
#   REF_GTF=gencode.v29.annotation.gtf.gz
  TX2SYMBOL=gencode.v29.metadata.HGNC.gz
  
else
  echo No reference speice!
  exit
fi

COWSAY=cowsay
PREFETCH=prefetch
PFASTQ_DUMP=pfastq-dump
FASTQ_DUMP=fastq-dump
FASTQC=fastqc
MULTIQC=multiqc
FASTXTRIMMER=fastx_trimmer
FASTQQUALITYTRIMMER=fastq_quality_trimmer
# TRIMMOMATIC=trimmomatic
SALMON=salmon
RSCRIPT_TXIMPORT=Rscript


if [[ "$RUNINDOCKER" -eq "1" ]]; then
  echo "RUNNING IN DOCKER"
  # docker を走らせ終わったらコンテナを削除。(-rm)ホストディレクトリをコンテナにマウントする。(-v)

  DRUN="$DOCKER run --rm -v $PWD:/home --workdir /home -i"
#   DRUN_SIMPLE="$DOCKER run --rm -v $PWD:/home --workdir /home -i"
  
  #--user=biodocker
  
  # 危険！
  chmod 777 .
  
  COWSAY_IMAGE=docker/whalesay
  SRA_TOOLKIT_IMAGE=inutano/sra-toolkit
  FASTQC_IMAGE=biocontainers/fastqc:v0.11.5_cv2
  MULTIQC_IMAGE=maxulysse/multiqc
  FASTXTOOLS_IMAGE=biocontainers/fastxtools:v0.0.14_cv2
#   TRIMMOMATIC_IMAGE=fjukstad/trimmomatic
#   TRIMMOMATIC_IMAGR=comics/trimmomatic
  SALMON_IMAGE=combinelab/salmon:latest
#   SALMON_IMAGE=fjukstad/salmon
  RSCRIPT_TXIMPORT_IMAGE=fjukstad/tximport
  
  $DOCKER pull $COWSAY_IMAGE
  $DOCKER pull $SRA_TOOLKIT_IMAGE
  $DOCKER pull $FASTQC_IMAGE
  $DOCKER pull $MULTIQC_IMAGE
  $DOCKER pull $FASTXTOOLS_IMAGE
#   $DOCKER pull $TRIMMOMATIC_IMAGE
  $DOCKER pull $SALMON_IMAGE
  $DOCKER pull $RSCRIPT_TXIMPORT_IMAGE

  COWSAY="$DRUN $COWSAY_IMAGE $COWSAY"
  PREFETCH="$DRUN -v $PWD:/root/ncbi/public/sra $SRA_TOOLKIT_IMAGE $PREFETCH"
  PFASTQ_DUMP="$DRUN $SRA_TOOLKIT_IMAGE $PFASTQ_DUMP"
  FASTQ_DUMP="$DRUN $SRA_TOOLKIT_IMAGE $FASTQ_DUMP"
  FASTQC="$DRUN $FASTQC_IMAGE $FASTQC"
#   MULTIQC="$DRUN_SIMPLE $MULTIQC_IMAGE $MULTIQC"
  MULTIQC="$DRUN $MULTIQC_IMAGE $MULTIQC"
  FASTXTRIMMER="$DRUN $FASTXTOOLS_IMAGE $FASTXTRIMMER"
  FASTQQUALITYTRIMMER="$DRUN $FASTXTOOLS_IMAGE $FASTQQUALITYTRIMMER"
#   TRIMMOMATIC="$DRUN $TRIMMOMATIC_IMAGE $TRIMMOMATIC"
#   TRIMMOMATIC="$DRUN $TRIMMOMATIC_IMAGE " # fjukstad/trimmomaticのentrypointのため
  SALMON="$DRUN $SALMON_IMAGE $SALMON"
#   SALMON="$DRUN $SALMON_IMAGE"
  RSCRIPT_TXIMPORT="$DRUN $RSCRIPT_TXIMPORT_IMAGE $RSCRIPT_TXIMPORT"
  
   # docker run --rm -v $PWD:/data -v $PWD:/root/ncbi/public/sra --workdir /data -it inutano/sra-toolkit bash
else
  echo "RUNNING LOCAL"
fi

# # 十分大きなものにする。
# MAXSIZE=20G
# SRA_ROOT=$HOME/ncbi/public/sra

# # テスト用。ダウンロードするread数。全部使うときは0に
# MAX_SPOT_ID=5000000

# if [ $MAX_SPOT_ID = 0 ]; then
#   MAX_SPOT_ID=""
# else
#   $COWSAY "test mode( MAX_SPOT_ID is set)"
#   MAX_SPOT_ID="-X $MAX_SPOT_ID"
# fi

echo ${1}
cat $1

# tximport_R.Rを取ってくる。
cp $SCRIPT_DIR/tximport_R.R ./

if [[ ! -f "multiqc_report_raw_reads.html" ]]; then
  $MULTIQC -n multiqc_report_raw_reads.html .
fi

# fastq_dump
for i in `tail -n +2  $1`
do
  name=`echo $i | cut -d, -f1`
  fq=`echo $i | cut -d, -f2`
  fqname_ext="${fq##*/}"
  # echo $fqname_ext

  # ファイル名を取り出す（拡張子なし）
  basename_fq="${fqname_ext%.*.*}"
  dirname_fq=`dirname $fq`


  
  # fastqc
  if [[ ! -f "${dirname_fq}/${basename_fq}_fastqc.zip" ]]; then
    $FASTQC -t $THREADS ${dirname_fq}/${basename_fq}.fastq.gz
  fi
  
  # fastx-toolkit
  if [[ ! -f "${dirname_fq}/${basename_fq}_trimmed.fastq.gz" ]]; then
    gunzip -c ${dirname_fq}/${basename_fq}.fastq.gz | $FASTXTRIMMER -Q33 -f 1 -l 220 | $FASTQQUALITYTRIMMER -z -Q33 -t 18 -l 20 -o ${dirname_fq}/${basename_fq}_trimmed.fastq.gz
  fi

  # fastqc
  if [[ ! -f "${dirname_fq}/${basename_fq}_trimmed_fastqc.zip" ]]; then
    $FASTQC -t $THREADS ${dirname_fq}/${basename_fq}_trimmed.fastq.gz
  fi
  
done

 # multiqc
# if [[ ! -f "multiqc_report_rawfastq.html" ]]; then
#   $MULTIQC -n multiqc_report_rawfastq.html .
# fi

# download $REF_TRANSCRIPT
if [[ ! -f "$REF_TRANSCRIPT" ]]; then
  wget $BASE_REF_TRANSCRIPT/$REF_TRANSCRIPT
fi

# # download $REF_GTF
# if [[ ! -f "$REF_GTF" ]]; then
#   wget $BASE_REF_TRANSCRIPT/$REF_GTF
# fi

# instance salmon index
if [[ ! -d "$SALMON_INDEX" ]]; then
  $SALMON index --threads $THREADS --transcripts $REF_TRANSCRIPT --index $SALMON_INDEX --type quasi -k 31 --gencode
fi

for i in `tail -n +2  $1`
do
  name=`echo $i | cut -d, -f1`
  fq=`echo $i | cut -d, -f2`
  fqname_ext="${fq##*/}"
  # echo $fqname_ext

  # ファイル名を取り出す（拡張子なし）
  basename_fq="${fqname_ext%.*.*}"
  dirname_fq=`dirname $fq`
  
  if [[ ! -f "salmon_output_${basename_fq}/quant.sf" ]]; then
    mkdir salmon_output_${basename_fq}
    # libtype auto detection mode
    $SALMON quant -i $SALMON_INDEX \
    -l A \
    -r ${dirname_fq}/${basename_fq}_trimmed.fastq.gz \
    -p $THREADS \
    -o salmon_output_${basename_fq} \
#   -g $REF_GTF
  fi
done

# multiqc
if [[ ! -f "multiqc_report.html" ]]; then
  $MULTIQC -n multiqc_report.html .
fi

# download $TX2SYMBOL
if [[ ! -f "$TX2SYMBOL" ]]; then
  wget $BASE_REF_TRANSCRIPT/$TX2SYMBOL
fi

# tximport
if [[ ! -f "counttable.tsv" ]]; then
  $RSCRIPT_TXIMPORT tximport_R.R $TX2SYMBOL $EX_MATRIX_FILE
fi


if [[ "$RUNINDOCKER" -eq "1" ]]; then

  chmod 755 .

fi
