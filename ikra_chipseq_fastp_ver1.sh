#!/bin/bash
set -xe

# オプション関連ここから
# 大部分は http://dojineko.hateblo.jp/entry/2016/06/30/225113 から引用させていただきました。

# 変数 EX_MATRIX_FILE, REF_SPECIES はここで定義
# if [[ $IF_TEST = true ]]; then でテストモード用の実行が可能

# 今まで$1 = EX_MATRIX_FILEだったのを変更している
# 以降の$1をEX_MATRIX_FILEで置き換える必要がある？(必要なら修正お願いします...)
#ikraのchip-seqバージョンに変更する。hisat2をbowtie2に書き換える。
#Pair-endの場合はhttps://github.com/yuifu/ngsdat2_epigenome_chipseq/blob/master/atacseq.mdのATAC-seqのやり方を参考にした。
#trim-galoreを代わりにfastpに変更。
#PEで使うpicardはなぜかdockerで使えなかったので,自分PCにインストールする必要あり。
#というかペアエンドではwithout docker(-w)じゃないとなぜかうまくいかない。

PROGNAME="$( basename $0 )"

VERSION="v2.0.1"

cat << "EOF"

# Example command in this experiment are blow.

# $cd chip_seq_test && bash ../ikra_chipseq_fastp_ver.sh chip_seq_test.csv mouse --align bowtie2 --threads 10 --fastq -w

#
#MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
#MMMMMMMMMMMMMMMMMY"=!``.``   _?TWMMMMMMMMMMMMMMMMM
#MMMMMMMMMMMMMB=.`.`.`....-,      ..?TMMMMMMMMMMMMM
#MMMMMMMMMM#^...`....7=_`.`.`.`.... .7G-WMMMMMMMMMM
#MMMMMMMM@~...`...Z!...`...`..."!.`     7JTMMMMMMMM
#MMMMMM#^...`...('.......`...P  .MMMNe.   ?JHMMMMMM
#MMMMMF........(~....`..`...d  JMMMMMMN  ,  4dMMMMM
#MMMMF........(]...........`W. MMMMMMMM~.f   j?MMMM
#MMMF.........(r............(e (MMMMMM^`.t`   jdMMM
#MM#_........(CW,............Jm,..?!` .,^``.   bMMM
#MMF........(C::?a-............(TTwz7"<.... ``.v,MM
#MM~......._P:::::?5J,............-_..`..`..Uy"`.MM
#M#....~...($:::::::+JTYY4uJ(-(--(((-.-J7T>~\````MM
#MN..__.~.~(P:::::+Y<>>>>>>>>>>>>>>>;;;;dG-J.`.``MM
#MM-.  ..~._E:::+Wz?>??+++>>>>>>>>>;>>++Z;?l`.`..MM
#MM].  ~....?p:+D1kjMNvMMP>>>>>>>>+z1>>+Cc;+t```.MM
#MMN__ _~.~~.dpd?>vogmvMMPj&&xu&xux>j&gA+J;+J-. dMM
#MMMb.  ...~.~(E>?>wMMIMMQMM8<dMMMEdMMMMM#C>j\`.MMM
#MMMMp.  _..~..j+?>jMMIMMMMNe>dM#>>MMN+dM#;;J`.MMMM
#MMMMMp..  ..~..4zIvMMIMM$?MM8dMB>>+TMMWMB;j3.MMMMM
#MMMMMMN,_  `.~..?x?z??>>?>>>>>>>>>?1d+;;+J1dMMMMMM
#MMMMMMMMm-_.  ~.~_7Gx>>?>>>>?>>>>>;>>j+Z<(MMMMMMMM
#MMMMMMMMMMm,........_?Tu&x+>>>>>+&&v"!.JMMMMMMMMMM
#MMMMMMMMMMMMNa,~............___....-(MMMMMMMMMMMMM
#MMMMMMMMMMMMMMMMNg..-........-((+MMMMMMMMMMMMMMMMM
#MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
#



EOF

# Usage
function usage() {
  cat << EOS >&2        
ikra ${VERSION} -RNAseq pipeline centered on Salmon-
Usage: ${PROGNAME} experiment_table.csv species [--test, --fastq, --help, --without-docker, --udocker, --protein-coding] [--threads [VALUE]][--output [VALUE]][--suffix_PE_1 [VALUE]][--suffix_PE_2 [VALUE]]
  args
    1.experiment matrix(csv)
    2.reference(human or mouse)
Options:
  --test  test mode(MAX_SPOT_ID=100000). (dafault : False)
  --fastq use fastq files instead of SRRid. The extension must be foo.fastq.gz (default : False)
  -u, --udocker
  -w, --without-docker
  -pc, --protein-coding use protein coding transcripts instead of comprehensive transcripts. (defalut : True)
  -ct, --comprehensive-transcripts use comprehensive transcripts instead of protein coding transcripts. (default : False) 
  -t, --threads
  -o, --output  output file. (default : output.tsv)  
  -l, --log  log file. (default : ikra.log)
  -a, --align carry out mapping onto a reference genome. bowtie2 or star (default : None)
  -g, --gencode specify the version of gencode. (defalut : Mouse=26, Human=37)
  -s1, --suffix_PE_1    suffix for PE fastq files. (default : _1.fastq.gz)
  -s2, --suffix_PE_2    suffix for PE fastq files. (default : _2.fastq.gz)
  -h, --help    Show usage.
  -v, --version Show version.
  -r, --remove-intermediates Remove intermediate files

Citation :
Hiraoka, Y., Yamada, K., Yamasaki, R., Kawasaki, Y., Kitabatake, R., Matsumoto, Y., Ishikawa, K., Umezu, Y., Hirose, H., & Yasumizu, Y. (2021). ikra v2.0: RNAseq pipeline centered on Salmon. https://doi.org/10.5281/zenodo.4718200

Github repo : https://github.com/yyoshiaki/ikra
EOS
  exit 1
}


# version
function version() {
  cat << EOS >&2
ikra ${VERSION} -RNAseqじゃなくてChip-seq pipeline-
EOS
  exit 1
}

# デフォルト値を先に定義しておく
RUNINDOCKER=1
DOCKER=docker
THREADS=1
IF_TEST=false
IF_FASTQ=false
IF_PC=True
SUFFIX_PE_1=_1.fastq.gz
SUFFIX_PE_2=_2.fastq.gz
OUTPUT_FILE=output.tsv
LOG_FILE=ikra.log
MAPPING_TOOL=None
IF_REMOVE_INTERMEDIATES=false
M_GEN_VER=26
H_GEN_VER=37

# オプションをパース
PARAM=()
for opt in "$@"; do
    case "${opt}" in
        #　モード選択など引数の無いオプションの場合
        '--test' )
            IF_TEST=true; shift
            ;;
        '--fastq' )
            IF_FASTQ=true; shift
            ;;
        '-pc'|'--protein-coding' )
            IF_PC=true; shift
            ;;
        '-ct'|'--comprehensive-transcripts' )
            IF_PC=true; shift
            ;;
        '-u'|'--udocker' )
            DOCKER=udocker; shift
            ;;
        '-w'|'--without-docker' )
            RUNINDOCKER=0; shift
            ;;
        #　引数が任意の場合
        '-t'|'--threads' )
            THREADS=4; shift
            if [[ -n "$1" ]] && [[ ! "$1" =~ ^-+ ]]; then
                THREADS="$1"; shift
            fi
            ;;
        '-s1'|'--suffix_PE_1' )
            if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
                echo "$PROGNAME: option requires an argument -- $1" 1>&2
                exit 1
            fi
            SUFFIX_PE_1="$2"
            shift 2
            ;;

        '-s2'|'--suffix_PE_2' )
            if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
                echo "$PROGNAME: option requires an argument -- $1" 1>&2
                exit 1
            fi
            SUFFIX_PE_2="$2"
            shift 2
            ;;

        '-o'|'--output' )
            if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
                echo "$PROGNAME: option requires an argument -- $1" 1>&2
                exit 1
            fi
              OUTPUT_FILE="$2"
              shift 2
              ;;
        '-l'|'--log' )
            if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
                echo "$PROGNAME: option requires an argument -- $1" 1>&2
                exit 1
            fi
              LOG_FILE="$2"
              shift 2
                ;;

        '-a'|'--align' )
            if [[ "$2" == "bowtie2" ]]; then
                MAPPING_TOOL=BOWTIE2
            elif [[ "$2" == "star" ]]; then
                MAPPING_TOOL=STAR
            elif [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
                echo "$PROGNAME: option requires an argument -- $1" 1>&2
                exit 1
            fi
              shift 2
                ;;

        '-g'|'--gencode' )
            if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
                echo "${PROGNAME}: option requires an argument -- $( echo $1 | sed 's/^-*//' )" 1>&2
                exit 1
            fi
            H_GEN_VER="$2"
            M_GEN_VER="$2"
            shift 2
            ;;

        '-h' | '--help' )
            usage
            ;;
        '-v' | '--version' )
            version
            ;;
        '-r' | '--remove' )
            IF_REMOVE_INTERMEDIATES=true ; shift
            ;;
        '--' | '-' )
            shift
            PARAM+=( "$@" )
            break
            ;;
        -* )
            echo "${PROGNAME}: illegal option -- '$( echo $1 | sed 's/^-*//' )'" 1>&2
            exit 1
            ;;
        * )
            if [[ -n "$1" ]] && [[ ! "$1" =~ ^-+ ]]; then
                PARAM+=( "$1" ); shift
            fi
            ;;
    esac
done

# オプション無しの値を使う場合はここで処理する
EX_MATRIX_FILE="${PARAM}"; PARAM=("${PARAM[@]:1}")
REF_SPECIES="${PARAM}"; PARAM=("${PARAM[@]:1}")

[[ -z "${EX_MATRIX_FILE}" ]] && usage
[[ -z "${REF_SPECIES}" ]] && usage

# 規定外のオプションがある場合にはusageを表示
if [[ -n "${PARAM[@]}" ]]; then
    usage
fi


cat << EOS | tee -a ${LOG_FILE}
ikra ${VERSION} -RNAseq pipeline centered on Salmon-
EOS

date >> ${LOG_FILE}
pwd >> ${LOG_FILE}
whoami >> ${LOG_FILE}
uname -n >> ${LOG_FILE}

# 結果を表示(オプションテスト用)
cat << EOS | column -t | tee -a ${LOG_FILE}
EX_MATRIX_FILE ${EX_MATRIX_FILE}
REF_SPECIES ${REF_SPECIES}
RUNINDOCKER ${RUNINDOCKER}
DOCKER ${DOCKER}
THREADS ${THREADS}
IF_TEST ${IF_TEST:-false}
IF_FASTQ ${IF_FASTQ:-false}
IF_PC ${IF_PC:-false}
IF_REMOVE_INTERMEDIATES ${IF_REMOVE_INTERMEDIATES:-false}
OUTPUT_FILE ${OUTPUT_FILE}
MAPPING_TOOL ${MAPPING_TOOL}
M_GEN_VER ${M_GEN_VER}
H_GEN_VER ${H_GEN_VER}
LOG_FILE ${LOG_FILE}
EOS

set -u

#　オプション関連ここまで

# 実験テーブル.csv

# 十分大きなものにする。
MAXSIZE=20G
SRA_ROOT=$HOME/ncbi/public/sra

SCRIPT_DIR=$(cd $(dirname $0); pwd)

if [[ $REF_SPECIES = mouse ]]; then
  BASE_REF_TRANSCRIPT=ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M${M_GEN_VER}
  REF_TRANSCRIPT=gencode.vM${M_GEN_VER}.transcripts.fa.gz
  if [ $IF_PC = false ]; then
    REF_TRANSCRIPT=gencode.vM${M_GEN_VER}.transcripts.fa.gz
  else
    REF_TRANSCRIPT=gencode.vM${M_GEN_VER}.pc_transcripts.fa.gz
  fi
  SALMON_INDEX=salmon_index_mouse
#   REF_GTF=gencode.vM${M_GEN_VER}.annotation.gtf.gz
  TX2SYMBOL=gencode.vM${M_GEN_VER}.metadata.MGI.gz

elif [[ $REF_SPECIES = human ]]; then
  BASE_REF_TRANSCRIPT=ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${H_GEN_VER}
  # REF_TRANSCRIPT=gencode.v${H_GEN_VER}.pc_transcripts.fa.gz

  if [ $IF_PC = false ]; then
    REF_TRANSCRIPT=gencode.v${H_GEN_VER}.transcripts.fa.gz
  else
    REF_TRANSCRIPT=gencode.v${H_GEN_VER}.pc_transcripts.fa.gz
  fi

  SALMON_INDEX=salmon_index_human
#   REF_GTF=gencode.v${H_GEN_VER}.annotation.gtf.gz
  TX2SYMBOL=gencode.v${H_GEN_VER}.metadata.HGNC.gz
else
  echo No reference speice!
  exit
fi

COWSAY=cowsay
FASTQ_DUMP=fastq-dump
FASTERQ_DUMP=fasterq-dump
FASTQC=fastqc
MULTIQC=multiqc
# TRIMMOMATIC=trimmomatic
TRIMGALORE=trim_galore
FASTP=fastp
#HISAT2=hisat2
BOWTIE2=bowtie2
PICARD=picard
SAMTOOLS=samtools
#STAR_MAPPING=STAR
#RSEM=rsem
SAMBAMBA=sambamba
BAMCOVERAGE=bamCoverage
SALMON=salmon
RSCRIPT_TXIMPORT=Rscript
WGET=wget
PIGZ=pigz
TAR=tar

if [[ "$RUNINDOCKER" -eq "1" ]]; then
  echo "RUNNING IN DOCKER"
  # docker を走らせ終わったらコンテナを削除。(-rm)ホストディレクトリをコンテナにマウントする。(-v)

  if [[ $DOCKER = docker ]]; then
    DRUN="$DOCKER run  -u `id -u`:`id -g` --rm -v $PWD:/home -e HOME=/home --workdir /home "
  elif [[ $DOCKER = udocker ]]; then
    DRUN="$DOCKER run --rm -v $PWD:/home --workdir /home "
  fi

  SCRIPT_DIR=`dirname "$0"`
  #--user=biodocker

  # 危険！
  # chmod 777 .

  COWSAY_IMAGE=docker/whalesay
  SRA_TOOLKIT_IMAGE=quay.io/biocontainers/sra-tools:2.10.9--pl526haddd2b5_0
  FASTQC_IMAGE=biocontainers/fastqc:v0.11.9_cv8
  MULTIQC_IMAGE=ewels/multiqc:latest
  TRIMGALORE_IMAGE=quay.io/biocontainers/trim-galore:0.6.7--hdfd78af_0
  FASTP_IMAGE=biocontainers/fastp:v0.20.1_cv1
  #HISAT2_IMAGE=quay.io/biocontainers/hisat2:2.2.1--h1b792b2_3
  BOWTIE2_IMAGE=biocontainers/bowtie2:v2.4.1_cv1
  PICARD_IMAGE=broadinstitute/picard:3.1.0
  SAMTOOLS_IMAGE=biocontainers/samtools:v1.7.0_cv4
  #STAR_IMAGE=quay.io/biocontainers/star:2.7.9a--h9ee0642_0
  #RSEM_IMAGE=biocontainers/rsem:v1.3.1dfsg-1-deb_cv1
  SAMBAMBA_IMAGE=quay.io/biocontainers/sambamba:0.8.0--h984e79f_0
  SALMON_IMAGE=combinelab/salmon:1.4.0
  RSCRIPT_TXIMPORT_IMAGE=fjukstad/tximport
  WGET_IMAGE=fjukstad/tximport
  PIGZ_IMAGE=genevera/docker-pigz
  TAR_IMAGE=fjukstad/tximport
  BAMCOVERAGE_IMAGE=quay.io/biocontainers/deeptools:3.5.1--py_0

  $DOCKER pull $COWSAY_IMAGE
  $DOCKER pull $SRA_TOOLKIT_IMAGE
  $DOCKER pull $FASTQC_IMAGE
  $DOCKER pull $MULTIQC_IMAGE
  $DOCKER pull $TRIMGALORE_IMAGE
  $DOCKER pull $FASTP_IMAGE
  #$DOCKER pull $HISAT2_IMAGE
  $DOCKER pull $BOWTIE2_IMAGE
  $DOCKER pull $PICARD_IMAGE
  $DOCKER pull $SAMTOOLS_IMAGE
  #$DOCKER pull $STAR_IMAGE
  #$DOCKER pull $RSEM_IMAGE
  $DOCKER pull $SAMBAMBA_IMAGE
  $DOCKER pull $BAMCOVERAGE_IMAGE
  $DOCKER pull $SALMON_IMAGE
  $DOCKER pull $RSCRIPT_TXIMPORT_IMAGE
  $DOCKER pull $PIGZ_IMAGE
  $DOCKER pull $TAR_IMAGE

  COWSAY="$DRUN $COWSAY_IMAGE $COWSAY"
  FASTQ_DUMP="$FASTQ_DUMP"
  FASTQC="$DRUN $FASTQC_IMAGE $FASTQC" 
  FASTQ_DUMP="$FASTQ_DUMP"
  FASTERQ_DUMP="$FASTERQ_DUMP"
  MULTIQC="$DRUN $MULTIQC_IMAGE $MULTIQC"
  TRIMGALORE="$DRUN $TRIMGALORE_IMAGE $TRIMGALORE"
  FASTP="$DRUN $FASTP_IMAGE $FASTP"
  #HISAT2="$DRUN $HISAT2_IMAGE $HISAT2"
  BOWTIE2="$DRUN $BOWTIE2_IMAGE $BOWTIE2"
  PICARD="$DRUN $PICARD_IMAGE"
  SAMTOOLS="$DRUN $SAMTOOLS_IMAGE $SAMTOOLS"
  #STAR_MAPPING="$DRUN $STAR_IMAGE $STAR_MAPPING"
  #RSEM="$DRUN $RSEM_IMAGE $RSEM"
  SAMBAMBA="$DRUN $SAMBAMBA_IMAGE $SAMBAMBA"
  BAMCOVERAGE="$DRUN $BAMCOVERAGE_IMAGE $BAMCOVERAGE"
  SALMON="$DRUN $SALMON_IMAGE $SALMON"
  RSCRIPT_TXIMPORT="$DRUN $RSCRIPT_TXIMPORT_IMAGE $RSCRIPT_TXIMPORT"
  WGET="$DRUN $WGET_IMAGE $WGET"
  PIGZ="$DRUN $PIGZ_IMAGE"
  TAR="$DRUN $TAR_IMAGE $TAR"

   # docker run --rm -v $PWD:/data -v $PWD:/root/ncbi/public/sra --workdir /data -it inutano/sra-toolkit bash
else
  echo "RUNNING LOCAL"
fi


# if [ $MAX_SPOT_ID = 0 ]; then
if [ $IF_TEST = true ]; then
  $COWSAY "test mode( MAX_SPOT_ID is set)"
  MAX_SPOT_ID="-X 100000"
else
  MAX_SPOT_ID=""
fi

echo $EX_MATRIX_FILE
cat $EX_MATRIX_FILE

# tximport
if [[  -f "tximport_R.R" ]]; then
  rm tximport_R.R
fi

# # tximport_R.Rを取ってくる。
# cp $SCRIPT_DIR/tximport_R.R ./

# 2019/06/09 devv1.3 tximport_R.Rを埋め込み

cat << 'EOF' > tximport_R.R
#! /usr/bin/Rscript
library(tximport)
library(readr)
library(stringr)
# Rscript tximport_R.R gencode.vM19.metadata.MGI.gz Illumina_PE_SRR.csv output.tsv
args1 = commandArgs(trailingOnly=TRUE)[1]
args2 = commandArgs(trailingOnly=TRUE)[2]
args3 = commandArgs(trailingOnly=TRUE)[3]
tx2knownGene <- read_delim(args1, '\t', col_names = c('TXNAME', 'GENEID'))
exp.table <- read.csv(args2, row.names=NULL)
files.raw <- exp.table[,2]
# files.raw <- c("SE/test/ttt30.fq.gz", "SE/test/ttt2.fq.gz")
files.raw <- gsub(".gz$", "", files.raw)
files.raw <- gsub(".fastq$", "", files.raw)
files.raw <- gsub(".fq$", "", files.raw)
split.vec <- sapply(files.raw, basename)
# print(paste(c("salmon_output_") , split.vec, c("/quant.sf"), sep=''))
# files <- paste(c("salmon_output_") , exp.table[,2], c("/quant.sf"), sep='')
files <- paste(c("salmon_output_") , split.vec, c("/quant.sf"), sep='')
names(files) <- exp.table[,1]
print(files)
# txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2knownGene)
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2knownGene, countsFromAbundance="scaledTPM")
write.table(txi.salmon$counts, file=args3, sep="\t",col.names=NA,row.names=T,quote=F,append=F)
write.table(exp.table[-c(2,3)], file="designtable.csv",row.names=F,quote=F,append=F)
EOF

if [ $IF_FASTQ = false ]; then
# fastq_dump
for i in `tail -n +2  $EX_MATRIX_FILE | tr -d '\r'`
do
name=`echo $i | cut -d, -f1`
SRR=`echo $i | cut -d, -f2`
LAYOUT=`echo $i | cut -d, -f3`
# ADAPTER=`echo $i | cut -d, -f4`

<<COMMENTOUT
There is no -N|--minSpotId and no -X|--maxSpotId option.
fasterq-dump version 2.9.1 processes always the whole accession,
although it may support partial access in future versions.
ということで条件分岐させる。
COMMENTOUT

# fasterq_dump
  # SE
  if [ $LAYOUT = SE ]; then
    # fastq_dump
    if [[ ! -f "$SRR.fastq.gz" ]]; then
      if [[ $MAX_SPOT_ID == "" ]]; then
        $FASTERQ_DUMP $SRR --threads $THREADS --force -p
        # gzip $SRR.fastq
        $PIGZ $SRR.fastq
      else
        $FASTQ_DUMP $SRR $MAX_SPOT_ID --gzip
      fi
    fi

    # fastqc
    if [[ ! -f "${SRR}_fastqc.zip" ]]; then
      $FASTQC -t $THREADS ${SRR}.fastq.gz
    fi

  # PE
  else
    # fastq_dump
    if [[ ! -f "${SRR}_1.fastq.gz" ]]; then
      if [[ $MAX_SPOT_ID == "" ]]; then
        $FASTERQ_DUMP $SRR --split-files --threads $THREADS --force -p
        # gzip ${SRR}_1.fastq
        # gzip ${SRR}_2.fastq
        $PIGZ ${SRR}_1.fastq
        $PIGZ ${SRR}_2.fastq
      else
        $FASTQ_DUMP $SRR $MAX_SPOT_ID --gzip --split-files
      fi
    fi

    # fastqc
    if [[ ! -f "${SRR}_1_fastqc.zip" ]]; then
      $FASTQC -t $THREADS ${SRR}${SUFFIX_PE_1}
      $FASTQC -t $THREADS ${SRR}${SUFFIX_PE_2}
    fi
  fi
done
fi

if [[ ! -f "multiqc_report_raw_reads.html" ]]; then
  $MULTIQC -n multiqc_report_raw_reads.html .
fi

# determin threads for trim galore.
# the sweet spot for TG is 4
if [ $THREADS -gt 4 ] ; then
  THREADS_TRIMGALORE=4
else
  THREADS_TRIMGALORE=$THREADS
fi

################################ FASTP ################################
for i in `tail -n +2  $EX_MATRIX_FILE | tr -d '\r'`
do
  if [ $IF_FASTQ = false ]; then
    # fasterq_dump
    name=`echo $i | cut -d, -f1`
    SRR=`echo $i | cut -d, -f2`
    LAYOUT=`echo $i | cut -d, -f3`
    dirname_fq="./"
  else
    name=`echo $i | cut -d, -f1`
    fq=`echo $i | cut -d, -f2`
    LAYOUT=`echo $i | cut -d, -f3`
    fqname_ext="${fq##*/}"
    # echo $fqname_ext

    # ファイル名を取り出す（拡張子なし）
    # basename_fq="${fqname_ext%.*.*}"
    basename_fq=${fqname_ext}
    dirname_fq=`dirname $fq`
    dirname_fq=${dirname_fq}/
    SRR=${basename_fq}
  fi

  # trim_galoreあらためfastp
  # SE
  if [ $LAYOUT = SE ]; then
    if [  -f "${dirname_fq}${SRR}.fq" ] && [ ! -f "${dirname_fq}${SRR}.fastq.gz" ]; then
      $PIGZ ${dirname_fq}${SRR}.fq
      ln -s ${SRR}.fq.gz ${dirname_fq}${SRR}.fastq.gz
    fi
    if [ -f "${dirname_fq}${SRR}.fastq" ] && [ ! -f "${dirname_fq}${SRR}.fastq.gz" ]; then
      $PIGZ ${dirname_fq}${SRR}.fastq
    fi
    if [ -f "${dirname_fq}${SRR}.fq.gz" ] && [ ! -f "${dirname_fq}${SRR}.fastq.gz" ]; then
      ln -s ${SRR}.fq.gz ${dirname_fq}${SRR}.fastq.gz
    fi
  #fastp -i ${dirname_fq}${SRR}.fastq.gz -o ${dirname_fq}${SRR}_trimmed.fq.gz --html {dirname_fq}${SRR}_trimmed.fq.gz.fastp.html
    if [[ ! -f "${dirname_fq}${SRR}_trimmed.fq.gz" ]]; then
      $FASTP -i ${dirname_fq}${SRR}.fastq.gz -o ${dirname_fq}${SRR}_trimmed.fq.gz --html ${dirname_fq}${SRR}_trimmed.fq.gz.fastp.html
    fi

    # fastqc
#    if [[ ! -f "${dirname_fq}${SRR}_trimmed_fastqc.zip" ]]; then
#      $FASTQC -t $THREADS ${dirname_fq}${SRR}_trimmed.fq.gz
#    fi

  # PE
  else
    if [ -f "${dirname_fq}${SRR}_1.fq" ] && [ ! -f "${dirname_fq}${SRR}_1.fastq.gz" ]; then
      ${PIGZ} ${dirname_fq}${SRR}_1.fq
      ${PIGZ} ${dirname_fq}${SRR}_2.fq
      ln -s ${SRR}_1.fq.gz ${dirname_fq}${SRR}_1.fastq.gz
      ln -s ${SRR}_2.fq.gz ${dirname_fq}${SRR}_2.fastq.gz
    fi
    if [ -f "${dirname_fq}${SRR}_1.fastq" ] && [ ! -f "${dirname_fq}${SRR}_1.fastq.gz" ]; then
      $PIGZ ${dirname_fq}${SRR}_1.fastq
      $PIGZ ${dirname_fq}${SRR}_2.fastq
    fi
    if [  -f "${dirname_fq}${SRR}.fq.gz" ] && [ ! -f "${dirname_fq}${SRR}_1.fastq.gz" ]; then
      cp ${dirname_fq}${SRR}.fq.gz ${dirname_fq}${SRR}.fastq.gz
    fi
    if [ -f "${dirname_fq}${SRR}${SUFFIX_PE_1}" ] && [ ! -f "${dirname_fq}${SRR}_1.fastq.gz" ]; then
      ln -s ${SRR}${SUFFIX_PE_1} ${dirname_fq}${SRR}_1.fastq.gz
      ln -s ${SRR}${SUFFIX_PE_2} ${dirname_fq}${SRR}_2.fastq.gz
    fi

    # trim-galoreじゃなくてfastp
    #FASTP --in1 data/fastq_atacseq/SRR891269_1.fastq.gz --in2 data/fastq_atacseq/SRR891269_2.fastq.gz \
    #--out1 fastp/SRR891269_1.trimmed.fastq.gz --out2 fastp/SRR891269_2.trimmed.fastq.gz \
    #--html fastp/SRR891269.fastp.html
 
    if [[ ! -f "${dirname_fq}${SRR}_1_val_1.fq.gz" ]]; then
      $FASTP --in1 ${dirname_fq}${SRR}_1.fastq.gz --in2 ${dirname_fq}${SRR}_2.fastq.gz \
      --out1 ${dirname_fq}${SRR}_1_val_1.fq.gz --out2 ${dirname_fq}${SRR}_2_val_2.fq.gz \
      --html $THREADS ${dirname_fq}${SRR}.fastp.html
    fi

    # fastqc
    #if [[ ! -f "${dirname_fq}${SRR}_1_val_1_fastqc.zip" ]]; then
     #$FASTQC -t $THREADS ${dirname_fq}${SRR}_1_val_1.fq.gz
    #$FASTQC -t $THREADS ${dirname_fq}${SRR}_2_val_2.fq.gz
    #fi
  fi
done

# download $REF_TRANSCRIPT
#if [[ ! -f "$REF_TRANSCRIPT" ]]; then
#  $WGET $BASE_REF_TRANSCRIPT/$REF_TRANSCRIPT
#fi

# # download $REF_GTF
# if [[ ! -f "$REF_GTF" ]]; then
#   wget $BASE_REF_TRANSCRIPT/$REF_GTF
# fi

################################ BOWTIE2 ################################

if [[ $MAPPING_TOOL = BOWTIE2 ]]; then

  # download reference genome index
  if [[ $REF_SPECIES = mouse ]]; then
    BASE_REF_GENOME=ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes
    REF_GENOME=mm10.zip
    SPECIES_NAME=mm10
  elif [[ $REF_SPECIES = human ]]; then
    BASE_REF_GENOME=ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines
    REF_GENOME=GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz
    SPECIES_NAME=hg38
  else
    echo No reference genome!
    exit
  fi

#  if [[ ! -f "$REF_GENOME" ]]; then
#   $WGET $BASE_REF_GENOME/$REF_GENOME
#    if [[ $REF_SPECIES = human ]]; then
#      tar xvzf $REF_GENOME
#    else
#      unzip $REF_GENOME
#    fi
#  else
#    if [[ $REF_SPECIES = human ]]; then
#      tar xvzf $REF_GENOME
#    else
#      tar xvzf $REF_GENOME
#    fi
#  fi
fi

for i in `tail -n +2  $EX_MATRIX_FILE | tr -d '\r'`
do
  echo "Processing: $i"
    if [ $IF_FASTQ = false ]; then
      # fasterq_dump
      name=`echo $i | cut -d, -f1`
      SRR=`echo $i | cut -d, -f2`
      LAYOUT=`echo $i | cut -d, -f3`
      dirname_fq=""
    else
      name=`echo $i | cut -d, -f1`
      fq=`echo $i | cut -d, -f2`
      LAYOUT=`echo $i | cut -d, -f3`
      fqname_ext="${fq##*/}"
      # echo $fqname_ext

      # ファイル名を取り出す（拡張子なし）
      # basename_fq="${fqname_ext%.*.*}"
      basename_fq=${fqname_ext}
      dirname_fq=`dirname $fq`
      dirname_fq=${dirname_fq}/
      SRR=${basename_fq}
    fi
    #SE
    if [ $LAYOUT = SE ]; then
      if [[ ! -f "bowtie2_output_${SRR}/${SRR}.bw" ]]; then
        mkdir bowtie2_output_${SRR}
        # libtype auto detection mode
        bowtie2\
        -p $THREADS \
        -x bowtie2_index/mm10 \
        -U ${dirname_fq}${SRR}_trimmed.fq.gz \
        > bowtie2_output_${SRR}/${SRR}.trim.sam
        #samtools
        $SAMTOOLS view -bhS -F 0x4 -q 42 bowtie2_output_${SRR}/${SRR}.trim.sam | \
        $SAMTOOLS sort -T bowtie2_output_${SRR}/${SRR}.trim - > bowtie2_output_${SRR}/${SRR}.trim.uniq.bam
        $SAMBAMBA index bowtie2_output_${SRR}/${SRR}.trim.uniq.bam -t $THREADS
        $BAMCOVERAGE -b bowtie2_output_${SRR}/${SRR}.trim.uniq.bam -o bowtie2_output_${SRR}/${SRR}.bw --numberOfProcessors=max
        rm bowtie2_output_${SRR}/${SRR}.trim.sam
      fi

    # PE
    else
      if [[ ! -f "bowtie2_output_${SRR}/${SRR}.bw" ]]; then
        mkdir bowtie2_output_${SRR}
        # libtype auto detection mode
        bowtie2\
        -p $THREADS \
        --no-mixed \
        --no-discordant \
        -X 2000 \
        -x bowtie2_index/mm10 \
        -1 ${dirname_fq}${SRR}_1_val_1.fq.gz \
        -2 ${dirname_fq}${SRR}_2_val_2.fq.gz \
        > bowtie2_output_${SRR}/${SRR}.trim.sam
        
        #samtools view&sort&index
        $SAMTOOLS view -bhS -F 0x4 bowtie2_output_${SRR}/${SRR}.trim.sam | \
        $SAMTOOLS sort -T bowtie2_output_${SRR}/${SRR}.trim - > bowtie2_output_${SRR}/${SRR}.trim.bam
        $SAMBAMBA  index bowtie2_output_${SRR}/${SRR}.trim.bam
        #適切なペアとしてマップされたリード（"read mapped in proper pair"）のみを抽出
        $SAMTOOLS view -f 0x2 -bh bowtie2_output_${SRR}/${SRR}.trim.bam \
        > bowtie2_output_${SRR}/${SRR}.trim.proper_pairs.bam
        $SAMTOOLS index bowtie2_output_${SRR}/${SRR}.trim.proper_pairs.bam
        cd bowtie2_output_${SRR}
        PICARD MarkDuplicates \
        I=${SRR}.trim.proper_pairs.bam \
        O=${SRR}.trim.proper_pairs.rmdup.bam \
        M=${SRR}.trim.proper_pairs.rmdup.bam.log.txt \
        REMOVE_DUPLICATES=true
        cd ..
        $SAMTOOLS index bowtie2_output_${SRR}/${SRR}.trim.proper_pairs.rmdup.bam
        $BAMCOVERAGE -b bowtie2_output_${SRR}/${SRR}.trim.proper_pairs.rmdup.bam \
        -o bowtie2_output_${SRR}/${SRR}.trim.proper_pairs.rmdup.bw --numberOfProcessors=max
        rm bowtie2_output_${SRR}/${SRR}.trim.sam
      fi
     # if [[ ! -f "picard_output_${SRR}/${SRR}.bw" ]]; then
        #Picardによるフラグメント長（インサートサイズ）の分布の可視化
        #mkdir picard_output_${SRR}
        #$PICARD CollectInsertSizeMetrics \
        #INPUT=bowtie2_output_${SRR}/${SRR}.trim.proper_pairs.rmdup.bam \
        #OUTPUT=picard_output_${SRR}/${SRR}_insert_size_metrics.txt \
        #HISTOGRAM_FILE=picard_output_${SRR}/${SRR}_hist.pdf \
        #MINIMUM_PCT=0
    #  fi
    fi

    #sambamba
   # if [[ ! -f "bowtie2_output_${SRR}/${SRR}.bam" ]]; then
   #   $SAMBAMBA view -S bowtie2_output_${SRR}/${SRR}.sam -f bam -o bowtie2_output_${SRR}/${SRR}.bam -t $THREADS
   #   $SAMBAMBA sort bowtie2_output_${SRR}/${SRR}.bam -m
   #   $SAMBAMBA index bowtie2_output_${SRR}/${SRR}.sorted.bam -t $THREADS
   #   $BAMCOVERAGE -b bowtie2_output_${SRR}/${SRR}.sorted.bam -o bowtie2_output_${SRR}/${SRR}.bw --numberOfProcessors=max
   # fi
    
   # rm bowtie2_output_${SRR}/${SRR}.trim.sam
    
done



################################ bowtie2終わり　################################

# if [[ "$RUNINDOCKER" -eq "1" ]]; then
#
#   chmod 755 .
#
# fi

cat << EOS | tee -a ${LOG_FILE}
RUN : success!
EOS
