#!/bin/bash

# preprocess_homalodisca.sh: Prepare Homalodisca RNA-seq data for Xylella detection
# 
# This script filters out host reads and prepares clean reads for pathogen detection.

set -e  # Exit if something goes south

# Default parameters
HOST_GENOME=""
OUTPUT_DIR="processed"
THREADS=4

# Function to show usage
usage() {
    echo "Usage: $0 -a <R1.fastq> -b <R2.fastq> -g <host_genome.fasta> -o <output_dir> -s <sample_name>"
    echo ""
    echo "Options:"
    echo "  -a    Forward reads (R1) FASTQ file"
    echo "  -b    Reverse reads (R2) FASTQ file" 
    echo "  -g     Host genome FASTA file (Homalodisca)"
    echo "  -o     Output directory (default: processed)"
    echo "  -s     Sample name"
    echo "  -t     Number of threads (default: 4)"
    echo "  -h     Show this help"
    echo ""
    echo "Example:"
    echo "  $0 -a sample1_R1.fastq.gz -b sample1_R2.fastq.gz -g homalodisca_genome.fasta -s sample1"
    exit 1
}

# Function to log messages
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" >&2
}

# Function to check if command exists
check_command() {
    if ! command -v "$1" &> /dev/null; then
        log "ERROR: $1 is not installed or not in PATH"
        exit 1
    fi
}

# Parse command line arguments
while getopts "a:b:g:o:s:t:h" opt; do
    case $opt in
        a) R1_FILE="$OPTARG" ;;
        b) R2_FILE="$OPTARG" ;;
        g) HOST_GENOME="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        s) SAMPLE_NAME="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Check required arguments
if [[ -z "$R1_FILE" || -z "$R2_FILE" || -z "$SAMPLE_NAME" ]]; then
    log "ERROR: Missing required arguments"
    usage
fi

# Check if files exist
if [[ ! -f "$R1_FILE" ]]; then
    log "ERROR: R1 file not found: $R1_FILE"
    exit 1
fi

if [[ ! -f "$R2_FILE" ]]; then
    log "ERROR: R2 file not found: $R2_FILE"
    exit 1
fi

# Check dependencies
log "Checking dependencies..."
check_command "STAR"
check_command "samtools"

# Check if host genome filtering is requested
FILTER_HOST=false
if [[ -n "$HOST_GENOME" ]]; then
    if [[ ! -f "$HOST_GENOME" && ! -d "$HOST_GENOME" ]]; then
        log "ERROR: Host genome path not found (checked file and directory): $HOST_GENOME"
        exit 1
    fi
    FILTER_HOST=true
    # Index folder based on the name of the genome
    GENOME_DIR="${HOST_GENOME}_star_index"
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

log "Starting preprocessing for sample: $SAMPLE_NAME"
log "R1 file: $R1_FILE"
log "R2 file: $R2_FILE"
log "Output directory: $OUTPUT_DIR"

# Count input reads
log "Counting input reads..."
if [[ "$R1_FILE" == *.gz ]]; then
    TOTAL_READS=$(zcat "$R1_FILE" | wc -l | awk '{print $1/4}')
else
    TOTAL_READS=$(wc -l < "$R1_FILE" | awk '{print $1/4}')
fi
log "Input reads: $TOTAL_READS"

# Bloque de trimming y  conteo
log "Starting quality trimming with fastp..."

# Definimos los nombres de los archivos trimmeados
TRIMMED_R1="${OUTPUT_DIR}/${SAMPLE_NAME}_trimmed_R1.fastq.gz"
TRIMMED_R2="${OUTPUT_DIR}/${SAMPLE_NAME}_trimmed_R2.fastq.gz"

fastp -i "$R1_FILE" -I "$R2_FILE" \
      -o "$TRIMMED_R1" -O "$TRIMMED_R2" \
      --html "${OUTPUT_DIR}/${SAMPLE_NAME}_fastp.html" \
      --thread "$THREADS" \
      --detect_adapter_for_pe \
      2> "${OUTPUT_DIR}/${SAMPLE_NAME}_fastp.log" # Guardamos el reporte de fastp

# Contamos cuántas lecturas sobrevivieron
# Dividimos por 4 porque cada lectura en FASTQ ocupa 4 líneas
log "Counting reads after trimming..."
TRIMMED_READS=$(zcat "$TRIMMED_R1" | echo $((`wc -l`/4)))

log "Input reads: $INPUT_READS"
log "Reads after trimming: $TRIMMED_READS"

# Actualizamos las variables para que STAR use los archivos limpios
R1_FILE="$TRIMMED_R1"
R2_FILE="$TRIMMED_R2"
# -----------------------------------------

# Host genome filtering (if genome provided)
if [[ "$FILTER_HOST" == true ]]; then
    log "Filtering host reads using STAR..."
    
    # Check if STAR index directory exists, create if not
    if [[ ! -d "$GENOME_DIR" ]] || [[ ! -f "$GENOME_DIR/SAindex" ]]; then
        log "Building STAR index for host genome in $GENOME_DIR..."
        mkdir -p "$GENOME_DIR"
        STAR --runThreadN "$THREADS" \
             --runMode genomeGenerate \
             --genomeDir "$GENOME_DIR" \
             --genomeFastaFiles "$HOST_GENOME"
    else
        log "STAR index found in $GENOME_DIR, skipping index build."
    fi
    
    # Align reads to host genome using STAR
    log "Aligning reads to host genome with STAR..."
    
    STAR --runThreadN "$THREADS" \
         --genomeDir "$GENOME_DIR" \
         --readFilesIn "$R1_FILE" "$R2_FILE" \
         --readFilesCommand zcat \
         --outFileNamePrefix "$OUTPUT_DIR/${SAMPLE_NAME}_" \
         --outSAMtype SAM \
	 --outFilterMismatchNmax 10 \
         --outReadsUnmapped Fastx \
	 --outTmpDir "$OUTPUT_DIR/${SAMPLE_NAME}_STARtmp"

    # Files to move unmapped non-host reads
    NONHOST_R1="${OUTPUT_DIR}/${SAMPLE_NAME}_nonhost_R1.fastq"
    NONHOST_R2="${OUTPUT_DIR}/${SAMPLE_NAME}_nonhost_R2.fastq"

    # Define unmapped files created by default
    UNMAPPED_R1="$OUTPUT_DIR/${SAMPLE_NAME}_Unmapped.out.mate1"
    UNMAPPED_R2="$OUTPUT_DIR/${SAMPLE_NAME}_Unmapped.out.mate2"

    # Robust file handling in case nothing maps
    log "Organizing non-host reads..."
    if [[ -f "$UNMAPPED_R1" && -f "$UNMAPPED_R2" ]]; then
        mv "$UNMAPPED_R1" "$NONHOST_R1"
        mv "$UNMAPPED_R2" "$NONHOST_R2"
    else
        log "WARNING: Unmapped files not found. Creating empty files for compatibility."
        touch "$NONHOST_R1"
        touch "$NONHOST_R2"
    fi
    
    # Definimos el SAM_FILE para el conteo y limpieza posterior (nombre por defecto de STAR)
    SAM_FILE="$OUTPUT_DIR/${SAMPLE_NAME}_Aligned.out.sam"
    
    # Count non-host reads
    NONHOST_READS=$(wc -l < "$NONHOST_R1" | awk '{print $1/4}')
    HOST_PERCENTAGE=$(awk "BEGIN {printf \"%.1f\", ($TOTAL_READS - $NONHOST_READS) / $TOTAL_READS * 100}")
    
    log "Non-host reads: $NONHOST_READS"
    log "Host reads filtered: ${HOST_PERCENTAGE}%"
    
    # Clean up SAM file
    if [[ -f "$SAM_FILE" ]]; then
        log "Cleaning up large SAM file..."
        rm "$SAM_FILE"
    fi
    
    # Clean up large SAM and STAR logs
    log "Cleaning up STAR intermediate files..."
    rm -f "$OUTPUT_DIR/${SAMPLE_NAME}_Aligned.out.sam"
    rm -f "$OUTPUT_DIR/${SAMPLE_NAME}_Log.out"
    rm -f "$OUTPUT_DIR/${SAMPLE_NAME}_Log.progress.out"
    rm -f "$OUTPUT_DIR/${SAMPLE_NAME}_Log.final.out"
    rm -f "$OUTPUT_DIR/${SAMPLE_NAME}_SJ.out.tab"
    rm -rf "$OUTPUT_DIR/${SAMPLE_NAME}_STARtmp"
    
    # Set output files for downstream analysis
    FINAL_R1="$NONHOST_R1"
    FINAL_R2="$NONHOST_R2"
    
else
    log "Skipping host filtering (no host genome provided)"
    
    # Copy or link input files
    FINAL_R1="$OUTPUT_DIR/${SAMPLE_NAME}_R1.fastq"
    FINAL_R2="$OUTPUT_DIR/${SAMPLE_NAME}_R2.fastq"
    
    if [[ "$R1_FILE" == *.gz ]]; then
        log "Decompressing input files..."
        zcat "$R1_FILE" > "$FINAL_R1"
        zcat "$R2_FILE" > "$FINAL_R2"
    else
        cp "$R1_FILE" "$FINAL_R1"
        cp "$R2_FILE" "$FINAL_R2"
    fi
fi

# Basic quality statistics
log "Generating quality statistics..."

# Count final reads
FINAL_READ_COUNT=$(wc -l < "$FINAL_R1" | awk '{print $1/4}')

# Simple quality metrics
AVG_LENGTH=$(awk 'NR%4==2{sum+=length($0); count++} END{if (count>0) print int(sum/count); else print 0}' "$FINAL_R1")

log "Final read count: $FINAL_READ_COUNT"
log "Average read length: $AVG_LENGTH bp"

# Create summary file
SUMMARY_FILE="$OUTPUT_DIR/${SAMPLE_NAME}_preprocessing_summary.txt"
cat > "$SUMMARY_FILE" << EOF
Preprocessing Summary for Sample: $SAMPLE_NAME
==============================================
Input R1 file: $R1_FILE
Input R2 file: $R2_FILE
Total input reads: $TOTAL_READS
Total reads post-trimming: $TRIMMED_READS
Host genome: ${HOST_GENOME:-"Not provided"}
Host filtering: $FILTER_HOST
Final read count: $FINAL_READ_COUNT
Average read length: $AVG_LENGTH bp
Output R1: $FINAL_R1
Output R2: $FINAL_R2
Processing date: $(date)
EOF

log "Preprocessing completed successfully!"
log "Summary saved to: $SUMMARY_FILE"
log "Processed files ready for pathogen detection:"
log "  R1: $FINAL_R1"
log "  R2: $FINAL_R2"

# Optional: combine R1 and R2 into single file for downstream analysis
COMBINED_FILE="$OUTPUT_DIR/${SAMPLE_NAME}_combined.fastq"
log "Creating combined FASTQ file..."
cat "$FINAL_R1" "$FINAL_R2" > "$COMBINED_FILE"
log "Combined file: $COMBINED_FILE"
