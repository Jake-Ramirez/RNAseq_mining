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
    echo "Usage: $0 -r1 <R1.fastq> -r2 <R2.fastq> -g <host_genome.fasta> -o <output_dir> -s <sample_name>"
    echo ""
    echo "Options:"
    echo "  -r1    Forward reads (R1) FASTQ file"
    echo "  -r2    Reverse reads (R2) FASTQ file" 
    echo "  -g     Host genome FASTA file (Homalodisca)"
    echo "  -o     Output directory (default: processed)"
    echo "  -s     Sample name"
    echo "  -t     Number of threads (default: 4)"
    echo "  -h     Show this help"
    echo ""
    echo "Example:"
    echo "  $0 -r1 sample1_R1.fastq.gz -r2 sample1_R2.fastq.gz -g homalodisca_genome.fasta -s sample1"
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
while getopts "r1:r2:g:o:s:t:h" opt; do
    case $opt in
        r1) R1_FILE="$OPTARG" ;;
        r2) R2_FILE="$OPTARG" ;;
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
'''
# Check dependencies
log "Checking dependencies..."
check_command "STAR"
check_command "samtools"

# Check if host genome filtering is requested
FILTER_HOST=false
if [[ -n "$HOST_GENOME" ]]; then
    if [[ ! -f "$HOST_GENOME" ]]; then
        log "ERROR: Host genome file not found: $HOST_GENOME"
        exit 1
    fi
    FILTER_HOST=true
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

# Host genome filtering (if genome provided)
if [[ "$FILTER_HOST" == true ]]; then
    log "Filtering host reads using BWA..."
    
    # Check if BWA index exists, create if not
    if [[ ! -f "$HOST_GENOME.bwt" ]]; then
        log "Building BWA index for host genome..."
        bwa index "$HOST_GENOME"
    fi
    
    # Align reads to host genome
    log "Aligning reads to host genome..."
    SAM_FILE="$OUTPUT_DIR/${SAMPLE_NAME}_host_aligned.sam"
    
    bwa mem -t "$THREADS" "$HOST_GENOME" "$R1_FILE" "$R2_FILE" > "$SAM_FILE"
    
    # Extract unaligned reads (non-host)
    log "Extracting non-host reads..."
    NONHOST_R1="$OUTPUT_DIR/${SAMPLE_NAME}_nonhost_R1.fastq"
    NONHOST_R2="$OUTPUT_DIR/${SAMPLE_NAME}_nonhost_R2.fastq"
    
    # Get unaligned reads (flag 4 = unaligned)
    samtools view -b -f 4 "$SAM_FILE" | \
    samtools fastq -1 "$NONHOST_R1" -2 "$NONHOST_R2" -0 /dev/null -s /dev/null -n -
    
    # Count non-host reads
    NONHOST_READS=$(wc -l < "$NONHOST_R1" | awk '{print $1/4}')
    HOST_PERCENTAGE=$(awk "BEGIN {printf \"%.1f\", ($TOTAL_READS - $NONHOST_READS) / $TOTAL_READS * 100}")
    
    log "Non-host reads: $NONHOST_READS"
    log "Host reads filtered: ${HOST_PERCENTAGE}%"
    
    # Clean up SAM file
    rm "$SAM_FILE"
    
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
AVG_LENGTH=$(awk 'NR%4==2{sum+=length($0); count++} END{print int(sum/count)}' "$FINAL_R1")

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
