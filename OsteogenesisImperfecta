#Final project, Osteogenesis Imperfecta 
# Authors: Usman Ahmed, Noha Elhaj, Tazmeen Gill

#anaconda2023
#conda activate /home/user/bioenv
#conda install -c bioconda blast # -----> This installs the NCBI's blast tool functionality into our environment
chmod +x final.sh

#Step 1 Genomic Data Retrieval: fetch and save relevant sequences ================================================

#Sequence 1 --> Healthy sequence COL1A1
# FASTA healthy genome 
# seq_start seq_stop = selecting the gene from the chromosome, as described on https://www.ncbi.nlm.nih.gov/gene/1277
# if healthy reference is empty, erorr message + exit the program 
echo "Retrieving NC_000017.11..."
if efetch -db nucleotide -id NC_000017.11 -seq_start 50184101 -seq_stop 50201631 -strand 1 -format fasta > Healthy_COL1A1.fasta; then 
    echo "Successfully retrieved NC_000017.11."
else
    echo "Failed to retrieve NC_000017.11."
    exit 1
fi

#Sequence 2 --> Healthy sequence COL1A2 
# selecting gene from the chromosome, as described on https://www.ncbi.nlm.nih.gov/gene/1278 
# if healthy reference is empty, erorr message + exit the program 
echo "Retrieving NC_000007.14..."
if efetch -db nucleotide -id NC_000007.14 -seq_start 94394895 -seq_stop 94431227 -strand 1 -format fasta > Healthy_COL1A2.fasta; then 
    echo "Successfully retrieved NC_000007.14."
else
    echo "Failed to retrieve NC_000007.14."
    exit 1 
fi

#Sequence 3 --> OI mutated sequence COL1A1
#FASTQ disease genome, just 1 is good; limited size based on storage capacity on cocalc
SRR1a1=SRR16482958 # 13yo male
echo "Retrieving sequence $SRR1a1..."
mkdir ./RAW
prefetch $SRR1a1 -O ./RAW
fastq-dump ./RAW/${SRR1a1}/${SRR1a1}.sra -O ./RAW # outputs SRR_Name.fastq
# if disease reference empty, error message + exit the program exit 1 = minor error exit
if [ -s ./RAW/${SRR1a1}.fastq ]; then
    echo "Successfully retrieved sequence $SRR1a1."
else
    echo "Failed to retrieve sequence $SRR1a1."
    exit 1
fi

#Sequence 4 --> OI mutated sequence COL1A2
SRR1a2=SRR16482959 # 12yo male
echo "Retrieving sequence $SRR1a2..."
prefetch $SRR1a2 -O ./RAW
fastq-dump ./RAW/${SRR1a2}/${SRR1a2}.sra -O ./RAW
# if disease reference empty, error message + exit the program exit 1 = minor error exit
if [ -s ./RAW/${SRR1a2}.fastq ]; then
    echo "Successfully retrieved sequence $SRR1a2."
else
    echo "Failed to retrieve sequence $SRR1a2."
    exit 1 
fi


# FASTQC Quality Check of the fastq file after download ; 
# .zip output; redirected with -o to desired directory
mkdir ./QC_results
QCDIR= ./QC_results
echo Quality checking COL1A1 raw read...
fastqc -o $QCDIR $RAW_DIR/${SRR1a1}.fastq
echo Finished FastQC for COL1A1

echo Quality checking COL1A2 raw read...
fastqc -o $QCDIR $RAW_DIR/${SRR1a2}.fastq
echo Finished FastQC for COL1A2
# we dont use these for anything tho.


# Step 2 Genetic Variation Analysis: Comparing sequences ============================================================
# "Compare the genetic sequences of affected individuals and healthy controls. 
# Identify any common or unique variations in the disease group."

#blast tool for aligning mutated and healthy COL1A1
blastn -query ./RAW/${SRR1a1}.fastq -subject Healthy_COL1A1.fasta -outfmt 6 -out COL1A1_mutvshealthy.txt

# Parse the BLAST output to extract and report details for each alignment
if test -s COL1A1_mutvshealthy.txt; then
    awk 'BEGIN {
        print "BLAST Results for mutated ('$SRR1a1') vs healthy (NC_000017.11) COL1A1:"
        print "-------------------------------"
    }
    {
        query_id = $1;
        subject_id = $2;
        percent_identity = $3;
        alignment_length = $4;
        mismatches = $5;
        gap_opens = $6;
        q_start = $7;
        q_end = $8;
        s_start = $9;
        s_end = $10;
        e_value = $11;
        bit_score = $12;

        print "Query ID: " query_id;
        print "Subject ID: " subject_id;
        print "Percentage of Identity: " percent_identity "%";
        print "Alignment Length: " alignment_length;
        print "Number of Mismatches: " mismatches;
        print "Number of Gap Opens: " gap_opens;
        print "Query Start: " q_start;
        print "Query End: " q_end;
        print "Subject Start: " s_start;
        print "Subject End: " s_end;
        print "E-value: " e_value;
        print "Bit Score: " bit_score;
        print "-------------------------------"
    }' COL1A1_mutvshealthy.txt
else
    echo "No BLAST results found or the output file is empty."
    exit 1
fi
# ^ which awk elements are 'common or unique variations'?

#blast tool for aligning mutated and healthy COL1A2
blastn -query ./RAW/${SRR1a2}.fastq -subject Healthy_COL1A2.fasta -outfmt 6 -out COL1A2_mutvshealthy.txt

if test -s COL1A2_mutvshealthy.txt; then
    awk 'BEGIN {
        print "BLAST Results for mutated ('$SRR1a2') vs healthy (NC_000007.14) COL1A2:"
        print "-------------------------------"
    }
    {
        query_id = $1;
        subject_id = $2;
        percent_identity = $3;
        alignment_length = $4;
        mismatches = $5;
        gap_opens = $6;
        q_start = $7;
        q_end = $8;
        s_start = $9;
        s_end = $10;
        e_value = $11;
        bit_score = $12;

        print "Query ID: " query_id;
        print "Subject ID: " subject_id;
        print "Percentage of Identity: " percent_identity "%";
        print "Alignment Length: " alignment_length;
        print "Number of Mismatches: " mismatches;
        print "Number of Gap Opens: " gap_opens;
        print "Query Start: " q_start;
        print "Query End: " q_end;
        print "Subject Start: " s_start;
        print "Subject End: " s_end;
        print "E-value: " e_value;
        print "Bit Score: " bit_score;
        print "-------------------------------"
    }' COL1A2_mutvshealthy.txt
else
    echo "No BLAST results found or the output file is empty."
    exit 1
fi




# Step 3 Functional Annotation: =======================================================================================
# "For the identified genetic variations, perform basic functional annotation using tools like Variant Effect Predictor (VEP) or SnpEff. 
# Determine the potential functional consequences (e.g., missense, nonsense, frameshift) of the variations."


# Step 4 Variant Frequency and Population Analysis:  =======================================================================
# "Investigate the frequency of the identified genetic variations in different populations using tools like ExAC, gnomAD, or 1000 Genomes Project. 
# Discuss any potential population-specific trends in genetic variation associated with the disease."



