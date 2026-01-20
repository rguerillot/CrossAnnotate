# CrossAnnotate: Sequence-Based Table Merger

**CrossAnnotate** is a tool designed to integrate disparate data tables by identifying 1-to-1 orthologs. It facilitates the transfer of curated functional annotations (e.g., UniProt) to newly sequenced clinical isolates or unannotated genomes using **Bidirectional Best Hit (BBH)** logic.


## Core Features
* **Homology-Driven Joining:** Merges tables based on sequence identity rather than fragile ID-matching.
* **Smart Autodetect:** Automatically samples sequences to determine if they are DNA or Protein, switching between DIAMOND `blastp` and `blastx` modes seamlessly.
* **Precision Metrics:** Calculates `identity_pct`, `len_ratio`, and `coverage` to help identify truncations, pseudogenes, or assembly artifacts.
* **Bulk Annotation Transfer:** Use the `--passthrough` flag to map entire suites of functional data in one pass.


## Installation

### Prerequisites
1. **DIAMOND Aligner**: Must be installed and in your PATH.
   ```bash
   conda install -c bioconda diamond
2. **Python Environment**: Requires Python 3.7+ and the pandas library.
   ```bash
   pip install pandas
   
## Usage
The script requires two tables (CSV, TSV, or TXT) and the names of the columns containing the sequences and unique identifiers.

**Command Example**
```bash
python3 table_sequence_merger.py \
    --t1 reference_uniprot.csv \
    --t2 clinical_isolate_prokka.csv \
    --i1 "Entry" \
    --i2 "feature_id" \
    --s1 "Sequence" \
    --s2 "aa_sequence" \
    --passthrough "*" \
    --threads 8 \
    --out "cross_annotated_results.csv"
```

