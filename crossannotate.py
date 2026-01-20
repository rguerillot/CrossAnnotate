import pandas as pd
import subprocess
import os
import argparse
import io
import sys

def detect_seq_type(df, col_name):
    """
    Analyzes the character frequency of the sequence column to 
    determine if it is DNA or Protein.
    """
    # Sample the first few sequences to check alphabet
    sample = "".join(df[col_name].dropna().head(10).astype(str)).upper()
    if not sample:
        return "prot"
    
    dna_chars = set("ATGCN")
    dna_count = sum(1 for char in sample if char in dna_chars)
    # If more than 85% of characters are ATGCN, treat as DNA
    return "dna" if (dna_count / len(sample)) > 0.85 else "prot"

def run_diamond_bbh(args):
    """
    Identifies 1-to-1 orthologs between two datasets using 
    DIAMOND-accelerated Bidirectional Best Hit (BBH) logic.
    """
    print(f"\n" + "="*75)
    print("   OrthoMap-DIAMOND: Universal Sequence-Based Table Merger")
    print("="*75)
    
    # 1. Loading Tables
    try:
        df1 = pd.read_csv(args.t1, sep=None, engine='python')
        df2 = pd.read_csv(args.t2, sep=None, engine='python')
    except Exception as e:
        print(f"Error reading input files: {e}"); sys.exit(1)

    # FIXED: Use explicit ID columns if provided, otherwise default to first column
    id_a = args.i1 if args.i1 else df1.columns[0]
    id_b = args.i2 if args.i2 else df2.columns[0]
    
    print(f"Mapping via IDs: Table 1 [{id_a}] <-> Table 2 [{id_b}]")
    
    # 2. Sequence Type Autodetection
    type_a = detect_seq_type(df1, args.s1)
    type_b = detect_seq_type(df2, args.s2)
    
    print(f"Step 1/5: Autodetected formats -> T1: {type_a.upper()}, T2: {type_b.upper()}")
    
    # Create temporary FASTA files for DIAMOND
    with open("ref_a.fasta", "w") as f:
        for _, r in df1.dropna(subset=[args.s1]).iterrows():
            f.write(f">{r[id_a]}\n{r[args.s1]}\n")
    with open("query_b.fasta", "w") as f:
        for _, r in df2.dropna(subset=[args.s2]).iterrows():
            f.write(f">{r[id_b]}\n{r[args.s2]}\n")

    # 3. DIAMOND Execution
    print("Step 2/5: Indexing and searching (Bidirectional)...")
    
    subprocess.run(["diamond", "makedb", "--in", "ref_a.fasta", "-d", "db_a", "--quiet"], check=True)
    subprocess.run(["diamond", "makedb", "--in", "query_b.fasta", "-d", "db_b", "--quiet"], check=True)

    fields = "qseqid sseqid pident length qlen slen evalue bitscore qcovhsp"
    common = ["--threads", str(args.threads), "--evalue", str(args.evalue), 
              "--id", str(args.id), "--query-cover", str(args.cov), 
              "--max-target-seqs", "1", "--outfmt", "6", *fields.split(), "--quiet"]
    
    # Forward Search (A -> B)
    mode_fwd = "blastp" if type_a == "prot" else "blastx"
    print(f"Step 3/5: Forward search ({mode_fwd})...")
    res_ab = subprocess.run(["diamond", mode_fwd, "-q", "ref_a.fasta", "-d", "db_b"] + common, capture_output=True, text=True).stdout
    
    # Reverse Search (B -> A)
    mode_rev = "blastp" if type_b == "prot" else "blastx"
    print(f"Step 4/5: Reverse search ({mode_rev})...")
    res_ba = subprocess.run(["diamond", mode_rev, "-q", "query_b.fasta", "-d", "db_a"] + common, capture_output=True, text=True).stdout

    # 4. Orthology Resolution
    cols = fields.split()
    df_ab = pd.read_csv(io.StringIO(res_ab), sep='\t', names=cols).sort_values('bitscore', ascending=False).drop_duplicates('qseqid')
    df_ba = pd.read_csv(io.StringIO(res_ba), sep='\t', names=cols).sort_values('bitscore', ascending=False).drop_duplicates('qseqid')

    # Reciprocal Merge (Crucial for 1-to-1 Orthology)
    ortho = pd.merge(df_ab, df_ba, left_on=['qseqid', 'sseqid'], right_on=['sseqid', 'qseqid'])
    
    if ortho.empty:
        print("Warning: No orthologs found. Check if ID columns and sequence types are correct.")
        sys.exit(0)

    # 5. Quality Metrics
    final = ortho[['qseqid_x', 'sseqid_x', 'pident_x', 'evalue_x', 'bitscore_x', 'qlen_x', 'slen_x', 'qcovhsp_x']].copy()
    final.columns = [id_a, id_b, 'identity_pct', 'evalue', 'bitscore', 'len_a', 'len_b', 'cov_a_pct']
    final['len_ratio'] = final['len_a'] / final['len_b']

    # 6. Intelligent Passthrough
    if args.passthrough:
        if args.passthrough == '*':
            anno_cols = df1.columns.tolist()
        else:
            anno_cols = [id_a] + [c.strip() for c in args.passthrough.split(',') if c.strip() in df1.columns]
        final = pd.merge(final, df1[list(set(anno_cols))], on=id_a, how='left')

    final.to_csv(args.out, index=False)

    # 7. Final Report
    print("\n" + "*"*75)
    print("                      EXECUTION SUMMARY")
    print("*"*75)
    print(f"Total Orthologs mapped: {len(final)}")
    print(f"Mean Identity:         {final['identity_pct'].mean():.2f}%")
    print(f"Table 1 Format:        {type_a.upper()}")
    print(f"Table 2 Format:        {type_b.upper()}")
    print("*"*75)
    print(f"Output: {args.out}\n")

    # Cleanup
    for f in ["ref_a.fasta", "query_b.fasta", "db_a.dmnd", "db_b.dmnd"]:
        if os.path.exists(f): os.remove(f)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="OrthoMap-DIAMOND: Sequence-based table merger with autodetection.",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    # Tables
    parser.add_argument('--t1', required=True, help='Reference Table (UniProt)')
    parser.add_argument('--t2', required=True, help='Query Table (RAST)')
    
    # ID Columns (The Fix)
    parser.add_argument('--i1', help='ID column name for T1 (Default: 1st col)')
    parser.add_argument('--i2', help='ID column name for T2 (Default: 1st col)')
    
    # Sequence Columns
    parser.add_argument('--s1', required=True, help='Sequence column in T1')
    parser.add_argument('--s2', required=True, help='Sequence column in T2')
    
    # Thresholds
    parser.add_argument('--id', default=90.0, type=float, help='Min Identity %')
    parser.add_argument('--cov', default=90.0, type=float, help='Min Coverage %')
    parser.add_argument('--evalue', default=1e-10, type=float)
    
    # Performance & Output
    parser.add_argument('--threads', default=4, type=int)
    parser.add_argument('--passthrough', help='Columns to transfer from T1 (use "*" for all)')
    parser.add_argument('--out', default='ortho_results.csv')

    run_diamond_bbh(parser.parse_args())
