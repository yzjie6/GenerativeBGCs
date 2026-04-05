import os

# High-GC optimal codons (e.g. widely used for Streptomyces/E.coli general synthetic expression)
# Maps 1 amino acid to its most frequent high-GC codon for stable chemical synthesis.
CODON_TABLE = {
    'A': 'GCG', 'R': 'CGC', 'N': 'AAC', 'D': 'GAC', 'C': 'TGC',
    'Q': 'CAG', 'E': 'GAG', 'G': 'GGC', 'H': 'CAC', 'I': 'ATC',
    'L': 'CTG', 'K': 'AAG', 'M': 'ATG', 'F': 'TTC', 'P': 'CCG',
    'S': 'AGC', 'T': 'ACC', 'W': 'TGG', 'Y': 'TAC', 'V': 'GTG',
    '*': 'TGA'  # Stop codon
}

def reverse_translate(aa_sequence):
    """Convert an amino acid sequence to DNA using the predefined codon table."""
    dna = []
    for aa in aa_sequence:
        # Default to neutral 'NNN' if valid AA code not found
        dna.append(CODON_TABLE.get(aa.upper(), 'NNN'))
    return "".join(dna).lower()

def write_gbk(chimera, protein_dict, output_path):
    """
    Format and write a standardized GenBank file.
    chimera: The chimera dict from main.py
    protein_dict: A lookup table Mapping protein_id -> sequence
    """
    p_ids = chimera["protein_ids"]
    
    features = []
    current_bp = 1
    full_dna = ""
    
    linker_dna = ""
    linker_aa = ""
    if chimera["rescued"] and chimera["rescue_linker"]:
        if "GGS" in chimera["rescue_linker"]:
            linker_aa = "GGSGGSGG" if "Extended" in chimera["rescue_linker"] else "GGSGG"
        elif "Alpha" in chimera["rescue_linker"]:
            linker_aa = "EAAAK"
        linker_dna = reverse_translate(linker_aa)
    
    for i, p_id in enumerate(p_ids):
        aa_seq = protein_dict.get(p_id, "")
        if not aa_seq: continue
            
        dna_seq = reverse_translate(aa_seq)
        if i == len(p_ids) - 1:
            # Last protein gets a stop codon
            dna_seq += CODON_TABLE['*'].lower()
            
        end_bp = current_bp + len(dna_seq) - 1
        
        # Determine source (host or donor)
        note = "Host backbone"
        if p_id == chimera["donor_protein"]:
            note = f"Inserted Donor part from {chimera['donor_organism']} ({chimera['donor_bgc']})"
            
        feat = f"""     CDS             {current_bp}..{end_bp}
                     /protein_id="{p_id}"
                     /note="{note}"
                     /translation="{aa_seq}"
"""
        features.append(feat)
        full_dna += dna_seq
        current_bp = end_bp + 1
        
        # Insert linker if not the last protein and rescued
        if i < len(p_ids) - 1 and linker_dna:
            end_bp_linker = current_bp + len(linker_dna) - 1
            feat_link = f"""     misc_feature    {current_bp}..{end_bp_linker}
                     /note="Computational Linker Insertion ({chimera['rescue_linker']})"
                     /translation="{linker_aa}"
"""
            features.append(feat_link)
            full_dna += linker_dna
            current_bp = end_bp_linker + 1
            
    # Now write the native auxiliary genes (regulators, transporters, tailoring)
    for aux_p in chimera.get("aux_proteins", []):
        aa_seq = aux_p["sequence"]
        if not aa_seq: continue
        dna_seq = reverse_translate(aa_seq) + CODON_TABLE['*'].lower()
        end_bp = current_bp + len(dna_seq) - 1
        
        # Check if it was swapped 
        note = "Native auxiliary/tailoring gene"
        
        # Simple detection if it comes from donor
        if "donor_aux" in chimera and aux_p["protein_id"] in chimera["donor_aux"]:
            note = f"Swapped tailoring gene from {chimera['donor_bgc']}"
            
        feat = f"""     CDS             {current_bp}..{end_bp}
                     /protein_id="{aux_p['protein_id']}"
                     /note="{note}"
                     /product="{aux_p['product']}"
                     /translation="{aa_seq}"
"""
        features.append(feat)
        full_dna += dna_seq
        current_bp = end_bp + 1
            
    total_bp = len(full_dna)
    
    # Format origin block
    origin_lines = []
    for i in range(0, total_bp, 60):
        chunk = full_dna[i:i+60]
        # format: line number right-justified padding length 9, then chunks of 10 separated by space
        chunk_chunks = [chunk[j:j+10] for j in range(0, len(chunk), 10)]
        line = f"{(i+1):>9} " + " ".join(chunk_chunks)
        origin_lines.append(line)
        
    origin_block = "\n".join(origin_lines)
    features_block = "".join(features).rstrip()
    
    gbk_content = f"""LOCUS       {chimera['chimera_id']}               {total_bp} bp    DNA     linear   SYN 01-JAN-2026
DEFINITION  Synthetic Chimeric Assembly Line {chimera['chimera_id']}
ACCESSION   {chimera['chimera_id']}
VERSION     {chimera['chimera_id']}.1
KEYWORDS    .
SOURCE      synthetic construct
  ORGANISM  synthetic construct
            other sequences; artificial sequences.
FEATURES             Location/Qualifiers
     source          1..{total_bp}
                     /organism="synthetic construct"
                     /mol_type="other DNA"
                     /note="Generated by GenerativeBGCs v4 Framework"
                     /host="{chimera['host_organism']}"
                     /chimera_djcs="{chimera['mean_djcs']}"
                     /donor_bgc="{chimera['donor_bgc']}"
{features_block}
ORIGIN
{origin_block}
//
"""
    with open(output_path, 'w') as f:
        f.write(gbk_content)


def export_chimeras_to_gbk(chimeras, all_entries, output_dir="results/gbk", top_n=10):
    os.makedirs(output_dir, exist_ok=True)
    
    # Build dictionary for O(1) protein lookup
    protein_dict = {}
    for entry in all_entries:
        for p in entry.get("core_proteins", []):
            protein_dict[p["protein_id"]] = p["sequence"]
        for p in entry.get("aux_proteins", []):
            protein_dict[p["protein_id"]] = p["sequence"]
            
    exported_count = 0
    exported_files = []
    
    # Process the top N chimeras
    for c in chimeras[:top_n]:
        out_file = os.path.join(output_dir, f"{c['chimera_id']}.gbk")
        write_gbk(c, protein_dict, out_file)
        exported_files.append(out_file)
        exported_count += 1
        
    return exported_count, exported_files
