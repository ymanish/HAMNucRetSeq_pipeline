from Bio import SeqIO
from src.config import path


def contains_non_canonical(sequence):
    """
    Check if a sequence contains 'N' or 'n'.
    Returns True if it contains, otherwise False.
    """
    return 'N' in sequence or 'n' in sequence


def remove_non_canonical_seq(input_fasta, output_fasta):
    """
    Remove sequences with 'N' or 'n' from a FASTA file.
    """
    with open(output_fasta, 'w') as out_handle:
        for record in SeqIO.parse(input_fasta, 'fasta'):
            # Check if "N" or "n" is absent in the sequence
            if 'N' not in record.seq and 'n' not in record.seq:
                SeqIO.write(record, out_handle, 'fasta')

def assert_fasta_record_count(file1, file2):
        count1 = sum(1 for _ in SeqIO.parse(file1, 'fasta'))
        count2 = sum(1 for _ in SeqIO.parse(file2, 'fasta'))
        assert count1 == count2, f"Record counts do not match: {count1} != {count2}"

        print(f"Record counts match: {count1}")


if __name__ == "__main__":
    # # Filter bound_top10000.fa
    # filter_canonical_seq(path.DATA_DIR / 'processed/bound_top10000.fa', path.DATA_DIR / 'processed/bound_top10000_filtered.fa')

    # # Filter unbound_regions.fa
    # filter_canonical_seq(path.DATA_DIR / 'processed/unbound_regions.fa', path.DATA_DIR / 'processed/unbound_regions_filtered.fa')

    
    # # Assert that the filtered files have the same number of records
    # assert_fasta_record_count(path.DATA_DIR / 'processed/bound_top10000_filtered.fa', path.DATA_DIR / 'processed/unbound_regions_filtered.fa')

    print(contains_non_canonical('ATgcccgggttttttaaaatgcc'))