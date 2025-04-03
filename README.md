# Gene2Go Expander

Expands a given gene2go file by using a Gene Ontology .obo file to find all child terms for each GO ID, then appending child entries to the original Gene2Go file, with the GO ID being replaced with the parent GO ID.

### Usage
```python
python Gene2Go-Expander.py
        [-h]
        --obo-file OBO_FILE
        --gene2go-file GENE2GO_FILE
        [--max-child-num MAX_CHILD_NUM]
        [--include-evidence-codes [INCLUDE_EVIDENCE_CODES ...]]
        [--exclude-evidence-codes [EXCLUDE_EVIDENCE_CODES ...]]
        [--include-tax-ids [INCLUDE_TAX_IDS ...]]
        [--exclude-tax-ids [EXCLUDE_TAX_IDS ...]]
```

### Required Arguments:
- `--obo-file` Path to the Gene Ontology .obo file
- `--gene2go-file` Path to the Gene2Go file (.CSV format) 

### Optional Arguments:
- `-h | --help` Bring up the help documentation
- `--max-child-num` Maximum number of related child GO IDs; any GO ID that exceeds this value will be discarded when saving to `related_go_ids_all.csv`. Default is 2,000.
Default for the following are all an empty list, which doesn't perform any filtering. Use either the include or exclude arguments for evidence codes and/or taxonomy IDs.
- `--include-evidence-codes` Keep any rows with matching evidence codes.
- `--exclude-evidence-codes` Discard any rows with matching evidence codes.
- `--include-tax-ids` Keep any rows with matching taxonomy IDs.
- `--exclude-tax-ids` Discard any rows with matching taxnomy IDs.

Example Usage:
```python
python Gene2Go-Expander.py --obo-file go-basic.obo --gene2go-file gene2go_filterby_taxID.csv --max-child-num 4000 --exclude-evidence-codes IEA ND IKR --include-tax-ids 9606 10090 10116
```
