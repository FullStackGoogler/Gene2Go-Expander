import pandas as pd
import csv
import argparse
import obonet
import ast

def process_ontology(obo_file_path, gene2go_file, max_child_num, include_evidence_codes, exclude_evidence_codes, include_tax_ids, exclude_tax_ids):
    # Step 1: Convert .obo file into a .csv file for later parsing
    csv_file_path = 'obo-data.csv'
    graph = obonet.read_obo(obo_file_path)

    # Get all column names
    headers = ['Node ID']
    for node, data in graph.nodes(data=True):
        for key in data.keys():
            if key not in headers:
                headers.append(key)

    with open(csv_file_path, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(headers)

        for node, data in graph.nodes(data=True):
            row = [node]
            for header in headers[1:]:
                row.append(data.get(header, ''))
            writer.writerow(row)

    print(f".obo file has been converted and saved to {csv_file_path}")

    # Step 2: Find all child terms for each GO ID in the provided .OBO file

    try:
        df = pd.read_csv(csv_file_path, quoting=csv.QUOTE_MINIMAL)
    except pd.errors.ParserError as e:
        print("Parsing error:", e)
        exit()

    # Convert 'is_a' and 'relationship' into a list
    df['is_a'] = df['is_a'].apply(lambda x: ast.literal_eval(x) if pd.notna(x) else [])

    # Function to find all related GO IDs
    def find_child_goids(GO_ID, visited=None):
        if visited is None:
            visited = set()

        if GO_ID in visited:
            return visited

        visited.add(GO_ID)

        # Find all child rows where the current go_id is in their 'is_a' list
        child_rows = df[df['is_a'].apply(lambda x: GO_ID in x)]

        # Recursively find all related child GO IDs
        for child_id in child_rows['Node ID']:
            if child_id not in visited:
                find_child_goids(child_id, visited)

        return visited

    relations = {}

    for index, row in df.iterrows():
        ID = row['Node ID']
        related_goids = find_child_goids(ID)
        relations[(ID, row['name'])] = sorted(list(related_goids))

    output_data = []

    for (main_go_id, name), go_ids in relations.items():
        related_go_ids = go_ids
        if len(related_go_ids) < max_child_num:
            related_go_ids_str = ";".join(related_go_ids)
            output_data.append({'go_id': main_go_id, 'related_goids': related_go_ids_str})

    output_df = pd.DataFrame(output_data)
    output_df.to_csv("related_go_ids_all.csv", index=False)

    print("GO ID to child GO ID relationships have been found and saved to 'related_go_ids_all.csv'.")

    # Step 3: Expand the Gene2Go file using the found child relationships from the Gene Ontology .OBO file

    # Initialize an empty DataFrame with the specified columns
    tempDF = pd.DataFrame(columns=['#tax_id', 'GeneID', 'GO_ID', 'GO_term'])

    # Load the CSV files into DataFrames
    gene2goDF = pd.read_csv(gene2go_file)

    # Include specified evidence codes
    if include_evidence_codes:
        gene2goDF = gene2goDF[gene2goDF['Evidence'].isin(include_evidence_codes)]
    # Exclude specified evidence codes
    if exclude_evidence_codes:
        gene2goDF = gene2goDF[~gene2goDF['Evidence'].isin(exclude_evidence_codes)]

    # Include specified tax IDs
    if include_tax_ids:
        gene2goDF = gene2goDF[gene2goDF['#tax_id'].isin(include_tax_ids)]
    # Exclude specified tax IDs
    if exclude_tax_ids:
        gene2goDF = gene2goDF[~gene2goDF['#tax_id'].isin(exclude_tax_ids)]

    relatedDF = pd.read_csv('related_go_ids_all.csv')

    # Filter the relatedDF based on the condition given
    relatedDF = relatedDF[relatedDF['related_goids'].apply(lambda x: len(x.split(';')) < max_child_num)]

    # Update the related_goids column to exclude the current go_id
    def update_related_goids(row):
        go_list = row['related_goids'].split(';')
        go_list = [go_id for go_id in go_list if go_id != row['go_id']]
        return ';'.join(go_list)

    relatedDF['related_goids'] = relatedDF.apply(update_related_goids, axis=1)

    # Further filter to keep non-empty related_goids and take the head(10)
    relatedDF = relatedDF[relatedDF['related_goids'] != '']
    relatedDF = relatedDF.head(10)

    # Iterate over each row in the relatedDF
    for _, row in relatedDF.iterrows():
        currGO_ID = row['go_id']
        currRelatedGO_IDs = row['related_goids'].split(';')

        existingRows = gene2goDF[gene2goDF['GO_ID'] == currGO_ID]

        newRows = []

        if not existingRows.empty:
            for go_id in currRelatedGO_IDs:
                currRows = gene2goDF[gene2goDF['GO_ID'] == go_id]

                if not currRows.empty:
                    for _, currRow in currRows.iterrows():
                        if currRow['GeneID'] not in existingRows['GeneID'].values:
                            newRow = currRow.copy()
                            newRow['GO_ID'] = currGO_ID
                            newRows.append(newRow)

        if newRows:
            newRowsDF = pd.DataFrame(newRows)
            tempDF = pd.concat([tempDF, newRowsDF], ignore_index=True)

    # Append the tempDF to gene2goDF
    gene2goDF = pd.concat([gene2goDF, tempDF], ignore_index=True)

    # Save the expanded DataFrame to a CSV file
    gene2goDF.to_csv('gene2go_expanded.csv', index=False)
    print("Inputted Gene2Go file has been expanded and saved to gene2go_expanded.csv!")

    # Group by the required columns and aggregate as described
    gene2goSummary = gene2goDF.groupby(['#tax_id', 'GeneID']).agg({
        'GO_ID': lambda x: ';'.join(sorted(set(x))),
        'Evidence': lambda x: ';'.join(sorted(set(x))),
        'Qualifier': lambda x: ';'.join(sorted(set(x))),
        'GO_term': lambda x: ';'.join(sorted(set(x))),
        'PubMed': lambda x: ';'.join(sorted(set(x))),
        'Category': lambda x: ';'.join(sorted(set(x)))
    }).reset_index()

    # Save the summary DataFrame to a CSV file
    gene2goSummary.to_csv('gene2go_summary.csv', index=False)
    print("Gene2Go expansion summary has been saved to gene2go_summary.csv!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Expands a given gene2go file by using a Gene Ontology .obo file to find all child terms for each GO ID, then appending child entries to the original Gene2Go file, with the GO ID being replaced with the parent GO ID.")
    parser.add_argument('--obo-file', required=True, help='Path to the Gene Ontology .obo file')
    parser.add_argument('--gene2go-file', required=True, help='Path to the gene2go .csv file')
    parser.add_argument('--max-child-num', type=int, default=2000, help='Maximum number of related child GO IDs. Default is 2,000; anything with more than the specified limit will be discarded when saving to related_go_ids_all.csv')
    parser.add_argument('--include-evidence-codes', nargs='*', default=[], help='List of evidence codes to include')
    parser.add_argument('--exclude-evidence-codes', nargs='*', default=[], help='List of evidence codes to exclude')
    parser.add_argument('--include-tax-ids', nargs='*', default=[], help='List of tax IDs to include')
    parser.add_argument('--exclude-tax-ids', nargs='*', default=[], help='List of tax IDs to exclude')

    args = parser.parse_args()

    # Convert tax IDs from strings to integers, as DataFrame might store them as integers
    include_tax_ids = [int(tid) for tid in args.include_tax_ids] if args.include_tax_ids else []
    exclude_tax_ids = [int(tid) for tid in args.exclude_tax_ids] if args.exclude_tax_ids else []

    process_ontology(
        args.obo_file,
        args.gene2go_file,
        args.max_child_num,
        args.include_evidence_codes,
        args.exclude_evidence_codes,
        include_tax_ids,
        exclude_tax_ids
    )
