
i = p1 !!

# Define dataframe
df = pd.read_csv(...)
# TODO Group by the same UniProt accession
# df1 = df.loc[:, ['id1', 'p1', 'd1']]
# df2 = df.loc[:, ['id2', 'p2', 'd2']]
# df3 = pd.concat([d1, df2])
# df3.columns = ['accession', 'position', 'disordered']
# df3 = df3.drop_duplicates()
# df3 = df3.groupby(by=['id1', 'p1', 'd1']).agg().reset_index()
# Define number of rows
n, _ = df.shape
# Initialize dictionary of regions
regions_dict = dict()
# Initialize current key
curr_key = None
# Loop through each row in dataframe
for i, row in df.iterrows():
    # Case we are on the same region
    if (curr_key is not None) and (row.id1, row.id2) == (curr_key[0], curr_key[1]): #change to accession, position...
        # Do only if current residue is disordered
        if row.d1:
            # Update end  of current region
            regions_dict[curr_key] = i
        # Otherwise, reset key
        else:
            curr_key = None
        # # Skip
        # continue
    # Case key is not initialized and not disordered
    elif (curr_key is None) or ((row.id1, row.id2) != (curr_key[0], curr_key[1])):
        # Unset key
        curr_key = None
        # Do only if disordered
        if row.d1:
            # Start new region
            curr_key = row.id1, row.id2, i
            # Store key, set end equal to start
            regions_dict[curr_key] = i
            # # Skip, go to next iteration
            # continue
