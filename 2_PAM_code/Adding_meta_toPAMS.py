import os
import sys
import pandas as pd

introns_meta = sys.argv[1]
pams_cdt = sys.argv[2]
output_file = sys.argv[3]
s = sys.argv[4] #Size of fragments, in this case 70 bp fragments or smaller!

intron_df = pd.read_csv(introns_meta, sep='\t')


max_order = intron_df.groupby(['gene_id', 'transcript_id'])['order_of_appearance'].transform('max')
intron_df['order_of_appearance'] = intron_df['order_of_appearance'].where(
            intron_df['order_of_appearance'] != max_order, 'Last'
            )

intron_df['order_of_appearance'] = intron_df['order_of_appearance'].apply(lambda x: 'First' if x == 1 else x)

Start_ID_med = intron_df.loc[
            intron_df['size'] >s,
                ['ID_Starting_median', 'order_of_appearance', 'strand',
                         'transcript_id', 'splicing_event', 'gene_id', 'gene_symbol',
                              'ID', 'Chr', 'start', 'end'] #'absolute_start', 'absolute_end'
                ].copy()

End_ID_med= intron_df.loc[
            intron_df['size']>s,
                ['ID_Ending_median','order_of_appearance', 'strand',
                         'transcript_id','splicing_event','gene_id', 'gene_symbol',
                              'ID','Chr', 'start', 'end'] #'absolute_start', 'absolute_end'
                ].copy()
ID_med = intron_df.loc[
            intron_df['size'] <=s,
                ['ID','order_of_appearance', 'strand',
                         'transcript_id', 'splicing_event', 'gene_id', 'gene_symbol',
                               'Chr', 'start', 'end'] #'absolute_start', 'absolute_end'
                ].copy()

def process_row(row):
    type_value=None
    id_parts = row['ID'].replace(".0","").split('_')
    id_start_med = str(row['ID_Starting_median']).replace(".0","")
    id_start_parts = id_start_med.split('_')

    if row['strand'] == '+':
        if 'nan' not in id_start_parts:
            if id_start_parts[2] == id_parts[1]:
                type_value = 'Donor 5\'ss end'
            elif ((id_start_parts[3] == id_parts[2])) :
                type_value = 'Acceptor 3\'ss end'

    elif row['strand'] == '-':
        if 'nan' not in id_start_parts:
            if ((id_start_parts[2]) == (id_parts[1])):
                type_value = 'Acceptor 3\'ss end'
            elif ((id_start_parts[3]) == (id_parts[2])):
                type_value = 'Donor 5\'ss end'

    return pd.Series([type_value])

Start_ID_med['Type'] = Start_ID_med.apply(process_row, axis=1)
Start_ID_med = Start_ID_med.rename(columns={'ID':'EntireIntron','ID_Starting_median':'ID'})
Start_ID_med = Start_ID_med.set_index('ID')

def process_row_end(row):
    id_parts = row['ID'].replace(".0","").split('_')
    id_end_med = str(row['ID_Ending_median']).replace(".0","")
    id_start_parts = id_end_med.split('_')
    type=None
    if row['strand'] == '+':
        if 'nan' not in id_start_parts:
            if id_start_parts[3] == id_parts[2]:
                type = 'Acceptor 3\'ss end'
            elif id_start_parts[2] == id_parts[1]:
                type = 'Donor 5\'ss end'

    elif row['strand'] == '-':
        if 'nan' not in id_start_parts:
            if id_start_parts[3] == id_parts[2]:
                type = 'Donor 5\'ss end'
            elif id_start_parts[2] == id_parts[1]:
                type = 'Acceptor 3\'ss end'

    return pd.Series([type])



End_ID_med['Type'] = End_ID_med.apply(process_row_end, axis=1)
End_ID_med = End_ID_med.rename(columns={'ID':'EntireIntron', 'ID_Ending_median':'ID'})
End_ID_med.set_index('ID', inplace=True)


ID_med['Type'] = ID_med.apply(lambda row: "Donor 5'ss and Acceptor 3'ss" , axis=1)
ID_med.set_index('ID', inplace=True)
ID_med.index = ID_med.index.str.replace(".0", "", regex=False)
Start_ID_med.index = Start_ID_med.index.str.replace(".0", "", regex=False)
End_ID_med.index = End_ID_med.index.str.replace(".0", "", regex=False)

concat_df = pd.concat([Start_ID_med, End_ID_med, ID_med], ignore_index=False)
concat_df.reset_index(inplace=True)

concat_df['transcript_id'] = concat_df['transcript_id'].str.replace('transcript_id=','')
concat_df['trid_coords'] = (concat_df['transcript_id'] + concat_df['ID'].astype(str).str.split(r'[+-]').str[1])
concat_df['gene_id']=concat_df['gene_id'].str.replace('gene_id=','')
concat_df['gene_symbol']=concat_df['gene_symbol'].str.replace('gene_symbol=','')
concat_df['meta'] = (
                concat_df['order_of_appearance'].astype(str) + ' | ' +
                        concat_df['Type'] + ' | ' +
                                concat_df['splicing_event'] + ' | ' +
                                        concat_df['transcript_id'] + ' | ' +
                                                concat_df['start'].astype(str) + ':' +
                                                        concat_df['end'].astype(str) + ' | ' +
                                                                concat_df['gene_id'].astype(str) + '|' + concat_df['gene_symbol']
                                                                ) #'absolute_start', 'absolute_end'

concat_df

concat_df_sizesmall = concat_df[concat_df['trid_coords'].isna()]
concat_df_sizesmall['trid_coords'] = (
                concat_df_sizesmall['transcript_id'] + '_' +
                        concat_df_sizesmall['ID'].astype(str).str.split('_').str[1:].str.join('_')
                        )
concat_df_sizesmall['ID2'] = (concat_df_sizesmall['Chr'] + '_' + concat_df_sizesmall['strand'] + '_' +
                        concat_df_sizesmall['ID'].astype(str).str.split('_').str[1:].str.join('_')
                        )

concat_df_sizesmall.drop('ID', axis=1, inplace=True)
# replacing 'ID2' with 'ID'
concat_df_sizesmall.rename(columns={'ID2': 'ID'}, inplace=True)

merged_df =pd.concat([concat_df, concat_df_sizesmall], ignore_index=True)


pams = pd.read_csv(pams_cdt, sep='\t')
pams['ID'] = pams['FILE'].astype(str).str.split('.').str[0]
pams['gene_symbol'] = pams['ID'].apply(lambda x: concat_df.loc[concat_df['ID'] == x, 'gene_symbol'].values[0] if any(concat_df['ID'] == x) else None
            ).fillna(pams['ID'].apply(lambda x: concat_df_sizesmall.loc[concat_df_sizesmall['ID'] == x, 'gene_symbol'].values[0] if any(concat_df_sizesmall['ID'] == x) else None))

pams['meta'] = pams['ID'].apply(
            lambda x: concat_df.loc[concat_df['ID'] == x, 'meta'].values[0] if any(concat_df['ID'] == x) else None
            ).fillna(pams['ID'].apply(lambda x: concat_df_sizesmall.loc[concat_df_sizesmall['ID'] == x, 'meta'].values[0] if any(concat_df_sizesmall['ID'] == x) else None))


pams.to_csv(output_file, sep='\t')
~                                                                                     
