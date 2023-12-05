from modules import data

print("1. Data Retrieval:")

df = data.download_expression_data('https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-PAAD.gistic.tsv.gz')

print(df.head())

print("2. Data Pre-processing:")


print("3. Timestamping Methodology:")


print("4. Identification of rhythmic genes:")


print("5. Validation and Laboratory Experiments:")

print("Finding drug targets")