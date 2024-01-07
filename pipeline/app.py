from modules import data, timestamp

print("1. Data Retrieval:")

df_path = data.download_and_extract_tsv_gz('https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-PAAD.gistic.tsv.gz')

print(f"Extracted TSV file saved to: {df_path}")

# Example usage
output, error = timestamp.run_runCIRCUST_R(df_path)

print("Output:")
print(output)

if error:
  print("Error:")
  print(error)

print('Daniels cool addition')

print("2. Data Pre-processing:")


print("3. Timestamping Methodology:")


print("4. Identification of rhythmic genes:")


print("5. Validation and Laboratory Experiments:")

print("Finding drug targets")