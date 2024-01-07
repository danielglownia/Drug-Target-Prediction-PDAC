import requests
import gzip
import pandas as pd
from io import BytesIO
import os
import tempfile
import shutil
from gzip import GzipFile

def download_expression_data(url):
    # Step 1: Download the file
    response = requests.get(url)
    
    if response.status_code == 200:
        # Step 2: Unzip the file
        with gzip.GzipFile(fileobj=BytesIO(response.content), mode='rb') as f:
            # Step 3: Read data into a pandas DataFrame
            df = pd.read_csv(f, sep='\t')
            return df
    else:
        print(f"Failed to download the file. Status code: {response.status_code}")
        return None

def download_and_extract_tsv_gz(url):

  # Check if the URL is valid for a tsv.gz file
  if not url.endswith(".tsv.gz"):
    raise ValueError("URL must be for a tsv.gz file")

  # Download the file
  try:
    response = requests.get(url, stream=True)
    response.raise_for_status()
  except requests.exceptions.RequestException as e:
    raise e

  # Create a temporary directory to store the downloaded file
  with tempfile.TemporaryDirectory() as tmpdir:
    # Save the downloaded file to a temporary file
    tmp_file_path = os.path.join(tmpdir, os.path.basename(url))
    with open(tmp_file_path, "wb") as f:
      for chunk in response.iter_content(chunk_size=1024):
        f.write(chunk)

    # Extract the downloaded file
    extracted_file_path = tmp_file_path[:-3]  # remove ".gz" extension
    with GzipFile(tmp_file_path, "rb") as gzip_file:
      with open(extracted_file_path, "wb") as out_file:
        shutil.copyfileobj(gzip_file, out_file)

  # Return the path to the extracted tsv file
  return extracted_file_path

    