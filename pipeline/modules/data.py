import requests
import gzip
import pandas as pd
from io import BytesIO

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
    
    