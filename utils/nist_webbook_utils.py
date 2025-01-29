import requests
import pandas as pd
from urllib.parse import urlencode

def generate_url(
    substance_id, 
    temperature, 
    p_low, 
    p_high, 
    p_increment, 
    digits=5, 
    ref_state="DEF", 
    t_unit="K", 
    p_unit="bar", 
    d_unit="mol/m3", 
    h_unit="kJ/mol", 
    w_unit="m/s", 
    vis_unit="uPa*s", 
    st_unit="N/m"
):
    """
    Generates the NIST Chemistry WebBook query URL.
    """
    base_url = "https://webbook.nist.gov/cgi/fluid.cgi"
    params = {
        "Action": "Data",
        "Wide": "on",
        "ID": substance_id,
        "Type": "IsoTherm",
        "Digits": digits,
        "PLow": p_low,
        "PHigh": p_high,
        "PInc": p_increment,
        "T": temperature,
        "RefState": ref_state,
        "TUnit": t_unit,
        "PUnit": p_unit,
        "DUnit": d_unit,
        "HUnit": h_unit,
        "WUnit": w_unit,
        "VisUnit": vis_unit,
        "STUnit": st_unit
    }
    return f"{base_url}?{urlencode(params)}"

def query_nist_webbook(url):
    """
    Sends an HTTP GET request to the NIST Chemistry WebBook.
    """
    response = requests.get(url)
    if response.status_code == 200:
        return response.text
    else:
        response.raise_for_status()

def parse_to_dataframe(response_text, sep="\t"):
    """
    Parses the response text into a pandas DataFrame with defined data types.
    Assumes the last column is a string, and all others are floating-point numbers.
    """
    # Extract tabular data (assumes header and tabular data are present)
    lines = response_text.splitlines()
    header, *data = [line.split(sep) for line in lines if line.strip()]

    # Create the DataFrame without specifying dtypes initially
    df = pd.DataFrame(data, columns=header)
    
    # Attempt to convert all columns except the last one to float, handling errors
    for col in header[:-1]:
        df[col] = pd.to_numeric(df[col], errors='coerce')  # Converts invalid entries to NaN

    return df

def get_nist_data(
    substance_id, 
    temperature, 
    p_low, 
    p_high, 
    p_increment
):
    """
    Combines URL generation, querying, and parsing into a single function.
    """
    url = generate_url(substance_id, temperature, p_low, p_high, p_increment)
    response_text = query_nist_webbook(url)
    dataframe = parse_to_dataframe(response_text)
    return dataframe

# Example Usage:
# df = get_nist_data("C7732185", 300, 0.9, 1.1, 0.01)
# print(df)
