import pandas as pd
import pickle

def save_ebsd_df(dataframe, file_path):
    """Save a DataFrame to a file."""
    dataframe.to_csv(file_path, index=False)

def load_ebsd_df(file_path):
    """Load a DataFrame from a file."""
    return pd.read_csv(file_path)

def save_phase_name(lst, file_path):
    """Save a list to a file using pickle."""
    with open(file_path, 'wb') as f:
        pickle.dump(lst, f)

def load_phase_name(file_path):
    """Load a list from a file using pickle."""
    with open(file_path, 'rb') as f:
        return pickle.load(f)
