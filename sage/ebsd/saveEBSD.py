import pandas as pd
import pickle

def saveEBSDdf(dataframe, file_path):
    """Save a DataFrame to a file."""
    dataframe.to_csv(file_path, index=False)

def loadEBSDdf(file_path):
    """Load a DataFrame from a file."""
    return pd.read_csv(file_path)

def savePhaseName(lst, file_path):
    """Save a list to a file using pickle."""
    with open(file_path, 'wb') as f:
        pickle.dump(lst, f)

def loadPhaseName(file_path):
    """Load a list from a file using pickle."""
    with open(file_path, 'rb') as f:
        return pickle.load(f)
