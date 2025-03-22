
import pandas as pd


def txt2list_of_lines(s_filename):
    """
    the output is a list where each item is a string of one line
    """
    with open(s_filename, 'r') as file:
        l_lines = file.readlines()
    return l_lines



def ll2df(ll_input):
    # Check that the input is a list.
    if not isinstance(ll_input, list):
        raise ValueError(f"Expected input to be a list of lists, but got {type(ll_input).__name__}.")
    
    # Check that every element in the list is itself a list.
    if any(not isinstance(row, list) for row in ll_input):
        raise ValueError("Expected input to be a list of lists. Every element in the input list must also be a list (representing a row)")
    
    return pd.DataFrame(ll_input)


def dl2df(dl_input):
    
    # Check that the input is a dictionary.
    if not isinstance(dl_input, dict):
        raise ValueError(f"Expected input to be a dictionary of lists, but got {type(dl_input).__name__}.")
    
    # Check that every value in the dictionary is a list.
    if any(not isinstance(val, list) for val in dl_input.values()):
        raise ValueError("Expected input to be a dictionary of lists. Every value in the dictionary must be a list (representing a column)")
    
    return pd.DataFrame(dl_input)



def df2ll(df_input):
    """
    Converts a pandas DataFrame into a list of lists where each nested list is a row from the input DataFrame.
    
    Parameters:
    df (pandas.DataFrame): The input DataFrame.
    
    Returns:
    list of lists: Each sublist represents a row in the DataFrame.
    
    Example:
    input df:
    1  3
    2  4

    output ll:
    [[1, 3], [2, 4]]


    example usage:
    ll_hello = cl.df2ll(df_hi)
    
    """
    # Validate that input is a DataFrame
    if not isinstance(df_input, pd.DataFrame):
        raise ValueError("The input must be a pandas DataFrame.")
    
    return df_input.values.tolist()


def df2dl(df_input):
    """
    Converts a pandas DataFrame into a dict of lists where each nested list is a column from the input DataFrame.
    
    Parameters:
    df (pandas.DataFrame): The input DataFrame.
    
    Returns:
    dic of lists: Each sublist represents a row in the DataFrame.
    
    Example:
    input df:
    A  B
    1  3
    2  4

    output dl:
    {'A': [1, 3], 'B': [2, 4]}


    example usage:
    ll_hello = cl.df2dl(df_hi)
    
    """

    # Verify that the input is a pandas DataFrame.
    if not isinstance(df_input, pd.DataFrame):
        raise ValueError("The input must be a pandas DataFrame.")
    
    return df_input.to_dict('list')