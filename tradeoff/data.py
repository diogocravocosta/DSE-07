import pandas as pd


def read_data(filepath: str) -> pd.DataFrame:
    """
    Creates a Pandas DataFrame from an Excel file.
    """
    file = pd.ExcelFile(filepath)

    concepts_df = file.parse('Concepts')
    criteria_df = file.parse('Criteria')
    matrix_df = file.parse('Trade-off Matrix')

    return concepts_df, criteria_df, matrix_df



