import os

def dir_exists(output_path: str):
    """
    Create the output directory if it doesn't exist.
    """
    if not os.path.exists(output_path):
        os.makedirs(output_path, exist_ok=True)
