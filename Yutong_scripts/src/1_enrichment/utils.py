import os
import sys
import logging
import subprocess
import pybedtools
from typing import Optional, Union, List

CMD=Union[str, List[str]]

def dir_exists(output_path: str):
    """
    Create the output directory if it doesn't exist.
    """
    if not os.path.exists(output_path):
        os.makedirs(output_path, exist_ok=True)

def ensure_parent_dir(file_path: str):
    """
    Ensure that the parent directory of the given file path exists.
    """
    parent_dir = os.path.dirname(file_path)
    if parent_dir and not os.path.exists(parent_dir):
        os.makedirs(parent_dir, exist_ok=True)
    return parent_dir

def run_command(command: CMD, capture: bool = False) -> str:
    """
    Run a shell command and optionally capture its output.
    """
    if isinstance(command, str):
        argv = ["/bin/bash", "-o", "pipefail", "-c", command]
    else:
        argv = command
    if capture:
        cp = subprocess.run(argv, check=True, capture_output=True, text=True)
        return cp.stdout
    else:
        subprocess.run(argv, check=True, text=True)
        return ""
