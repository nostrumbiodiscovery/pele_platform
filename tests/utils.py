import math
import os
import re


def truncate(f, n):
    return math.floor(f * 10 ** n) / 10 ** n


def check_file(
    folder, filename, values, errors, delimiter=None, truncate_digits_to=4
):
    filename = os.path.join(folder, filename)

    with open(filename, "r") as f:
        lines = f.readlines()

        if delimiter is None:
            for value in values:
                if value not in "".join(lines):
                    errors.append(value)

        elif isinstance(delimiter, str):
            for value in values:
                split_values = value.split(delimiter)
            for split_value in split_values:
                if split_value.replace(".", "").replace("-", "").isnumeric():
                    if str(
                        truncate(float(split_value), truncate_digits_to)
                    ) not in "".join(lines):
                        errors.append(value)
                else:
                    if split_value not in "".join(lines):
                        errors.append(split_value)
        else:
            raise TypeError("Wrong delimiter type")

    return errors


def check_file_regex(directory, file, patterns, errors):
    """
    Checks file for strings matching regex patterns.

    Parameters
    -----------
    directory : str
        Path to folder.
    file : str
        Name of the file to check.
    patterns : List[str]
        List of regex patterns to check.
    errors : List[str]
        List of existing errors (or empty list).
    """
    path = os.path.join(directory, file)

    with open(path) as f:
        content = f.read()

    for pattern in patterns:
        if not re.search(pattern, content):
            errors.append(pattern)

    return errors
