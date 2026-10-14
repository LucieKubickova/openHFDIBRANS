# -- Python script containing additional functions for OF_caseClass
import numpy as np
import re

def isFloat(val):
    """function to determine if val is float"""
    try:
        float(val)
        return True
    except:
        return False
    
def readDataFromLogFile(file_path):
    skiprows = -1
    endLine = -1
    with open(file_path, "r") as file:
        lines = file.readlines()
    for lineI in range(len(lines)):
        if 'Time = ' in lines[lineI]:
            skiprows = lineI
            break
    for lineI in range(len(lines)):
        if 'End' in lines[lineI]:
            endLine = lineI
            break
    # Extract all numbers using regex
    numbers = []
    for lineI in range(skiprows, endLine):
        line = lines[lineI]
        matches = re.findall(r"[-+]?\d*\.\d+|[-+]?\d+", line)  # Matches integers and decimals
        numbers.extend(map(float, matches))  # Convert to float and add to the list

    # print(numbers)
    # Convert the list to a NumPy array
    moistureSim = np.array(numbers).reshape(-1,2)
    return moistureSim