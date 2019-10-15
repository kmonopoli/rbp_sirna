#!/anaconda2/bin/python
import pandas as pd
from pandas import ExcelWriter
from pandas import ExcelFile


# get list of siRNA human gene targets
# from:
#   BDNA
#   qPCR
#   luciferase
df = pd.read_excel('sirna_data/bdna.xlsx', sheetname='bdna')
print("Column headings:")
print(df.columns)


# get RBP binding data for each gene


# get sequences that correspond to chromosome coordinates

# compile into single file
