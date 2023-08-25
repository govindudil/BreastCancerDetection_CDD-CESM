# import necessary libraries
import pandas as pd
import numpy as np

# read the csv file
fileName = 'cc_cm.csv'
folderName = '../dataset/cc_cm/'
readFileName = folderName + fileName
df = pd.read_csv(readFileName)

# shuffle the data
df = df.sample(frac=1, replace=False).reset_index(drop=True)

# print the minimum and maximum values for normalization
print("Min: ")
for i in df.iloc[:, 1:-4].min():
    print(i, end=', ')
print()

print("Max: ")
for i in df.iloc[:, 1:-4].max():
    print(i, end=', ')
print()

# normalize the data except the first and last four columns from -1 to +1
df.iloc[:, 1:-4] = df.iloc[:, 1:-4].apply(lambda x: (x - x.min()) / (x.max() - x.min()) * 2 - 1)

# save the normalized data
saveFileName = folderName + 'normalized_' + fileName
df.to_csv(saveFileName, index=False)

# print the number of samples in each class
print(df['isCancerous'].value_counts())
print(df.shape[0])
