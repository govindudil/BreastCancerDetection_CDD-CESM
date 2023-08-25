# read a csv and balance the data according to the number of samples in each class
# output a csv with the same number of samples in each class

import pandas as pd
import numpy as np

# read the csv
fileName = 'dataset.csv'
folderName = '../'
readFileName = folderName + fileName
df = pd.read_csv(readFileName)

# shuffle the data
df = df.sample(frac=1, replace=False).reset_index(drop=True)

# Print the values used for normalization
print("Min: ")
for i in df.iloc[:, 1:-4].min():
    print(i, end=', ')
print()

print("Max: ")
for i in df.iloc[:, 1:-4].max():
    print(i, end=', ')
print()

# Normalize the data except the first and last four columns from -1 to +1
# df.iloc[:, 1:-4] = df.iloc[:, 1:-4].apply(lambda x: (x - x.min()) / (x.max() - x.min()) * 2 - 1)

# save the data
saveFileName = folderName + 'preprocessed_' + fileName
df.to_csv(saveFileName, index=False)

# split into train test
training_samples = int(df.shape[0] * 0.8)//5*5
train_df = df.iloc[:training_samples, :]
test_df = df.iloc[training_samples:, :]

# split the train data into 5 folds
train_df1 = train_df.iloc[:int(train_df.shape[0] * 0.2), :]
train_df2 = train_df.iloc[int(train_df.shape[0] * 0.2):int(train_df.shape[0] * 0.4), :]
train_df3 = train_df.iloc[int(train_df.shape[0] * 0.4):int(train_df.shape[0] * 0.6), :]
train_df4 = train_df.iloc[int(train_df.shape[0] * 0.6):int(train_df.shape[0] * 0.8), :]
train_df5 = train_df.iloc[int(train_df.shape[0] * 0.8):, :]

# save the data
trainFileName = folderName + 'train_' + fileName
testFileName = folderName + 'test_' + fileName
train_df.to_csv(trainFileName, index=False)
test_df.to_csv(testFileName, index=False)

# save the 5 folds into five folders concating other four folds
train_df_train = pd.concat([train_df2, train_df3, train_df4, train_df5])
train_df_train.to_csv(folderName + '1/train.csv', index=False)
train_df1.to_csv(folderName + '1/val.csv', index=False)

train_df_train = pd.concat([train_df1, train_df3, train_df4, train_df5])
train_df_train.to_csv(folderName + '2/train.csv', index=False)
train_df2.to_csv(folderName + '2/val.csv', index=False)

train_df_train = pd.concat([train_df1, train_df2, train_df4, train_df5])
train_df_train.to_csv(folderName + '3/train.csv', index=False)
train_df3.to_csv(folderName + '3/val.csv', index=False)

train_df_train = pd.concat([train_df1, train_df2, train_df3, train_df5])
train_df_train.to_csv(folderName + '4/train.csv', index=False)
train_df4.to_csv(folderName + '4/val.csv', index=False)

train_df_train = pd.concat([train_df1, train_df2, train_df3, train_df4])
train_df_train.to_csv(folderName + '5/train.csv', index=False)
train_df5.to_csv(folderName + '5/val.csv', index=False)

# print the number of samples in each class
print(df['isCancerous'].value_counts())
print(df.shape[0])

# print the number of samples in each class
print("Training Samples:")
print(train_df_train.shape[0])

print("Validation Samples:")
print(train_df1.shape[0])

print("Testing Samples:")
print(test_df.shape[0])


