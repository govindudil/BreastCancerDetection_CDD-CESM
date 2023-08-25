# import necessary libraries
import pandas as pd
import numpy as np

# read the csv file
fileName = 'cc_cm.csv'
folderName = '../dataset/cc_cm/'
readFileName = folderName + fileName
df = pd.read_csv(readFileName)

# get the number of samples in each class
class_0 = df[df['isCancerous'] == 0]
class_1 = df[df['isCancerous'] == 1]
print(class_0.shape[0])
print(class_1.shape[0])

# get the count of samples in each class
class_0_count = class_0.shape[0]
class_1_count = class_1.shape[0]

# get the difference between the two classes
diff = class_0_count - class_1_count

# if class 0 has more samples than class 1
if diff > 0:
    # Duplicate class 1 to balance the data
    class_1_dub = class_1.sample(n=(diff%class_1_count), replace=False)
    for i in range(diff//class_1_count):
        class_1_dub = pd.concat([class_1, class_1_dub])

    # Concatenate the two classes
    class_1_dub = pd.concat([class_1, class_1_dub])
    df = pd.concat([class_0, class_1_dub])
# if class 1 has more samples than class 0
elif diff < 0:
    diff = -diff
    # Duplicate class 0 to balance the data
    class_0_dub = class_0.sample(n=(diff%class_0_count), replace=False)
    for i in range(diff//class_0_count):
        class_0_dub = pd.concat([class_0, class_0_dub])
    # Concatenate the two classes
    df = pd.concat([class_1, class_0_dub])

# shuffle the data
df = df.sample(frac=1, replace=False).reset_index(drop=True)

# save the data
saveFileName = folderName + 'balanced_' + fileName
df.to_csv(saveFileName, index=False)

# print the number of samples in each class
print(df['isCancerous'].value_counts())
print(df.shape[0])
