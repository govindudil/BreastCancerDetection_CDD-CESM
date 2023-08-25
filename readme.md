# A breast cancer detection application with Grammatical Algorithms in C++ for Evolution (GRACE) and Gray Level Co-Occurrence Matrix (GLCM).

---
## Structure of the repository
---

1. **Feature_Extraction**: This directory contains methods and data related to feature extraction using GLCM.
   - **Dataset**: Contains datasets used for feature extraction.
     - **Input**: Contains input datasets.
       - **CDD-CESM**: Dataset used in the study.
   - **Images**: Contains image names used in the project.
     - **test**: Contains test image names.
       - **cc_cm**: Images related to this category.
       - **cc_dm**: Images related to this category.
       - **mlo_cm**: Images related to this category.
       - **mlo_dm**: Images related to this category.
     - **train**: Contains training image names.
       - **cc_cm**: Images related to this category.
       - **cc_dm**: Images related to this category.
       - **mlo_cm**: Images related to this category.
       - **mlo_dm**: Images related to this category.
   - **markedImages**: Contains images that have been annotated.

2. **Grammatical_Evolution**: This directory is related to the grammatical evolution aspect of the project.
   - **best_individuals**: Contains data about the best individuals in the grammatical evolution process.
   - **grammar**: Contains grammar files.
   - **src**: Source files related to grammatical evolution.
   - **template**: Contains template files.
   - **tmp**: Temporary files.

3. **Results**: Contains results from various experiments.
   - **cc_cm**: Results related to this category.
   - **cc_dm**: Results related to this category.
   - **mlo_cm**: Results related to this category.
   - **mlo_dm**: Results related to this category.

4. **dataset**: Contains dataset after extracting features via GLCM used in the project.
   - **cc_cm**: Dataset and scripts related to this category.
     - **Scripts**: Scripts related to this dataset.
   - **cc_dm**: Dataset and scripts related to this category.
     - **Scripts**: Scripts related to this dataset.
   - **mlo_cm**: Dataset and scripts related to this category.
     - **Scripts**: Scripts related to this dataset.
   - **mlo_dm**: Dataset and scripts related to this category.
     - **Scripts**: Scripts related to this dataset.

This gives a high-level overview of the repository's structure. Each directory and sub-directory contains specific files and data related to its purpose in the project.



---
## Building and running the project
---

### 1. Feature Extraction:
1. Navigate to the `Feature_Extraction` directory.
2. Before building, you need to edit the `main.cpp` file:
   - Open `main.cpp` in a text editor.
   - Locate lines 26 and 27.
   - Modify the filenames depending on which file you want to generate (options: `cc_cm`, `cc_dm`, `mlo_cm`, `mlo_dm`).
   - Save and close the file.
3. Build the program using CMake:
   ```bash
   cmake .
   ```
4. After the build configuration is generated, compile the program using:
   ```bash
   make
   ```
5. Once the compilation is successful, run the program using:
   ```bash
   ./Dissertation
   ```

### 2. Dataset Preparation:
1. After running the feature extraction, a file will be generated in the `Feature_Extraction/Dataset/Output` directory.
2. Copy this generated file to the appropriate folder inside the `dataset` directory based on the file generated. For example, if you generated `cc_cm`, then navigate to the `dataset/cc_cm` directory.
3. Rename the copied file to `dataset.csv`:
   ```bash
   mv [generated_filename] dataset.csv
   ```
4. Run the `preprocess.sh` script to preprocess the dataset:
   ```bash
   ./preprocess.sh
   ```

### 3. Grammatical Evolution:
1. Navigate to the `Grammatical_Evolution` directory.
2. Build the GE repository using:
   ```bash
   make
   ```
3. After successful compilation, run the GEGCC binary to evolve the population:
   ```bash
   ./GEGCC
   ```

---
Make sure you have all the necessary dependencies installed and the required permissions to execute scripts and binaries. 
If you have any questions or need further assistance, please feel free to contact me at:
- Email: [govindudil@icloud.com](mailto:govindudil@icloud.com)
- University Email: [22206094@studentmail.ul.ie](mailto:22206094@studentmail.ul.ie)

I'll be happy to help!

---