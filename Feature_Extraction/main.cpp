//
// Created by govindu on 5/26/23.
//

#include <iostream>
#include <opencv2/core/cvstd.hpp>
#include <opencv2/opencv.hpp>
#include <dirent.h>
#include <list>
#include "glcmMain.h"
#include <cstdio>
#include <ctime>
#include <time.h>
#include <dirent.h>
#include <fstream> 
#include <vector>
#include <string>
#include <json/json.h>
#include <omp.h>
#include <mutex>

#define WINDOW_SIZE 250
#define LABEL_THRESHOLD 5*255/100 // 5% of the maximum value of the image
#define OVERLAPS 5

#define FILE_NAME "cc_cm"
#define FILE_PATH "Images/train"
#define SAVE_FILE_PATH "Dataset/Output/Train"

using namespace std;
using namespace cv;


// Define a struct to hold metadata information about an image
typedef struct MetaData {
    string patientID;
    string side;
    string energy;
    string view;
} MetaData;

// Define some file paths to be used later in the code
const char * savePath = "outputTmp";
const char * CSVpath = "Dataset/Input/CDD-CESM/Radiology_hand_drawn_segmentations_v2.csv";
const char * pathCM = "Dataset/Input/CDD-CESM/Subtracted images of CDD-CESM";
const char * pathDM = "Dataset/Input/CDD-CESM/Low energy images of CDD-CESM";
const char * pathRef = "markedImages";
const char * image = "P2_L_DM_CC";

// Define two mutexes for synchronization
std::mutex mtxCounter;
std::mutex mtxSaveResults;

// Define a function to get a list of files in a directory
vector<char *> files_(const char *  path){
    DIR *dir; struct dirent *diread;
    vector<char *> files;

    if ((dir = opendir(path)) != nullptr) {
        while ((diread = readdir(dir)) != nullptr) {
            files.push_back(diread->d_name);
        }
        closedir (dir);
    } else {
        perror ("opendir");
    }
    return files;
}

// This function reads an image from a given path and calculates 15 GLCM features on it
int read(const string& path){

    // Read the image from the given path
    Mat img;
    img = imread(path,0 );

    // Start the timer for the overall process time
    clock_t START,END;
    double cpu_time_used_MAIN;
    START = clock();

    // Define an array of GLCM features to calculate
    char *c[] = {"asm", "contrast", "corr", "var", "idm", "savg", "svar", "sentropy", "entropy", "dvar", "dentropy", "icorr1", "icorr2", "maxcorr", "sum"};

    // Start the timer for the GLCM calculation time
    clock_t start,end;
    double cpu_time_used;
    start = clock();

    // Calculate the GLCM features on the image
    double *result = calculator(img,c,15,0,1);

    // Release the image from memory
    img.release();

    // End the timer for the GLCM calculation time
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    // End the timer for the overall process time
    END = clock();
    cpu_time_used_MAIN = ((double) (END - START)) / CLOCKS_PER_SEC;

    // cout<<"TIME - OVERALL IMAGE------------------------ : "<< cpu_time_used_MAIN<<endl;
    // cout << ".";

    return 0;
}

// This function reads an image and calculates 15 GLCM features on it
// Input: Mat object representing the image, angle and distance for GLCM calculation (optional)
// Output: Array of 15 GLCM features
double *readImage(Mat img, int angle = 0, int distance = 1){
    // Start the timer for the overall process time
    clock_t START,END;
    double cpu_time_used_MAIN;
    START = clock();

    // Define an array of GLCM features to calculate
    char *c[] = {"asm", "contrast", "corr", "var", "idm", "savg", "svar", "sentropy", "entropy", "dvar", "dentropy", "icorr1", "icorr2", "maxcorr", "sum"};

    // Start the timer for the GLCM calculation time
    clock_t start,end;
    double cpu_time_used;
    start = clock();

    // Calculate the GLCM features on the image
    double *result = calculator(img,c,15, angle, distance);

    // Release the image from memory
    img.release();

    // End the timer for the GLCM calculation time
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    // End the timer for the overall process time
    END = clock();
    cpu_time_used_MAIN = ((double) (END - START)) / CLOCKS_PER_SEC;

    // cout<<"TIME - OVERALL IMAGE------------------------ : "<< cpu_time_used_MAIN<<endl;
    // cout << ".";

    // Return the array of GLCM features
    return result;
}

// Function to read csv and get unique image names from a column
// Input: CSV file path, column name (#filename)
// Output: Vector of unique image names
vector<string> getUniqueImageNames(string csvPath, string columnName) {
    // Read the CSV file
    ifstream csv(csvPath);

    // Get the column names
    string line;
    getline(csv, line);
    stringstream ss(line);
    vector<string> columnNames;
    while (ss.good()) {
        string substr;
        getline(ss, substr, ',');
        columnNames.push_back(substr);
    }

    // Get the index of the column
    int columnIndex = -1;
    for (int i = 0; i < columnNames.size(); i++) {
        if (columnNames[i] == columnName) {
            columnIndex = i;
            break;
        }
    }

    // Get the unique image names
    vector<string> imageNames;
    while (getline(csv, line)) {
        stringstream ss(line);
        for (int i = 0; i < columnIndex; i++) {
            getline(ss, line, ',');
        }
        getline(ss, line, ',');
        imageNames.push_back(line);
    }
    sort(imageNames.begin(), imageNames.end());
    imageNames.erase(unique(imageNames.begin(), imageNames.end()), imageNames.end());

    return imageNames;
}

// Function to read image names in a folder and get image names
// Input: Image folder path
// Output: Vector of unique image names
vector<string> getImageNames(string folderPath) {
    // Read the folder
    DIR *dir;
    struct dirent *ent;
    vector<string> imageNames;

    // Open the directory and read its contents
    string fullFolderPath = folderPath;
    if ((dir = opendir(fullFolderPath.c_str())) != NULL) {
        while ((ent = readdir(dir)) != NULL) {
            // Ignore hidden files and directories
            if (ent->d_name[0] != '.') {
                imageNames.push_back(ent->d_name);
            }
        }
        closedir(dir);
    } else {
        // Print an error message if the directory cannot be opened
        perror("");
        sort(imageNames.begin(), imageNames.end());
        imageNames.erase(unique(imageNames.begin(), imageNames.end()), imageNames.end());
        return imageNames;
    }

    // Sort and remove duplicates from the vector of image names
    sort(imageNames.begin(), imageNames.end());
    imageNames.erase(unique(imageNames.begin(), imageNames.end()), imageNames.end());

    return imageNames;
}

// Function to find the nipple location of the mammogram
// Input: Image
// Output: Nipple location
Point findNipple(Mat img, string imageName) {
    // Consider the furthest right non-background pixel that was neither on the top 15% and not in the bottom 15% of the image and its coordinate as nipple.
    // The top 15% and bottom 15% of the image are ignored to avoid the case where the nipple is not visible in the image.
    int y = 0;
    int x = img.cols - 1;
    Point nipple = Point(-1, -1);
    // Loop through the rows of the image
    for (; y < img.rows; y++) {
        // Loop through the columns of the image from right to left
        for (; x >= 0; x--) {
            // If the pixel is not background and is further right than the current nipple, set it as the new nipple
            if (img.at<uchar>(y, x) > 0 && nipple.x < x) {
                nipple = Point(x, y);
                break;
            // If the pixel is not background but is not further right than the current nipple, break out of the loop
            } else if (img.at<uchar>(y, x) > 0) {
                break;
            }
        }
        // Reset the column index to the rightmost column for the next row
        x = img.cols - 1;
    }

    // Convert the grayscale image to a color image
    cvtColor(img, img, COLOR_GRAY2BGR);

    // Draw a circle at the nipple location
    circle(img, nipple, 50, Scalar(0, 0, 255), -1);

    // Save the image with the nipple circle drawn on it
    char outputImageName[200];
    sprintf(outputImageName, "%s/nipple/%s", savePath, imageName.c_str());
    cv::imwrite(outputImageName, img);

    // Return the nipple location
    return nipple;
}

// Function to isolate background
// Input: Image
Mat isolateBackground(Mat img, string imageName) {
    // Add 1 to all pixels without overflow
    for (int y = 0; y < img.rows; y++) {
        for (int x = 0; x < img.cols; x++) {
            if (img.at<uchar>(y, x) < 255) {
                img.at<uchar>(y, x) = img.at<uchar>(y, x) + 1;
            }
        }
    }

    // Find the nipple location
    int y = 0;
    int x = img.cols - 1;
    Point nipple = Point(-1, -1);
    for (; y < img.rows; y++) {
        for (; x >= 0; x--) {
            // If the pixel is not background and is further right than the current nipple, set it as the new nipple
            if (img.at<uchar>(y, x) > 0 && nipple.x < x) {
                nipple = Point(x, y);
                break;
            // If the pixel is not background but is not further right than the current nipple, break out of the loop
            } else if (img.at<uchar>(y, x) > 0) {
                break;
            // If the pixel is background, set its value to 0
            } else {
                img.at<uchar>(y, x) = 0;
            }
        }
        // Reset the column index to the rightmost column for the next row
        x = img.cols - 1;
    }

    // Save the image with the background isolated
    char outputImageName[200];
    sprintf(outputImageName, "%s/no_background/%s", savePath, imageName.c_str());
    cv::imwrite(outputImageName, img);

    // Return the image with the background isolated
    return img;
}

// This function initializes the results saved to a CSV file
void initSaveResults() {
    // Lock the mutex to ensure thread safety
    std::lock_guard<std::mutex> lock(mtxSaveResults);

    // Define an array of GLCM feature names
    char *c[] = {"asm", "contrast", "corr", "var", "idm", "savg", "svar", "sentropy", "entropy", "dvar", "dentropy", "icorr1", "icorr2", "maxcorr", "sum"};

    // Define file paths for the two CSV files to be created
    char folderPath1[300];
    char folderPath2[300];
    sprintf(folderPath1, "%s/%s.csv", SAVE_FILE_PATH, FILE_NAME);
    sprintf(folderPath2, "%s/%s2.csv", SAVE_FILE_PATH, FILE_NAME);

    // Open the first file for writing
    ofstream file;
    file.open(folderPath1, ios::out);

    // Open the second file for writing
    ofstream file2;
    file2.open(folderPath2, ios::out);

    // Write the header row for the files
    file << "imageName" << ",";
    file2 << "imageName" << ",";
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 15; j++) {
            file << c[j] << "_" << i << ",";
            file2 << c[j] << "_" << i << ",";
        }
    }
    file << "CC,";
    file2 << "CC,";
    file << "L,";
    file2 << "L,";
    file << "CM,";
    file2 << "CM,";

    file << "isCancerous" << endl;
    file2 << "isCancerous" << endl;

    // Close the files
    file.close();
    file2.close();
}

// This function saves the results of the feature extraction process to a CSV file.
// Input: Image name, array of 15 GLCM features, distance and angle for GLCM calculation, and a boolean indicating whether the image is cancerous or not.
void saveResults(string imageName, double *result, int distance, int angle, bool isCancerous) {
    // Lock the mutex to ensure thread safety
    std::lock_guard<std::mutex> lock(mtxSaveResults);

    // Define the file path for the CSV file
    char folderPath1[300];
    sprintf(folderPath1, "%s/%s.csv", SAVE_FILE_PATH, FILE_NAME);

    // Open the file for appending
    ofstream file;
    file.open(folderPath1, ios::app);

    // Write the results to the file
    file << imageName << ",";
    for (int i = 0; i < 15; i++) {
        file << result[i] << ",";
    }
    file << distance << ",";
    file << angle << ",";
    file << isCancerous << ",";
    file << endl;

    // Close the file
    file.close();
}


// This function saves the results of the feature extraction process to a CSV file.
// Input: Image name, array of 8 sets of 15 GLCM features, metadata information about the image, distance and angle for GLCM calculation, and a boolean indicating whether the image is cancerous or not.
void saveResults2(string imageName, double **resultArray, MetaData metaData, int distance, int angle, bool isCancerous) {
    
    char folderPath[300];

    // Lock the mutex to ensure thread safety
    std::lock_guard<std::mutex> lock(mtxSaveResults);

    // Open the file for appending
    ofstream file;

    // Check if the image has any non-zero GLCM features
    bool flag = false;
    for (int i = 0; i < 8; i++) {
        if ((resultArray[i])[14] != 0) {
            flag = true;
            break;
        }
    }

    // Choose the appropriate file path based on whether the image has non-zero GLCM features or not
    if (!flag) {
        sprintf(folderPath, "%s/%s2.csv", SAVE_FILE_PATH, FILE_NAME);
    } else {
        sprintf(folderPath, "%s/%s.csv", SAVE_FILE_PATH, FILE_NAME);
    }

    // Open the file for appending
    file.open(folderPath, ios::app);

    // Write the results to the file
    file << imageName << ",";

    // Write the 8 sets of 15 GLCM features to the file
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 15; j++) {
            file << (resultArray[i])[j] << ",";
        }
    }

    // Write metadata information to the file
    if (metaData.view == "CC") {
        file << "1,";
    } else {
        file << "0,";
    }

    if (metaData.side == "L") {
        file << "1,";
    } else {
        file << "0,";
    }

    if (metaData.energy == "CM") {
        file << "1,";
    } else {
        file << "0,";
    }

    // Write the boolean indicating whether the image is cancerous or not to the file
    file << isCancerous << endl;

    // Close the file
    file.close();
}

// Function to generate texture features
// Input: Image
void generateTextureFeatures(Mat img, Mat refImg, MetaData metaData, string imageName, int angle = 0, int distance = 1) {
    // Define the window size
    int windowSize = WINDOW_SIZE / OVERLAPS; // Adjust the size according to your requirements

    // Create a matrix to hold the result of feature extraction
    Mat result = Mat::zeros((img.rows - windowSize)/windowSize, (img.cols - windowSize)/windowSize, CV_32FC1);

    // Create an array to hold the feature extraction results for all images
    double *resultArray[8];

    // Create an array to hold the labels for each window
    bool labels[int((img.rows - windowSize)/windowSize) * int((img.cols - windowSize)/windowSize)];

    // Create a character array to hold the output image name
    char outputImageName[200];

    // Iterate over the image using a sliding window
    for (int y = 0; y <= img.rows - WINDOW_SIZE; y = y + windowSize) {
            
        for (int x = 0; x <= img.cols - WINDOW_SIZE; x = x + windowSize) {
            
            // Calculate the index of the current window
            int count = (int((img.cols - windowSize)/windowSize) * y/windowSize) + x/windowSize;

            // Extract the current window
            cv::Mat window = img(cv::Rect(x, y, WINDOW_SIZE, WINDOW_SIZE));
            cv::Mat refWindow = refImg(cv::Rect(x, y, WINDOW_SIZE, WINDOW_SIZE));

            // Calculate mean
            cv::Scalar mean = cv::mean(refWindow);

            // If the mean is below the threshold, set the label to false
            if (mean[0] < LABEL_THRESHOLD) {
                labels[count] = false;
                if (mean[0] > 0) {
                    continue;
                }
            } else {
                labels[count] = true;
            }
            
            // Extract features for the current window
            resultArray[0] = readImage(window, 0, 1);
            resultArray[1] = readImage(window, 0, 2);
            resultArray[2] = readImage(window, 45, 1);
            resultArray[3] = readImage(window, 45, 2);
            resultArray[4] = readImage(window, 90, 1);
            resultArray[5] = readImage(window, 90, 2);
            resultArray[6] = readImage(window, 135, 1);
            resultArray[7] = readImage(window, 135, 2);

            // Save the label and feature extraction results
            result.at<float>(y/windowSize, x/windowSize) = labels[count]*255;
            saveResults2(imageName, resultArray, metaData, distance, angle, labels[count]);
        }
    }
}

// Generate reference image with cancerous region
Mat generateReferenceImage(Mat img, string imageName, string side) {
    // Read the CSV file
    ifstream csv(CSVpath);

    // Get the column names
    string line;
    getline(csv, line);
    stringstream ss(line);
    vector<string> columnNames;
    while (ss.good()) {
        string substr;
        getline(ss, substr, ',');
        columnNames.push_back(substr);
    }

    // Get the index of the column
    int columnIndex = -1;
    for (int i = 0; i < columnNames.size(); i++) {
        if (columnNames[i] == "#filename") {
            columnIndex = i;
            break;
        }
    }

    // Get the region_shape_attributes where #filename is the current image name
    vector<string> lines;
    while (getline(csv, line)) {
        stringstream ss(line);

        for (int i = 0; i < columnIndex; i++) {
            getline(ss, line, ',');
        }
        getline(ss, line, ',');
        if (line == imageName) {
            lines.push_back(ss.str());
        }
    }

    // Create image
    cv::Mat image(img.rows, img.cols, CV_8UC3, cv::Scalar(0, 0, 0));
    // cv::Mat image = img.clone();
    // Make imagecolor
    // cvtColor(image, image, COLOR_GRAY2BGR);

    // Iterate over each line in the CSV file
    for (int i = 0; i < lines.size(); i++) {
        // Convert the line to a stringstream
        stringstream ss(lines[i]);

        // Initialize variables
        string token;
        int commaCount = 0;
        bool reachedFifthComma = false;
        string json = "";
        string jsonBuff = "";
        char target = ',';
        char cBuff = ' ';

        // Iterate over each character in the line
        for (char c : lines[i]) {
            // Check if the character matches the target
            if (c == target) {
                commaCount++;
                if (commaCount == 5) {
                    reachedFifthComma = true;
                    continue;
                }
            }
            if (reachedFifthComma) {
                // If reached the 5th comma, now parsing the JSON data
                if (c == target) {
                    // If reached another comma, finished parsing a JSON value
                    json += jsonBuff;
                    jsonBuff = cBuff;
                } else if (c == '"') {
                    // If reached a double quote, now parsing a string value
                    if (cBuff != '"') {
                        json += jsonBuff;
                        jsonBuff = cBuff;
                    }
                } else if (c == '{') {
                    // If reached an opening brace, now parsing a JSON object
                    jsonBuff = ' ';
                } else {
                    // Otherwise, still parsing a JSON value
                    jsonBuff += cBuff;
                }
            }
            cBuff = c;
        }

        // Parse the JSON data
        Json::Value root;
        Json::CharReaderBuilder builder;
        std::istringstream iss(json);
        std::string errs;
        Json::parseFromStream(builder, iss, &root, &errs);

        // Extract the shape name from the JSON data
        const Json::Value& name = root["name"];

        // Draw the shape on the image based on its name
        if (name.asString() == "polygon") {
            // Extract polygon points
            const auto& allPointsX = root["all_points_x"];
            const auto& allPointsY = root["all_points_y"];

            // Check if the number of x and y points match
            if (allPointsX.size() != allPointsY.size()) {
                std::cerr << "Mismatched point coordinates\n";
                return image;
            }

            // Draw polygon
            std::vector<cv::Point> points;
            for (int i = 0; i < allPointsX.size(); i++) {
                points.push_back(cv::Point(allPointsX[i].asInt(), allPointsY[i].asInt()));
            }

            // Draw the polygon
            cv::polylines(image, points, true, cv::Scalar(255, 255, 255), 5);
            // Fill the polygon
            cv::fillPoly(image, points, cv::Scalar(255, 255, 255));

        } else if (name.asString() == "circle") {
            // Extract circle parameters
            const auto& cx = root["cx"];
            const auto& cy = root["cy"];
            const auto& r = root["r"];

            // Draw circle
            cv::circle(image, cv::Point(cx.asInt(), cy.asInt()), r.asInt(), cv::Scalar(255, 255, 255), 2);
            // Fill the circle
            cv::circle(image, cv::Point(cx.asInt(), cy.asInt()), r.asInt(), cv::Scalar(255, 255, 255), -1);

        } else if (name.asString() == "ellipse") {
            // Extract ellipse parameters
            const auto& cx = root["cx"];
            const auto& cy = root["cy"];
            const auto& rx = root["rx"];
            const auto& ry = root["ry"];

            // Draw ellipse
            cv::ellipse(image, cv::Point(cx.asInt(), cy.asInt()), cv::Size(rx.asInt(), ry.asInt()), 0, 0, 360, cv::Scalar(255, 255, 255), 2);
            // Fill the ellipse
            cv::ellipse(image, cv::Point(cx.asInt(), cy.asInt()), cv::Size(rx.asInt(), ry.asInt()), 0, 0, 360, cv::Scalar(255, 255, 255), -1);
            
        } else if (name.asString() == "polyline") {
            // Extract polyline points
            const auto& allPointsX = root["all_points_x"];
            const auto& allPointsY = root["all_points_y"];

            // Check if the number of x and y points match
            if (allPointsX.size() != allPointsY.size()) {
                std::cerr << "Mismatched point coordinates\n";
                return image;
            }

            // Draw polyline
            std::vector<cv::Point> points;
            for (int i = 0; i < allPointsX.size(); i++) {
                points.push_back(cv::Point(allPointsX[i].asInt(), allPointsY[i].asInt()));
            }
            cv::polylines(image, points, false, cv::Scalar(255, 255, 255), 2);
            // Fill the polygon
            cv::fillPoly(image, points, cv::Scalar(255, 255, 255));
            // Select the outermost points and fill the polygon
            cv::fillConvexPoly(image, points, cv::Scalar(255, 255, 255));

        } else if (name.asString() == "point") {
            // Extract point coordinates
            const auto& cx = root["cx"];
            const auto& cy = root["cy"];

            // Draw point
            cv::circle(image, cv::Point(cx.asInt(), cy.asInt()), 2, cv::Scalar(255, 255, 255), -1);

        } else {
            // Invalid region shape
            std::cerr << "Invalid region shape: " << name.asString() << endl;
            return image;
        }
    }

    // Flip the image horizontally if it's on the right side
    if (side == "R") {
        flip(image, image, 1);
    }

    // Return the annotated image
    return image;
}
// This function reads a CSV file and saves the results to a double *result array
// Input: Path to the CSV file
// Output: Pointer to the array containing the results
double* readCSV(char* csvPath){
    // TODO: Implement the function to read the CSV file and save the results to a double *result array


}

int main() {

    // Read the image from a file
    char folderPath[300];
    sprintf(folderPath, "%s/%s", FILE_PATH, FILE_NAME);
    vector<string> files = getImageNames(folderPath);
    initSaveResults();
    int counter = 0;

    // Parallelize the loop to process multiple images at once
    #pragma omp parallel for
    for (int i = 0; i < files.size(); i++) {

        // Lock the counter to prevent race conditions
        mtxCounter.lock();
        ++counter;

        // Get the image path and metadata
        char inputImagePath[300];
        cout << "----- Processing [" << counter << "] (" << i+1 << "/" << files.size() << "): " << files[i] << endl;
        mtxCounter.unlock();

        // Split the name between two '_' and get 2nd parameter
        MetaData metaData;
        metaData.patientID = files[i].substr(0, files[i].find("_"));
        metaData.side = files[i].substr(files[i].find("_") + 1);
        metaData.energy = metaData.side.substr(metaData.side.find("_") + 1);
        metaData.view = metaData.energy.substr(metaData.energy.find("_") + 1);

        metaData.side = metaData.side.substr(0, metaData.side.find("_"));
        metaData.energy = metaData.energy.substr(0, metaData.energy.find("_"));
        metaData.view = metaData.view.substr(0, metaData.view.find("."));

        // Get the input image path based on the metadata
        if (metaData.energy == "CM") {
            sprintf(inputImagePath, "%s/%s", pathCM, files[i].c_str());
        } else if (metaData.energy == "DM") {
            sprintf(inputImagePath, "%s/%s", pathDM, files[i].c_str());
        } else {
            cout << "Error: Invalid image name" << files[i] << endl;
            continue;
        }

        // Read the input image and flip it if necessary
        Mat img = imread(inputImagePath, IMREAD_GRAYSCALE);
        if (metaData.side == "R") {
            flip(img, img, 1);
        }

        // Find the nipple location and isolate the background
        isolateBackground(img, files[i]);

        // If files[i] file exists in pathRef folder read image from there else create a black image with same size as img
        sprintf(inputImagePath, "%s/%s", pathRef, files[i].c_str());
        Mat refImg = imread(inputImagePath, IMREAD_GRAYSCALE);
        if (refImg.empty()) {
            refImg = Mat::zeros(img.rows, img.cols, CV_8UC1);
            imwrite(inputImagePath, refImg);
        }

        // Generate texture features for the image
        generateTextureFeatures(img, refImg, metaData, files[i], 0, 1);

    }
}