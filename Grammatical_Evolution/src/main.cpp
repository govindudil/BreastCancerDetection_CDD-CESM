#ifndef _MAIN_CPP_
#define _MAIN_CPP_

#include <fstream>
#include <grace/grace.hpp>
#include <iostream>
#include <chrono>
#include <array>

int fold_counter = 1;
int gen_counter = 0;
double best_score = 3;

class CancerPCAECBEvaluator : public Evaluator {
public:
  bool evaluate(Population &population) {
    // Measure time taken
    auto timeStart = std::chrono::high_resolution_clock::now();

    // Initialize variables
    best_score = 3;
    auto bestIndividual = std::dynamic_pointer_cast<NSGA2Genome>(population.individuals.front());

    // Open file to write out phenotypes
    std::ofstream file;
    std::string middle_filename = "tmp/cancerMiddle_" + std::to_string(fold_counter) + ".cpp";
    file.open(middle_filename);
    file.clear();
    bool firstRun = true;
    int populationCounter = 0;

    // Write out phenotypes to file
    for (auto individual : population.individuals) {
      auto castIndividual = std::dynamic_pointer_cast<NSGA2Genome>(individual);

      // Check if individual is valid
      if (castIndividual->isPhenotypeValid == false) {
        for (int i = 0; i < castIndividual->objectives.size(); ++i) {
          castIndividual->objectives[i] = 1;
        }
        continue;
      }

      if (!firstRun) {
        file << "," << std::endl;
      } else {
        firstRun = false;
      }
      file << "    [](double** Dim, int i) { return " << castIndividual->phenotype << "; }";
      ++populationCounter;
    }
    file.close();

    // Write out phenotype to file
    std::ifstream start;
    std::ifstream middle;
    std::ifstream end;
    
    start.open("template/cancerStart.c");
    middle.open(middle_filename);
    end.open("template/cancerEnd.c");

    std::string filename = "tmp/evaluate_population_" + std::to_string(fold_counter) + ".cpp";
    file.open(filename);
    file.clear();

    std::string contents;
    while (getline(start, contents)) {
      file << contents << std::endl;
    }

    file << populationCounter << "> expressions = {" << std::endl;

    while (getline(middle, contents)) {
      file << contents << std::endl;
    }

    while (getline(end, contents)) {
      file << contents << std::endl;
    }

    file.close();
    start.close();
    middle.close();
    end.close();

    auto checkPoint1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = checkPoint1 - timeStart;
    std::cout << "Checkpoint 1: " << elapsed.count() << " s\n";
    // Compile and run the files
    std::string resultfile = "tmp/result_" + std::to_string(fold_counter);
    std::string command = "g++ " + filename + " obj/GEcancer.o -lm -o tmp/a_" + std::to_string(fold_counter) + ".out && tmp/a_" + std::to_string(fold_counter) + ".out " + std::to_string(fold_counter) + " > " + resultfile;
    if (system(command.c_str()) == -1) {
      std::cerr << "Compilation or execution failed.\n";
      std::cerr << "Execution aborted.\n";
      exit(0);
    }

    // Measure time taken
    auto checkPoint2 = std::chrono::high_resolution_clock::now();
    elapsed = checkPoint2 - timeStart;
    std::cout << "Checkpoint 2: " << elapsed.count() << " s\n";
    elapsed = checkPoint2 - checkPoint1;
    std::cout << "Time taken to analyse fitness: " << elapsed.count() << " s\n";

    // Read in results and update objectives
    std::ifstream result;
    result.open(resultfile);
    for (auto individual : population.individuals) {
      auto castIndividual = std::dynamic_pointer_cast<NSGA2Genome>(individual);

      // Check if individual is valid
      if (castIndividual->isPhenotypeValid == false) {
        continue;
      }

      int i = 0;
      double sum = 0;
      for (int j = 0; j < 3; j++) {
        getline(result, contents);
        castIndividual->objectives[i] = stod(contents);
        sum += castIndividual->objectives[i] * castIndividual->objectives[i];
        ++i;
      }

      // Update best individual
      if (sum < best_score) {
        best_score = sum;
        bestIndividual = std::dynamic_pointer_cast<NSGA2Genome>(individual);
      }
    }
    result.close();

    // Measure time taken
    auto checkPoint3 = std::chrono::high_resolution_clock::now();
    elapsed = checkPoint3 - timeStart;
    std::cout << "Checkpoint 3: " << elapsed.count() << " s\n";
    // Validate with best individual
    if (bestIndividual->isPhenotypeValid == true) {
          std::ifstream startVal;
          std::ifstream endVal;
          std::ofstream fileVal;

          // Open the validation start and end templates
          startVal.open("template/validateStart.c");
          endVal.open("template/validateEnd.c");

          // Open the temporary file for validation
          fileVal.open("tmp/individual_val.cpp");
          fileVal.clear();

          std::string contentsVal;

          // Write the contents of the validation start template to the temporary file
          while (getline(startVal, contentsVal)) {
            fileVal << contentsVal << std::endl;
          }

          // Write the phenotype of the best individual to the temporary file
          fileVal << bestIndividual->phenotype << std::endl;

          // Write the contents of the validation end template to the temporary file
          while (getline(endVal, contentsVal)) {
            fileVal << contentsVal << std::endl;
          }

          // Close the validation files
          fileVal.close();
          startVal.close();
          endVal.close();

          // Open the training start and end templates
          startVal.open("template/trainStart.c");
          endVal.open("template/trainEnd.c");

          // Open the temporary file for training
          fileVal.open("tmp/individual_train.cpp");
          fileVal.clear();

          // Write the contents of the training start template to the temporary file
          while (getline(startVal, contentsVal)) {
            fileVal << contentsVal << std::endl;
          }

          // Write the phenotype of the best individual to the temporary file
          fileVal << bestIndividual->phenotype << std::endl;

          // Write the contents of the training end template to the temporary file
          while (getline(endVal, contentsVal)) {
            fileVal << contentsVal << std::endl;
          }

          // Close the training files
          fileVal.close();
          startVal.close();
          endVal.close();

          // Measure time taken
          auto checkPoint4 = std::chrono::high_resolution_clock::now();
          elapsed = checkPoint4 - timeStart;
          std::cout << "Checkpoint 4: " << elapsed.count() << " s\n";

          // Compile and run the files for training data
          std::string commandTra = "echo 0 > tmp/result_train && g++ tmp/individual_train.cpp obj/GEcancer.o -lm -o tmp/a_train.out && tmp/a_train.out " + std::to_string(fold_counter) + " > tmp/result_train";
          if (system(commandTra.c_str()) == -1) {
            std::cerr << "Compilation or execution failed.\n";
            std::cerr << "Execution aborted.\n";
            exit(0);
          }

          // Read back scores
          std::ifstream resultTra;
          resultTra.open("tmp/result_train");

          // Open the log file for the best individual against training data
          std::ofstream bestIndvidualLog;
          std::string filename = "log/best_train_indvidual_" + std::to_string(fold_counter) + ".csv";
          bestIndvidualLog.open(filename, std::ios::app);

          // Write the generation number and scores to the log file
          bestIndvidualLog << gen_counter;
          while (getline(resultTra, contentsVal)) {
            bestIndvidualLog << ", " << contentsVal;
          }
          bestIndvidualLog << ", \"" << bestIndividual->phenotype << "\"" << std::endl;

          // Close the log file
          bestIndvidualLog.close();

          // Measure time taken
          auto checkPoint5 = std::chrono::high_resolution_clock::now();
          elapsed = checkPoint5 - timeStart;
          std::cout << "Checkpoint 5: " << elapsed.count() << " s\n";
          elapsed = checkPoint5 - checkPoint4;
          std::cout << "Time taken to analyse the best individual against train data: " << elapsed.count() << " s\n";

          // Compile and run the files for validation data
          std::string commandVal = "echo 0 > tmp/result_val && g++ tmp/individual_val.cpp obj/GEcancer.o -lm -o tmp/a_val.out && tmp/a_val.out " + std::to_string(fold_counter) + " > tmp/result_val";
          if (system(commandVal.c_str()) == -1) {
            std::cerr << "Compilation or execution failed.\n";
            std::cerr << "Execution aborted.\n";
            exit(0);
          }

          // Read back scores
          std::ifstream resultVal;
          resultVal.open("tmp/result_val");

          // Open the log file for the best individual against validation data
          filename = "log/best_val_indvidual_" + std::to_string(fold_counter) + ".csv";
          bestIndvidualLog.open(filename, std::ios::app);

          // Write the generation number and scores to the log file
          bestIndvidualLog << gen_counter;
          while (getline(resultVal, contentsVal)) {
            bestIndvidualLog << ", " << contentsVal;
          }
          bestIndvidualLog << ", \"" << bestIndividual->phenotype << "\"" << std::endl;

          // Close the log file
          bestIndvidualLog.close();

          // Measure time taken
          auto checkPoint6 = std::chrono::high_resolution_clock::now();
          elapsed = checkPoint6 - timeStart;
          std::cout << "Checkpoint 6: " << elapsed.count() << " s\n";
          elapsed = checkPoint6 - checkPoint5;
          std::cout << "Time taken to analyse the best individual against validation data: " << elapsed.count() << " s\n";
          resultVal.close();

        } else {
          // If the best individual is invalid, log the results as such
          std::ofstream bestIndvidualLog;
          std::string filename = "log/best_val_indvidual_" + std::to_string(fold_counter) + ".csv";
          bestIndvidualLog.open(filename, std::ios::app);

          // Write the generation number and "Invalid Individual" to the log file
          bestIndvidualLog << gen_counter;
          for (int i = 0; i < 11; ++i) {
            bestIndvidualLog << ", " << bestIndividual->objectives[i];
          }
          bestIndvidualLog << ", " << "Invalid Individual" << std::endl;

          // Close the log file
          bestIndvidualLog.close();

          filename = "log/best_train_indvidual_" + std::to_string(fold_counter) + ".csv";
          bestIndvidualLog.open(filename, std::ios::app);

          // Write the generation number and "Invalid Individual" to the log file
          bestIndvidualLog << gen_counter;
          for (int i = 0; i < 11; ++i) {
            bestIndvidualLog << ", " << bestIndividual->objectives[i];
          }
          bestIndvidualLog << ", " << "Invalid Individual" << std::endl;

          // Close the log file
          bestIndvidualLog.close();
        }

        // Increment the generation counter
        ++gen_counter;

        // Measure time taken
        auto finish = std::chrono::high_resolution_clock::now();
        elapsed = finish - timeStart;
        std::cout << "Time taken to analyse a generation: " << elapsed.count() << " s\n";

        // Return true to indicate successful analysis of the generation
        return true;
      }
};

// main function
int main(int argc, char **argv) {
  // create a backup folder with a unique name
  if (system("unique_folder=\"$(date +'\%Y\%m\%d_\%H\%M')_$(uuidgen | cut -d'-' -f1)\"; mkdir log.bak/\"$unique_folder\"; mv log/* log.bak/\"$unique_folder\"/") == -1) {
    std::cerr << "Compilation or execution failed.\n";
    std::cerr << "Execution aborted.\n";
    exit(0);
  }

  // Iterate five times
  for (fold_counter = 1; fold_counter < 6; ++fold_counter) {
    gen_counter = 0;
    best_score = 3;
    // create and open log files for best individual
    std::ofstream bestIndvidualLog;
    std::string filename = "log/best_train_indvidual_" + std::to_string(fold_counter) + ".csv";
    bestIndvidualLog.open(filename);
    bestIndvidualLog << "Generation, FPR, 1 - TPR, 1 - AUC, F1, BA, MCC, Accuracy, TP, FP, TN, FN, Phenotype" << std::endl;
    bestIndvidualLog.close();
    filename = "log/best_val_indvidual_" + std::to_string(fold_counter) + ".csv";
    bestIndvidualLog.open(filename);
    bestIndvidualLog << "Generation, FPR, 1 - TPR, 1 - AUC, F1, BA, MCC, Accuracy, TP, FP, TN, FN, Phenotype" << std::endl;
    bestIndvidualLog.close();

    // Create the Genetic Algorithm with NSGA2 Algorithm
    NSGA2Algorithm<NSGA2Population, GEInitialiser, GEMapper,
                  CancerPCAECBEvaluator, NSGA2Selection, GECrossover, GEMutation,
                  NSGA2Replacement, GETermination, NSGA2Statistics>
        ga(argc, argv, "./settings.ini", 3);

    // Run the evolution
    ga.evolve();

    // create a folder for the current fold and move the log files and statistics to it
    std::string command = "mkdir log/fold_" + std::to_string(fold_counter) + "; mv log/generation* log/fold_" + std::to_string(fold_counter) + "; mv log/best_train_indvidual_" + std::to_string(fold_counter) + ".csv log/best_val_indvidual_" + std::to_string(fold_counter) + ".csv log/statistics.csv log/fold_" + std::to_string(fold_counter);
    if (system(command.c_str()) == -1) {
      std::cerr << "Compilation or execution failed.\n";
      std::cerr << "Execution aborted.\n";
      exit(0);
    }
  }

  // Return Success
  return EXIT_SUCCESS;
}

#endif
