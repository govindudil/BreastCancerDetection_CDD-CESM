require "csv"

$pivot = 0
$inputs = 123
$outputs = 1

class ProgramSynthesisCsvStore
  def initialize(filename)
    @csv = CSV.new(File.read(filename), headers: false)
  end

  def read()
    fitness_cases = []
    i = 0
    @csv.map do |row|
      if i > 0
        row = row.to_a
        fitness_case = {}
        $inputs.times do |input|
          fitness_case['i' + input.to_s] = row[input]
        end
        $outputs.times do |output|
          fitness_case['o' + output.to_s] = row[$inputs + output]
        end
        fitness_cases << fitness_case
      end
      i += 1
    end
    fitness_cases
  end
end

# training_file = File.open('TrainingSampleData.txt', 'a')
validation_file = File.open('ValidationSampleData.txt', 'a')
store = ProgramSynthesisCsvStore.new('val.csv')
fitness_cases = store.read()
cases_for_training = 0
cases_for_validation = 0
fitness_cases.each_with_index do |fitness_case, i|
  case_row = ""
  $inputs.times do |input|
    case_row = case_row + fitness_case['i' + input.to_s] + "\t"
  end
  $outputs.times do |output|
    case_row = case_row + fitness_case['o' + output.to_s]
    if output == $outputs - 1
      case_row = case_row + "\n"
    else
      case_row = case_row + "\t"
    end
  end
  random_pivot = (1..100).to_a.sample
  training_or_validation = random_pivot < $pivot ? 0 : 1
  if training_or_validation == 0
    cases_for_training += 1
    # training_file.write(case_row)
  else
    cases_for_validation += 1
    validation_file.write(case_row)
  end
end
# training_file.close
validation_file.close
# puts 'cases_for_training'
# puts cases_for_training
puts 'cases_for_validation'
puts cases_for_validation
