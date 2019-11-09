#include "report.h"

#include "load.h"
#include "matrix_factory.h"
#include "print_utils.h"
#include "relaxation.h"
#include "square_matrix.h"
#include "square_matrix_manager.h"
#include "vector_utils.h"
#include "time_measurer.h"

#include <chrono>
#include <experimental/filesystem>
#include <fstream>

namespace report {

namespace filesystem = std::experimental::filesystem;

// This functions are not accessible anywhere
void BeginDocument(std::ofstream& output) {
  output << "\\documentclass[a4paper,12pt]{article}\n"
            "\n"
            "\\usepackage[left=1cm,right=1cm,\n"
            "top=1.5cm,bottom=1.5cm,bindingoffset=0cm]{geometry}\n"
            "\\usepackage{cmap}\t\t\t\t\t% поиск в PDF\n"
            "\\usepackage[T2A]{fontenc}\t\t\t% кодировка\n"
            "\\usepackage[utf8]{inputenc}\t\t\t% кодировка исходного текста\n"
            "\n"
            "\\usepackage{natbib}\n"
            "\\usepackage{esvect}\n"
            "\\usepackage{graphicx}\n"
            "\\usepackage{amsmath}\n"
            "\\usepackage{dsfont}\n"
            "\\usepackage{amsmath}\n"
            "\n"
            "\\usepackage[russian]{babel}\n"
            "\n"
            "\\begin{document}\n"
            "\n"
            "\\author{Филипович Фёдор}\n"
            "\\title{Решение примеров из лабы №1 \\ "
            "Решение систем линейных алгебраических уравнений}\n"
            "\\date{\\today}\n"
            "\n"
            "\\maketitle\n\n";
}

void EndDocument(std::ofstream& output) {
  output << "\\end{document}\n";
}

void BeginLatexGatherEnvironment(std::ofstream& output) {
  output << "\\begin{gather*}\n";
}

void EndLatexGatherEnvironment(std::ofstream& output) {
  output << "\\end{gather*}\n";
}

void LatexNewLine(std::ofstream& output, int32_t count = 1) {
  for (int32_t i = 0; i < count; ++i) {
    output << "\\\\\n";
  }
}

void LatexSolutions(const std::string& output_tex_file) {
  if (filesystem::path(output_tex_file).extension().string() != ".tex") {
    throw std::runtime_error(".tex file must be specified");
  }
  filesystem::create_directories(filesystem::path(output_tex_file).parent_path());

  // TODO: create a singleton storage for directories
  const std::string lowdiag_task_directory("../data/task1");
  const std::string lu_systems_task_directory("../data/task2");
  const std::string sweep_method_task_directory("../data/task4");

  const std::vector<std::string> tasks_names = {
      R"(\item[\textbf{\textit{№1.}}] \textbf{\textit{Обратная матрица к lower-diagonal матрице (нижнетреугольная + диагональ).}})",
      R"(\item[\textbf{\textit{№2.}}] \textbf{\textit{LU-разложение с выбором по столбцу.}})",
      R"(\item[\textbf{\textit{№3.}}] \textbf{\textit{Сравнение производительности алгоритмов, основанных на LU  и LDLT  разложениях.}})",
      R"(\item[\textbf{\textit{№4.}}] \textbf{\textit{Решение СЛАУ методом прогонки.}})",
      R"(\item[\textbf{\textit{№5.}}] \textbf{\textit{Метод релаксации для матрицы специфичной структуры.}})"
  };

  std::string bold_italics = R"(\textbf{\textit{)";
  const int32_t precision = 3;

  std::ofstream output(output_tex_file);
  SquareMatrixManager<DoubleType> manager(SquareMatrix<DoubleType>(0));

  BeginDocument(output);
  output << "\\begin{enumerate}\n";

  // Task1
  int32_t index = 0;
  output << tasks_names[0] << "\n";
  LatexNewLine(output);
  BeginLatexGatherEnvironment(output);
  for (const auto& matrix : loader::LoadMatrices<DoubleType>(lowdiag_task_directory)) {
    manager.SetMatrix(matrix);
    std::string matrix_name = {static_cast<char>('A' + index)};
    PrintMatrix(output, matrix, true, precision, matrix_name);
    output << "\\implies ";
    if (matrix.GetDim() > 4) {
      LatexNewLine(output);
    }
    try {
      auto inverse = manager.InverseLowdiag(true);
      PrintMatrix(output, inverse, true, precision, matrix_name + "^{-1}");
    } catch (std::logic_error& ex) {
      output << "det\\ A = 0";
    }
    LatexNewLine(output, 1);

    ++index;
  }
  EndLatexGatherEnvironment(output);

  // Task2
  output << tasks_names[1] << "\n";
  BeginLatexGatherEnvironment(output);
  index = 0;
  for (const auto& system : loader::LoadSystems<DoubleType>(lu_systems_task_directory)) {
    manager.SetMatrix(system.first);
    PrintMatrix(output, system.first, true, 0, {static_cast<char>('A' + index)});
    PrintColumn(output, system.second, true, 0, {static_cast<char>('a' + index)});
    LatexNewLine(output);

    auto dlu = manager.PerformDLU(true);
    PrintMatrix(output, dlu.low_up, true, ((index == 0) ? 2 : 5), "G_{L, U}^" + std::to_string(index));
    LatexNewLine(output);
    PrintColumn(output, dlu.rows_permutations, true, 0, "P");
    output << ";\\ ";
    PrintColumn(output, SquareMatrixManager<DoubleType>::SolveSystemDLU(dlu, system.second),
                true, precision, "X");
    LatexNewLine(output, 1);

    ++index;
  }
  EndLatexGatherEnvironment(output);

  // Task4
  index = 0;
  output << tasks_names[3] << "\n";
  BeginLatexGatherEnvironment(output);
  for (const auto& system : loader::LoadThreeDiagonalSystems<DoubleType>(sweep_method_task_directory)) {
    PrintColumn(output, system.first.SolveSystem(system.second),
                true, 5, "b_" + std::to_string(index));
    ++index;
  }
  EndLatexGatherEnvironment(output);

  // Task5
  output << tasks_names[4] << "\n";
  output << std::fixed << std::setprecision(10);
  BeginLatexGatherEnvironment(output);
  for (int32_t size : {5, 10, 100, 500, 1000, 2000, 4000}) {
    output << "X_{" + std::to_string(size) + '}' << " = \\left(\\begin{array}{c}\n";
    auto solution = SolveSystem<long double>(size, 0.42);
    output << std::get<0>(solution) << " \\\\\n" << std::get<1>(solution) << " \\\\\n" << "\\ldots \\\\\n"
           << std::get<1>(solution) << "\\\\\n" << std::get<2>(solution) << "\n\\end{array}\\right)";
    LatexNewLine(output);
  }
  EndLatexGatherEnvironment(output);

  output << "\\end{enumerate}\n";
  EndDocument(output);
}

void DumpMeasurements(const std::string& measurements_dump_file) {
  filesystem::create_directories(filesystem::path(measurements_dump_file).parent_path());

  const std::vector<std::string> task_names = {
      "Task1. Lowdiag inverse",
      "Task2. LU decomposition",
      "Task3. Compare LU and LDL^T",
      "Task4. Sweep method",
      "Task5. Relaxation method"
  };

  std::ofstream measurements_dump(measurements_dump_file);
  std::vector<int64_t> measurements;

  TimeMeasurer<std::chrono::milliseconds> measurer;
  MatrixFactory<DoubleType> matrix_factory
      (-RANDOM_GENERATOR_RADIUS, RANDOM_GENERATOR_RADIUS);
  SquareMatrixManager<DoubleType> manager(SquareMatrix<DoubleType>(0));

  // Explore tasks (i.e. measure execution time and other characteristics)
  // TODO: test with long double and float types

  // Task1. Lowdiag inverse
  {
    int64_t execution_time;
    int32_t size = SIZE_LEAP;

    std::cerr << task_names[0] << std::endl;
    measurements_dump << task_names[0] << std::endl;
    do {
      manager.SetMatrix(matrix_factory.CreateLowdiagMatrix(size));
      measurer.UpdateTimeStamp();
      (void) manager.InverseLowdiag();
      measurements.push_back(execution_time = measurer.GetElapsedTime());
      size += SIZE_LEAP;
      std::cerr << size << " : " << execution_time << std::endl;
    } while (execution_time < MILLISECONDS_THRESHOLD);
    PrintVector(measurements_dump, measurements, ", ");
  }

  // Task2. LU decomposition
  {
    std::cerr << task_names[1] << std::endl;
    measurements_dump << task_names[1] << std::endl;

    auto system = loader::LoadSystem<DoubleType>("../data/task2/sampleB.data");
    manager.SetMatrix(system.first);
    auto original_solution = manager.SolveSystem(system.second);

    measurements_dump << "Original solution: ";
    PrintRow(measurements_dump, original_solution);
    measurements_dump << std::endl;

    measurements_dump << std::fixed << std::setprecision(6);
    for (DoubleType noise_radius : {0.001, 0.01, 0.1, 0.5, 1.0, 5.0, 10.0, 20.0}) {
      measurements_dump << "noise_radius = " << noise_radius << " => ";
      auto noise_vector = NoiseVector(system.second, noise_radius);
      measurements_dump << "noise_vector_norm = " << GetDiffEuclidNorm(system.second, noise_vector);
      auto solution = manager.SolveSystem(noise_vector);
      measurements_dump << " => solution_noise = " << GetDiffEuclidNorm(original_solution, solution) << std::endl;
    }
  }

  // Task3. Compare LU and LDL^T
  {
    std::cerr << task_names[2] << std::endl;
    measurements_dump << task_names[2] << std::endl;

    measurements.clear();
    std::vector<int64_t> symmetric_measurements;
    for (int32_t size = 100; size <= 2000; size += 100) {
      auto matrix = matrix_factory.CreateRandomMatrix(size, true);
      manager.SetMatrix(matrix);
      auto vector = matrix_factory.CreateRandomVector(size);

      measurer.UpdateTimeStamp();
      (void) manager.SolveSystem(
          vector, SquareMatrixManager<DoubleType>::MatrixType::ORDINARY);
      measurements.push_back(measurer.GetElapsedTime());

      measurer.UpdateTimeStamp();
      (void) manager.SolveSystem(
          vector, SquareMatrixManager<DoubleType>::MatrixType::SYMMETRIC);
      symmetric_measurements.push_back(measurer.GetElapsedTime());
    }
    PrintVector(measurements_dump, measurements, ", ");
    PrintVector(measurements_dump, symmetric_measurements, ", ");
  }

  // Task5. Relaxation method
  {
    const int32_t testing_size = 100'000'000;

    std::cerr << task_names[4] << std::endl;
    measurements_dump << task_names[4] << std::endl;

    std::vector<int32_t> itetaions_number;
    std::vector<long double> omega_values;

    long double omega = 0.0;
    while (omega <= 2.0) {
      omega_values.push_back(omega);
      itetaions_number.push_back(std::get<3>(SolveSystem(testing_size, omega)));
      omega += 0.2;
    }
    PrintVector(measurements_dump, omega_values, ", ");
    measurements_dump << std::endl;
    PrintVector(measurements_dump, itetaions_number, ", ");
  }
}

} // namespace report