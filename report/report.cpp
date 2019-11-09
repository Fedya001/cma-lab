#include "report.h"

#include "load.h"
#include "matrix_factory.h"
#include "print_utils.h"
#include "relaxation.h"
#include "square_matrix.h"
#include "square_matrix_manager.h"
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

} // namespace report