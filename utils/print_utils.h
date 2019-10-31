#pragma once

template<class M>
void PrintVector(std::ostream& out, const std::vector<M>& vector, const std::string& separator = " ") {
  bool first = true;
  for (const auto& element : vector) {
    if (!first) {
      out << separator;
    }
    out << element;
    first = false;
  }
}

template<class M>
void PrintMatrix(std::ostream& out, const SquareMatrix<M>& matrix,
                 bool latex_syntax_on = false, uint8_t precision = 2,
                 const std::string& label = "A") {
  out << std::fixed << std::setprecision(precision);
  out << label << " = " <<
      (latex_syntax_on ? "\\left(\\begin{array}{ *{" + std::to_string(matrix.GetDim()) + "}{c} }\n" : "{\n");

  bool first_row = true;
  for (const auto& row : matrix.data_) {
    if (!first_row) {
      out << (latex_syntax_on ? " \\\\\n" : "\n");
    }

    PrintVector(out, row, (latex_syntax_on ? " & " : "\t"));
    first_row = false;
  }
  out << (latex_syntax_on ? "\n\\end{array}\\right)" : "\n}");
}

template<class M>
void PrintColumn(std::ostream& out, const std::vector<M>& column,
                 bool latex_syntax_on = false, uint8_t precision = 2,
                 const std::string& label = "b") {
  out << std::fixed << std::setprecision(precision);
  out << label << " = "
      << (latex_syntax_on ? "\\left(\\begin{array}{c}\n" : "{\n");

  PrintVector(out, column, (latex_syntax_on ? " \\\\\n" : "\n"));
  out << (latex_syntax_on ? "\\end{array}\\right)" : "}");
}

template<class M>
void PrintRow(std::ostream& out, const std::vector<M>& row,
              uint8_t precision = 2, const std::string& label = "b") {
  out << std::fixed << std::setprecision(precision);
  out << label << " = " << "(";
  PrintVector(out, row, " ");
  out << ")";
}

