#pragma once

#include <cstdint>
#include <string>

namespace report {
using DoubleType = double;

const inline int64_t MILLISECONDS_THRESHOLD = 80'000;
const inline int32_t SIZE_LEAP = 250;
const inline double RANDOM_GENERATOR_RADIUS = 1000.0;

// Prints all tasks solutions (in data/task_i folders)
// in specified .tex file
void LatexSolutions(const std::string& output_tex_file);

// Dumps all measurements (to build graphics in jupyter notebook)
void DumpMeasurements(const std::string& measurements_dump_file);

} // namespace report