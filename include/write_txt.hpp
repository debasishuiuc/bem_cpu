// write_txt.hpp

#pragma once
#include "geometry_analyzer.hpp"
#include "trgl6_icos.hpp"
#include <string>

void write_txt(const Mesh& mesh, const GeometryAnalyzer& analyzer, const std::string& outdir);
