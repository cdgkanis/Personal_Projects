//
// Created by User on 1/5/2026.
//

#ifndef IMPORT_HANDLER_H
#define IMPORT_HANDLER_H
#pragma once
#include <Physics.h>
#include <fstream>
#include <string>

[[nodiscard]] Geometry import_msh_2dgeometry(const std::string& filename,const std::shared_ptr<Material>& Mat,double Thickness);

#endif //IMPORT_HANDLER_H
