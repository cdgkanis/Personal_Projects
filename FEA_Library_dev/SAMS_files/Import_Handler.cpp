//
// Created by User on 1/5/2026.
//

#include "Import_Handler.h"

Geometry import_msh_2dgeometry(const std::string& filename,
                               const std::shared_ptr<Material>& Mat,
                               double Thickness) {
    std::ifstream in(filename);
    if (!in) {
        throw std::runtime_error("import_msh_2dgeometry(): Could not open mesh file: " + filename);
    }

    std::vector<std::shared_ptr<Node>> nodes;
    std::vector<std::unique_ptr<Element>> elements;
    std::unordered_map<unsigned int, std::size_t> node_id_to_index;

    std::string line;
    while (std::getline(in, line)) {
        if (line == "$Nodes") {
            std::getline(in, line);
            int nNodes = std::stoi(line);

            nodes.reserve(nNodes);

            for (int i = 0; i < nNodes; ++i) {
                std::getline(in, line);
                std::istringstream iss(line);

                unsigned int node_id_base1;
                double x, y, z;
                iss >> node_id_base1 >> x >> y >> z;

                auto node = std::make_shared<Node>();
                node->ID = static_cast<unsigned int>(nodes.size());
                node->x = {x, y, z};

                node_id_to_index[node_id_base1] = nodes.size();
                nodes.push_back(node);
            }

            std::getline(in, line); // $EndNodes
        }
        else if (line == "$Elements") {
            std::getline(in, line);
            int nElem = std::stoi(line);

            for (int i = 0; i < nElem; ++i) {
                std::getline(in, line);
                std::istringstream iss(line);

                unsigned int elem_id;
                int elem_type, num_tags;
                iss >> elem_id >> elem_type >> num_tags;

                std::vector<int> tags(num_tags);
                for (int t = 0; t < num_tags; ++t) {
                    iss >> tags[t];
                }

                if (elem_type == 2) {
                    std::vector<unsigned int> conn_gmsh(3);
                    iss >> conn_gmsh[0] >> conn_gmsh[1] >> conn_gmsh[2];

                    std::vector<unsigned int> conn(3);
                    std::vector<std::shared_ptr<Node>> elem_nodes;
                    elem_nodes.reserve(3);

                    for (int k = 0; k < 3; ++k) {
                        std::size_t idx = node_id_to_index.at(conn_gmsh[k]);
                        conn[k] = static_cast<unsigned int>(idx);
                        elem_nodes.push_back(nodes[idx]);   // same shared node
                    }

                    elements.push_back(
                        std::make_unique<Tri3>(
                            elem_id, std::move(conn), std::move(elem_nodes), Mat, Thickness
                        )
                    );
                }
                else if (elem_type == 3) {
                    std::vector<unsigned int> conn_gmsh(4);
                    iss >> conn_gmsh[0] >> conn_gmsh[1] >> conn_gmsh[2] >> conn_gmsh[3];

                    std::vector<unsigned int> conn(4);
                    std::vector<std::shared_ptr<Node>> elem_nodes;
                    elem_nodes.reserve(4);

                    for (int k = 0; k < 4; ++k) {
                        std::size_t idx = node_id_to_index.at(conn_gmsh[k]);
                        conn[k] = static_cast<unsigned int>(idx);
                        elem_nodes.push_back(nodes[idx]);   // same shared node
                    }

                    elements.push_back(
                        std::make_unique<Quad4>(
                            elem_id, std::move(conn), std::move(elem_nodes), Mat, Thickness
                        )
                    );
                }
            }

            std::getline(in, line); // $EndElements
        }
    }

    return Geometry{std::move(nodes), std::move(elements)};
}