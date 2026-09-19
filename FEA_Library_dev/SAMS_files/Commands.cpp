//
// Created by User on 1/5/2026.
//

#include "Commands.h"

#include <vector>

namespace fs = std::filesystem;

// Helper functions
nlohmann::json materialToJson(const Material& material){
    if (material.type() == MaterialType::Isotropic) {
        const auto* iso = dynamic_cast<const IsotropicMaterial*>(&material);
        if (!iso) {
            throw std::runtime_error("materialToJson(): Material type says Isotropic, but cast failed.");
        }

        return json{
                    {"name", iso->name()},
                    {"type", "isotropic"},
                    {"rho",  iso->rho()},
                    {"E",    iso->E()},
                    {"nu",   iso->nu()}
        };
    }

    if (material.type() == MaterialType::Orthotropic) {
        const auto* ortho = dynamic_cast<const OrthotropicMaterial*>(&material);
        if (!ortho) {
            throw std::runtime_error("materialToJson(): Material type says Orthotropic, but cast failed.");
        }

        return json{
                    {"name", ortho->name()},
                    {"type", "orthotropic"},
                    {"rho",  ortho->rho()},
                    {"Ex",   ortho->Ex()},
                    {"Ey",   ortho->Ey()},
                    {"Ez",   ortho->Ez()},
                    {"nuxy", ortho->nuxy()},
                    {"nuyz", ortho->nuyz()},
                    {"nuxz", ortho->nuxz()},
                    {"Gxy",  ortho->Gxy()},
                    {"Gyz",  ortho->Gyz()},
                    {"Gxz",  ortho->Gxz()}
        };
    }

    throw std::runtime_error("materialToJson(): Unsupported material type.");
}

static double readDouble(const std::string& prompt){
    double value;
    while (true) {
        std::cout << prompt;
        if (std::cin >> value) {
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            return value;
        }

        std::cout << "readDouble(): Invalid number. Try again.\n";
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
}

static std::string readString(const std::string& prompt){
    std::string value;
    do {
        std::cout << prompt;
        std::getline(std::cin, value);
    } while (value.empty());

    return value;
}

static MaterialType readMatType(){
    int choice;

    while (true) {
        std::cout << "Select material type:\n";
        std::cout << "1. Isotropic\n";
        std::cout << "2. Orthotropic\n";
        std::cout << "Choice: ";

        if (std::cin >> choice) {
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            if (choice == 1) return MaterialType::Isotropic;
            if (choice == 2) return MaterialType::Orthotropic;
        }

        std::cout << "readMatType(): Invalid choice. Try again.\n";
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
}

static Eigen::MatrixXd reduceMatrix(const Eigen::MatrixXd& A, const std::vector<int>& free_dofs) {
    Eigen::MatrixXd R(free_dofs.size(), free_dofs.size());
    for (int i = 0; i < static_cast<int>(free_dofs.size()); ++i) {
        for (int j = 0; j < static_cast<int>(free_dofs.size()); ++j) {
            R(i, j) = A(free_dofs[i], free_dofs[j]);
        }
    }
    return R;
}

static Eigen::VectorXd reduceVector(const Eigen::VectorXd& v, const std::vector<int>& free_dofs) {
    Eigen::VectorXd r(free_dofs.size());
    for (int i = 0; i < static_cast<int>(free_dofs.size()); ++i) {
        r(i) = v(free_dofs[i]);
    }
    return r;
}

static Eigen::VectorXd expandVector(const Eigen::VectorXd& reduced,const std::vector<int>& free_dofs,int total_dofs) {
    Eigen::VectorXd full = Eigen::VectorXd::Zero(total_dofs);
    for (int i = 0; i < static_cast<int>(free_dofs.size()); ++i) {
        full(free_dofs[i]) = reduced(i);
    }
    return full;
}

// Commands
std::shared_ptr<Material> loadMaterial(const std::string& materialName,const std::string& filePath){
    std::ifstream file(filePath);
    if (!file) {
        throw std::runtime_error("loadMaterial(): Could not open material library: " + filePath);
    }

    json data;
    file >> data;

    if (!data.contains("materials") || !data["materials"].is_array()) {
        throw std::runtime_error("loadMaterial(): JSON file must contain a 'materials' array.");
    }

    for (const auto& mat : data["materials"]) {
        if (!mat.contains("name") || !mat.contains("type") || !mat.contains("rho")) {
            continue;
        }

        if (mat["name"].get<std::string>() != materialName) {
            continue;
        }

        const auto type = mat["type"].get<std::string>();
        const auto rho = mat["rho"].get<double>();

        if (type == "isotropic") {
            if (!mat.contains("E") || !mat.contains("nu")) {
                throw std::runtime_error("loadMaterial(): Isotropic material is missing E or nu.");
            }

            return std::make_shared<IsotropicMaterial>(
                mat["name"].get<std::string>(),
                rho,
                mat["E"].get<double>(),
                mat["nu"].get<double>()
            );
        }

        if (type == "orthotropic") {
            const char* keys[] = {
                "Ex", "Ey", "Ez",
                "nuxy", "nuyz", "nuxz",
                "Gxy", "Gyz", "Gxz"
            };

            for (const auto* key : keys) {
                if (!mat.contains(key)) {
                    throw std::runtime_error(
                        "loadMaterial(): Orthotropic material is missing field: " + std::string(key)
                    );
                }
            }

            return std::make_shared<OrthotropicMaterial>(
                mat["name"].get<std::string>(),
                rho,
                mat["Ex"].get<double>(),
                mat["Ey"].get<double>(),
                mat["Ez"].get<double>(),
                mat["nuxy"].get<double>(),
                mat["nuyz"].get<double>(),
                mat["nuxz"].get<double>(),
                mat["Gxy"].get<double>(),
                mat["Gyz"].get<double>(),
                mat["Gxz"].get<double>()
            );
        }

        throw std::runtime_error("loadMaterial(): Unknown material type: " + type);
    }

    throw std::runtime_error("loadMaterial(): Material not found: " + materialName);
}

void addMaterialToLibrary(const Material& material,const std::string& filePath,bool overwrite){
    json data;

    {
        if (std::ifstream in(filePath); in.good()) {
            in >> data;
        } else {
            data = json{{"materials", json::array()}};
        }
    }

    if (!data.contains("materials") || !data["materials"].is_array()) {
        data["materials"] = json::array();
    }

    json newMat = materialToJson(material);

    auto& materials = data["materials"];
    const auto newName = newMat["name"].get<std::string>();

    for (auto& entry : materials) {
        if (entry.contains("name") && entry["name"].get<std::string>() == newName) {
            if (!overwrite) {
                throw std::runtime_error("addMaterialToLibrary(): Material already exists in library: " + newName);
            }
            entry = newMat;
            std::ofstream out(filePath);
            if (!out) {
                throw std::runtime_error("addMaterialToLibrary(): Could not write material library: " + filePath);
            }
            out << std::setw(4) << data << std::endl;
            return;
        }
    }

    materials.push_back(newMat);

    std::ofstream out(filePath);
    if (!out) {
        throw std::runtime_error("addMaterialToLibrary(): Could not write material library: " + filePath);
    }

    out << std::setw(4) << data << std::endl;
}

std::shared_ptr<Material> createMaterialFromUserInput(){
    MaterialType type = readMatType();

    std::string name = readString("Material name: ");
    double rho = readDouble("Density rho [kg/m^3]: ");

    if (type == MaterialType::Isotropic) {
        double E  = readDouble("Young's modulus E [Pa]: ");
        double nu = readDouble("Poisson ratio nu [-]: ");

        return std::make_unique<IsotropicMaterial>(name, rho, E, nu);
    }

    double Ex   = readDouble("Ex [Pa]: ");
    double Ey   = readDouble("Ey [Pa]: ");
    double Ez   = readDouble("Ez [Pa]: ");
    double nuxy = readDouble("nuxy [-]: ");
    double nuyz = readDouble("nuyz [-]: ");
    double nuxz = readDouble("nuxz [-]: ");
    double Gxy  = readDouble("Gxy [Pa]: ");
    double Gyz  = readDouble("Gyz [Pa]: ");
    double Gxz  = readDouble("Gxz [Pa]: ");

    return std::make_shared<OrthotropicMaterial>(
        name, rho,
        Ex, Ey, Ez,
        nuxy, nuyz, nuxz,
        Gxy, Gyz, Gxz
    );
}

void Project::assignMaterial(const std::shared_ptr<Material>& material) {
    if (!material) {
        throw std::invalid_argument("Project::assignMaterial(): null material.");
    }
    defaultMaterial = material;
    state = ProjectState::MaterialAssigned;
}

void Project::importGeometry(const std::string& mshFile,const std::shared_ptr<Material>& material,double thickness,unsigned int dofsPerNode,double meshDensity) {
    if (!material) {
        throw std::invalid_argument("Project::importGeometry(): null material.");
    }

    Geometry geom = import_msh_2dgeometry(mshFile, material, thickness);
    geometryFile = mshFile;
    state = ProjectState::GeometryImported;

    Mesh builtMesh = Mesh::buildMesh(std::move(geom),static_cast<int>(dofsPerNode),meshDensity, material, thickness);

    mesh = std::make_unique<Mesh>(std::move(builtMesh));
    state = ProjectState::MeshReady;
}

void Project::addNodalLoad(unsigned int node_id, const std::vector<double>& values) const {
    if (!mesh) {
        throw std::runtime_error("Project::addNodalLoad(): mesh not available.");
    }
    mesh->addNodalLoad(node_id, values);
}

void Project::setNodalLoad(unsigned int node_id, const std::vector<double>& values) const {
    if (!mesh) {
        throw std::runtime_error("Project::setNodalLoad(): mesh not available.");
    }
    mesh->setNodalLoad(node_id, values);
}

void Project::printSummary(std::ostream& os) const {
    os << "Project name : " << projectName << '\n';
    os << "Geometry file: " << geometryFile << '\n';
    os << "Mesh loaded  : " << (mesh ? "yes" : "no") << '\n';

    if (mesh) {
        os << "Nodes        : " << mesh->nodes.size() << '\n';
        os << "Elements     : " << mesh->elements.size() << '\n';
        os << "Total DOFs   : " << mesh->total_dofs << '\n';
    }
}

void Project::setupElementAssumptions(const std::vector<std::string>& assumptions) const {
    if (!mesh) {
        throw std::runtime_error("Project::setupElementAssumptions(): mesh not assigned.");
    }

    for (const auto& elem : mesh->elements) {
        for (const std::string& ass : assumptions) {
            elem->AddAssumption(ass);
        }
    }
}

void Project::removeElementAssumptions(AssumptionType assumption) const {
    if (!mesh) {
        throw std::runtime_error("Project::removeElementAssumptions(): mesh not assigned.");
    }

    for (const auto& elem : mesh->elements) {
        elem->RemoveAssumption(assumption);
    }
}

// Eigen Frequencies Small Deformations Protject =======================================================================
void EigenFrequencyProject::setNumModes(int n) {
    if (n <= 0) {
        throw std::invalid_argument("EigenFrequencyProject::setNumModes(): number of modes must be positive.");
    }
    numModes = n;
}

void EigenFrequencyProject::assemble() {
    if (!mesh) {
        throw std::runtime_error("EigenFrequencyProject::assemble(): mesh not assigned.");
    }

    mesh->globalStiffnessMatrix();
    mesh->globalMassMatrix();
}

void EigenFrequencyProject::solve() {
    if (!mesh) {
        throw std::runtime_error("EigenFrequencyProject::solve(): mesh not assigned.");
    }

    if (mesh->K.rows() == 0 || mesh->K.cols() == 0) {
        throw std::runtime_error("EigenFrequencyProject::solve(): stiffness matrix K is empty. Call assemble() first.");
    }

    if (mesh->M.rows() == 0 || mesh->M.cols() == 0) {
        throw std::runtime_error("EigenFrequencyProject::solve(): mass matrix M is empty. Call assemble() first.");
    }

    std::vector<int> fixed_dofs = mesh->getSupportDofs();

    Eigen::MatrixXd K_used = mesh->K;
    Eigen::MatrixXd M_used = mesh->M;

    std::vector<int> free_dofs;
    free_dofs.reserve(mesh->total_dofs);

    if (!fixed_dofs.empty()) {
        std::sort(fixed_dofs.begin(), fixed_dofs.end());
        fixed_dofs.erase(std::unique(fixed_dofs.begin(), fixed_dofs.end()), fixed_dofs.end());

        std::vector<bool> is_fixed(mesh->total_dofs, false);
        for (int dof : fixed_dofs) {
            if (dof < 0 || dof >= static_cast<int>(mesh->total_dofs)) {
                throw std::out_of_range("EigenFrequencyProject::solve(): fixed DOF out of range.");
            }
            is_fixed[dof] = true;
        }

        for (int i = 0; i < static_cast<int>(mesh->total_dofs); ++i) {
            if (!is_fixed[i]) {
                free_dofs.push_back(i);
            }
        }

        if (free_dofs.empty()) {
            throw std::runtime_error("EigenFrequencyProject::solve(): all DOFs are constrained.");
        }

        K_used.resize(free_dofs.size(), free_dofs.size());
        M_used.resize(free_dofs.size(), free_dofs.size());

        for (int i = 0; i < static_cast<int>(free_dofs.size()); ++i) {
            for (int j = 0; j < static_cast<int>(free_dofs.size()); ++j) {
                K_used(i, j) = mesh->K(free_dofs[i], free_dofs[j]);
                M_used(i, j) = mesh->M(free_dofs[i], free_dofs[j]);
            }
        }
    } else {
        for (int i = 0; i < static_cast<int>(mesh->total_dofs); ++i) {
            free_dofs.push_back(i);
        }
    }

    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> solver(K_used, M_used);

    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("EigenFrequencyProject::solve(): eigenvalue solver failed.");
    }

    const Eigen::VectorXd& lambda_all = solver.eigenvalues();
    const Eigen::MatrixXd& modes_reduced = solver.eigenvectors();

    std::vector<int> valid_ids;
    valid_ids.reserve(lambda_all.size());

    for (int i = 0; i < lambda_all.size(); ++i) {
        if (lambda_all(i) > 0.0) {
            valid_ids.push_back(i);
        }
    }

    if (valid_ids.empty()) {
        throw std::runtime_error("EigenFrequencyProject::solve(): no positive eigenvalues found.");
    }

    const int nm = std::min(numModes, static_cast<int>(valid_ids.size()));

    omega.resize(nm);
    frequenciesHz.resize(nm);
    modes = Eigen::MatrixXd::Zero(mesh->total_dofs, nm);

    for (int i = 0; i < nm; ++i) {
        const int id = valid_ids[i];

        omega(i) = std::sqrt(lambda_all(id));
        frequenciesHz(i) = omega(i) / (2.0 * M_PI);

        for (int r = 0; r < static_cast<int>(free_dofs.size()); ++r) {
            modes(free_dofs[r], i) = modes_reduced(r, id);
        }
    }

    state = ProjectState::Completed;
}

std::vector<double> EigenFrequencyProject::getEigenFrequencies() const {
    std::vector<double> frequencies(frequenciesHz.size());
    for (int i = 0; i < frequenciesHz.size(); ++i) {
        frequencies[i] = frequenciesHz(i);
    }
    return frequencies;
}

void EigenFrequencyProject::exportModeToCSV(const std::string& filePath, int modeIndex) const {
    if (!mesh) {
        throw std::runtime_error("EigenFrequencyProject::exportModeToCSV(): mesh not assigned.");
    }

    if (modes.cols() == 0 || frequenciesHz.size() == 0) {
        throw std::runtime_error("EigenFrequencyProject::exportModeToCSV(): no modal results available.");
    }

    if (modeIndex < 0 || modeIndex >= modes.cols()) {
        throw std::out_of_range("EigenFrequencyProject::exportModeToCSV(): mode index out of range.");
    }

    if (mesh->node_dofs.empty()) {
        throw std::runtime_error("EigenFrequencyProject::exportModeToCSV(): node DOF map is empty.");
    }

    std::ofstream out(filePath);
    if (!out) {
        throw std::runtime_error("EigenFrequencyProject::exportModeToCSV(): cannot open file: " + filePath);
    }

    out << std::setprecision(16);

    out << "# mode_index," << modeIndex << '\n';
    out << "# frequency_hz," << frequenciesHz(modeIndex) << '\n';
    out << "node_id,x,y,ux,uy,umag\n";

    const Eigen::VectorXd mode = modes.col(modeIndex);

    for (std::size_t i = 0; i < mesh->nodes.size(); ++i) {
        const auto& node = mesh->nodes[i];
        const auto& dofs = mesh->node_dofs[i];

        if (dofs.size() < 2) {
            throw std::runtime_error("EigenFrequencyProject::exportModeToCSV(): expected at least 2 DOFs per node.");
        }

        const double x = node->x[0];
        const double y = node->x[1];
        const double ux = mode(dofs[0]);
        const double uy = mode(dofs[1]);
        const double umag = std::sqrt(ux * ux + uy * uy);

        out << node->ID << ","
            << x << ","
            << y << ","
            << ux << ","
            << uy << ","
            << umag << '\n';
    }
}

void EigenFrequencyProject::exportModesToCSV(const std::string& folderPath, int nModes) const {
    if (nModes <= 0) {
        throw std::invalid_argument("EigenFrequencyProject::exportModesToCSV(): nModes must be positive.");
    }

    std::filesystem::create_directories(folderPath);

    const int available = static_cast<int>(modes.cols());
    const int n = std::min(nModes, available);

    for (int i = 0; i < n; ++i) {
        std::ostringstream filePath;
        filePath << folderPath << "\\mode_" << (i + 1) << ".csv";
        exportModeToCSV(filePath.str(), i);
    }
}

void EigenFrequencyProject::plotModes(const std::string& folderPath, int nModes) {
    namespace fs = std::filesystem;

    const fs::path pythonExe =
        R"(C:\Users\User\Desktop\8o_e3amhno\Analisis_of_machenical_structurs\SAMS\venv\Scripts\python.exe)";
    const fs::path scriptPath =
        R"(C:\Users\User\Desktop\8o_e3amhno\Analisis_of_machenical_structurs\SAMS\plot_mode.py)";
    const fs::path folder(folderPath);
    const fs::path connCsv = folder / "connectivity.csv";

    if (folderPath.empty()) {
        throw std::runtime_error("EigenFrequencyProject::plotModes(): folder path is empty.");
    }
    if (nModes <= 0) {
        throw std::invalid_argument("EigenFrequencyProject::plotModes(): nModes must be positive.");
    }
    if (!fs::exists(pythonExe)) {
        throw std::runtime_error("EigenFrequencyProject::plotModes(): python executable not found.");
    }
    if (!fs::exists(scriptPath)) {
        throw std::runtime_error("EigenFrequencyProject::plotModes(): plot_mode.py not found.");
    }
    if (!fs::exists(folder) || !fs::is_directory(folder)) {
        throw std::runtime_error("EigenFrequencyProject::plotModes(): results folder not found.");
    }
    if (!fs::exists(connCsv)) {
        throw std::runtime_error(
            "EigenFrequencyProject::plotModes(): connectivity file not found: " + connCsv.string()
        );
    }

    for (int i = 1; i <= nModes; ++i) {
        fs::path csvFile = folder / ("mode_" + std::to_string(i) + ".csv");

        if (!fs::exists(csvFile)) {
            throw std::runtime_error(
                "EigenFrequencyProject::plotModes(): missing file: " + csvFile.string()
            );
        }

        std::ostringstream cmd;
        cmd << "cmd /c \""
            << "\"" << pythonExe.string() << "\" "
            << "\"" << scriptPath.string() << "\" "
            << "\"" << csvFile.string() << "\" "
            << "\"" << connCsv.string() << "\""
            << "\"";

        std::cout << "Command: " << cmd.str() << '\n';

        int rc = std::system(cmd.str().c_str());

        if (rc != 0) {
            throw std::runtime_error(
                "EigenFrequencyProject::plotModes(): failed while plotting " +
                csvFile.string() + ". Return code = " + std::to_string(rc)
            );
        }
    }
}

void EigenFrequencyProject::exportConnectivityToCSV(const std::string& path) const {
    namespace fs = std::filesystem;

    if (!mesh) {
        throw std::runtime_error(
            "EigenFrequencyProject::exportConnectivityToCSV(): mesh not assigned."
        );
    }

    if (path.empty()) {
        throw std::invalid_argument(
            "EigenFrequencyProject::exportConnectivityToCSV(): path is empty."
        );
    }

    fs::path outPath(path);

    // If user passed a directory, place connectivity.csv inside it
    if (!outPath.has_extension()) {
        fs::create_directories(outPath);
        outPath /= "connectivity.csv";
    }
    // If user passed a file path, ensure parent directory exists
    else if (outPath.has_parent_path()) {
        fs::create_directories(outPath.parent_path());
    }

    std::ofstream out(outPath);
    if (!out) {
        throw std::runtime_error(
            "EigenFrequencyProject::exportConnectivityToCSV(): cannot open file: " +
            outPath.string()
        );
    }

    out << "elem_id,type,n1,n2,n3,n4\n";

    for (std::size_t e = 0; e < mesh->elements.size(); ++e) {
        const auto& elem = mesh->elements[e];
        const auto& conn = elem->connectivity();

        if (elem->name() == "Tri3") {
            if (conn.size() != 3) {
                throw std::runtime_error(
                    "exportConnectivityToCSV(): Tri3 element does not have 3 nodes."
                );
            }

            out << e << ",Tri3,"
                << conn[0] << ","
                << conn[1] << ","
                << conn[2] << ",-1\n";
        }
        else if (elem->name() == "Quad4") {
            if (conn.size() != 4) {
                throw std::runtime_error(
                    "exportConnectivityToCSV(): Quad4 element does not have 4 nodes."
                );
            }

            out << e << ",Quad4,"
                << conn[0] << ","
                << conn[1] << ","
                << conn[2] << ","
                << conn[3] << "\n";
        }
        else {
            throw std::runtime_error(
                "exportConnectivityToCSV(): unsupported element type: " + elem->name()
            );
        }
    }
}

void EigenFrequencyProject::printSummary(std::ostream& os) const {
    Project::printSummary(os);

    os << "Analysis type : EigenFrequencies\n";
    os << "Requested modes: " << numModes << '\n';
    os << "State         : ";

    switch (state) {
        case ProjectState::Created:            os << "Created"; break;
        case ProjectState::MaterialAssigned:   os << "MaterialAssigned"; break;
        case ProjectState::GeometryImported:   os << "GeometryImported"; break;
        case ProjectState::MeshReady:          os << "MeshReady"; break;
        case ProjectState::MatricesAssembled:  os << "MatricesAssembled"; break;
        case ProjectState::LoadsAssigned:      os << "LoadsAssigned"; break;
        case ProjectState::Ready:              os << "Ready"; break;
        case ProjectState::Completed:          os << "Completed"; break;
    }
    os << '\n';

    if (frequenciesHz.size() > 0) {
        os << "Computed frequencies [Hz]:\n";
        for (int i = 0; i < frequenciesHz.size(); ++i) {
            os << "  Mode " << i + 1
               << " : " << frequenciesHz(i) << '\n';
        }
    } else {
        os << "Computed frequencies [Hz]: none\n";
    }
}

// Transient Small Deformations Protject ===============================================================================
void TransientProject::setTimeSpan(double t_start, double t_final, double time_step) {
    if (t_final <= t_start) {
        throw std::invalid_argument(
            "TransientProject::setTimeSpan(): t_final must be greater than t_start."
        );
    }

    if (time_step <= 0.0) {
        throw std::invalid_argument(
            "TransientProject::setTimeSpan(): time_step must be positive."
        );
    }

    t0 = t_start;
    tf = t_final;
    dt = time_step;
}

void TransientProject::assemble() {
    if (!mesh) {
        throw std::runtime_error("TransientProject::assemble(): mesh not assigned.");
    }

    if(!mesh->K.size()){
        mesh->globalStiffnessMatrix();
    }
    if(!mesh->M.size()){
        mesh->globalMassMatrix();
    }

    setState(ProjectState::MatricesAssembled);

    if (mesh->C.rows() == 0 || mesh->C.cols() == 0) {
        mesh->C = Eigen::MatrixXd::Zero(mesh->total_dofs, mesh->total_dofs);
    }

    if (mesh->F.size() == 0) {
        mesh->F = Eigen::VectorXd::Zero(mesh->total_dofs);
    }

    setState(ProjectState::Ready);
}

void TransientProject::solve() {
    if (!mesh) {
        throw std::runtime_error("TransientProject::solve(): No mesh assigned.");
    }
    if (state != ProjectState::Ready) {
        throw std::runtime_error("TransientProject::solve(): Project not ready.");
    }
    if (dt <= 0.0 || tf <= t0) {
        throw std::runtime_error("TransientProject::solve(): Invalid time range.");
    }

    const int ndof_full = static_cast<int>(mesh->total_dofs);

    if (mesh->M.rows() != ndof_full || mesh->M.cols() != ndof_full ||
        mesh->C.rows() != ndof_full || mesh->C.cols() != ndof_full ||
        mesh->K.rows() != ndof_full || mesh->K.cols() != ndof_full) {
        throw std::runtime_error("TransientProject::solve(): Full system matrices have wrong dimensions.");
    }

    // -------------------------------------------------------------------------------------------------
    // Build free DOFs
    std::vector<int> fixed_dofs = mesh->getSupportDofs();

    std::vector<bool> is_fixed(ndof_full, false);
    for (int dof : fixed_dofs) {
        if (dof < 0 || dof >= ndof_full) {
            throw std::out_of_range("TransientProject::solve(): fixed DOF out of range.");
        }
        is_fixed[dof] = true;
    }

    std::vector<int> free_dofs;
    free_dofs.reserve(ndof_full);
    for (int i = 0; i < ndof_full; ++i) {
        if (!is_fixed[i]) {
            free_dofs.push_back(i);
        }
    }

    if (free_dofs.empty()) {
        throw std::runtime_error("TransientProject::solve(): all DOFs are constrained.");
    }

    const int ndof = static_cast<int>(free_dofs.size());

    // -------------------------------------------------------------------------------------------------
    // Reduced matrices
    Eigen::MatrixXd M = reduceMatrix(mesh->M, free_dofs);
    Eigen::MatrixXd C = reduceMatrix(mesh->C, free_dofs);
    Eigen::MatrixXd K = reduceMatrix(mesh->K, free_dofs);

    if (M.rows() != ndof || M.cols() != ndof ||
        C.rows() != ndof || C.cols() != ndof ||
        K.rows() != ndof || K.cols() != ndof) {
        throw std::runtime_error("TransientProject::solve(): Reduced matrix dimensions are invalid.");
    }

    // -------------------------------------------------------------------------------------------------
    // Cache variable-load DOFs in FULL coordinates
    std::vector<std::pair<std::function<Eigen::VectorXd(double)>, std::vector<int>>> activeLoads;
    activeLoads.reserve(loadFunctions.size());

    for (const auto& [f, x_target] : loadFunctions) {
        if (x_target.size() < 2) {
            throw std::runtime_error("TransientProject::solve(): Variable load target must contain at least x,y.");
        }

        int node_id = mesh->findNodeByCoordinates(x_target[0], x_target[1], 1e-9);
        if (node_id < 0) {
            throw std::runtime_error("TransientProject::solve(): No node found for variable load.");
        }

        const std::vector<int>& dofs = mesh->node_dofs[node_id];
        Eigen::VectorXd test = f(t0);

        if (test.size() != static_cast<int>(dofs.size())) {
            throw std::runtime_error("TransientProject::solve(): Load size does not match node DOFs.");
        }

        activeLoads.emplace_back(f, dofs);
    }

    // -------------------------------------------------------------------------------------------------
    // Build reduced force vector
    auto buildForce = [&](double t) -> Eigen::VectorXd {
        Eigen::VectorXd F_full = Eigen::VectorXd::Zero(ndof_full);

        if (mesh->F.size() == ndof_full) {
            F_full += mesh->F;
        }

        for (const auto& [f, dofs] : activeLoads) {
            Eigen::VectorXd f_local = f(t);
            for (int i = 0; i < f_local.size(); ++i) {
                F_full(dofs[i]) += f_local(i);
            }
        }

        return reduceVector(F_full, free_dofs);
    };

    // -------------------------------------------------------------------------------------------------
    // Time storage
    const int nSteps = static_cast<int>(std::floor((tf - t0) / dt)) + 1;

    time.clear();
    U_hist.clear();
    V_hist.clear();
    A_hist.clear();

    time.reserve(nSteps);
    U_hist.reserve(nSteps);
    V_hist.reserve(nSteps);
    A_hist.reserve(nSteps);

    // -------------------------------------------------------------------------------------------------
    // Initial conditions in reduced system
    Eigen::VectorXd u =
        (U0.size() == ndof_full) ? reduceVector(U0, free_dofs)
                                 : Eigen::VectorXd::Zero(ndof);

    Eigen::VectorXd v =
        (V0.size() == ndof_full) ? reduceVector(V0, free_dofs)
                                 : Eigen::VectorXd::Zero(ndof);

    Eigen::LDLT<Eigen::MatrixXd> Msolver(M);
    if (Msolver.info() != Eigen::Success) {
        throw std::runtime_error("TransientProject::solve(): Failed to factorize reduced mass matrix.");
    }

    Eigen::VectorXd F0 = buildForce(t0);
    Eigen::VectorXd a = Msolver.solve(F0 - C * v - K * u);

    if (Msolver.info() != Eigen::Success) {
        throw std::runtime_error("TransientProject::solve(): Failed to solve initial acceleration.");
    }

    if (!u.allFinite() || !v.allFinite() || !a.allFinite()) {
        throw std::runtime_error("TransientProject::solve(): Initial conditions produced non-finite values.");
    }

    time.push_back(t0);
    U_hist.push_back(expandVector(u, free_dofs, ndof_full));
    V_hist.push_back(expandVector(v, free_dofs, ndof_full));
    A_hist.push_back(expandVector(a, free_dofs, ndof_full));

    // -------------------------------------------------------------------------------------------------
    switch (method) {
        case Solver::Neumark: {
            const double b = beta;
            const double g = gamma;

            if (b <= 0.0) {
                throw std::runtime_error("TransientProject::solve(): Newmark beta must be > 0.");
            }

            const double c0 = 1.0 / (b * dt * dt);
            const double c1 = g / (b * dt);
            const double c2 = 1.0 / (b * dt);
            const double c3 = 1.0 / (2.0 * b) - 1.0;
            const double c4 = g / b - 1.0;
            const double c5 = dt * (g / (2.0 * b) - 1.0);

            Eigen::MatrixXd Keff = K + c0 * M + c1 * C;
            Eigen::LDLT<Eigen::MatrixXd> KeffSolver(Keff);

            if (KeffSolver.info() != Eigen::Success) {
                throw std::runtime_error("TransientProject::solve(): Failed to factorize effective stiffness matrix.");
            }

            for (int n = 0; n < nSteps - 1; ++n) {
                const double t_np1 = t0 + (n + 1) * dt;
                Eigen::VectorXd Fnp1 = buildForce(t_np1);

                Eigen::VectorXd Feff =
                    Fnp1
                    + M * (c0 * u + c2 * v + c3 * a)
                    + C * (c1 * u + c4 * v + c5 * a);

                Eigen::VectorXd u_new = KeffSolver.solve(Feff);

                if (KeffSolver.info() != Eigen::Success) {
                    throw std::runtime_error("TransientProject::solve(): Failed to solve Newmark step.");
                }

                Eigen::VectorXd a_new = c0 * (u_new - u) - c2 * v - c3 * a;
                Eigen::VectorXd v_new = v + dt * ((1.0 - g) * a + g * a_new);

                if (!u_new.allFinite() || !v_new.allFinite() || !a_new.allFinite()) {
                    throw std::runtime_error("TransientProject::solve(): Newmark diverged.");
                }

                time.push_back(t_np1);
                U_hist.push_back(expandVector(u_new, free_dofs, ndof_full));
                V_hist.push_back(expandVector(v_new, free_dofs, ndof_full));
                A_hist.push_back(expandVector(a_new, free_dofs, ndof_full));

                u = u_new;
                v = v_new;
                a = a_new;
            }
            break;
        }

        case Solver::CentralDifs: {
            Eigen::MatrixXd Aeff = (1.0 / (dt * dt)) * M + (1.0 / (2.0 * dt)) * C;
            Eigen::LDLT<Eigen::MatrixXd> AeffSolver(Aeff);

            if (AeffSolver.info() != Eigen::Success) {
                throw std::runtime_error("TransientProject::solve(): Failed to factorize central-difference matrix.");
            }

            Eigen::VectorXd u_minus_1 = u - dt * v + 0.5 * dt * dt * a;

            for (int n = 0; n < nSteps - 1; ++n) {
                const double t_n = t0 + n * dt;
                const double t_np1 = t0 + (n + 1) * dt;

                Eigen::VectorXd Fn = buildForce(t_n);

                Eigen::VectorXd rhs =
                    Fn
                    - (K - 2.0 / (dt * dt) * M) * u
                    - ((1.0 / (dt * dt)) * M - (1.0 / (2.0 * dt)) * C) * u_minus_1;

                Eigen::VectorXd u_new = AeffSolver.solve(rhs);

                if (AeffSolver.info() != Eigen::Success) {
                    throw std::runtime_error("TransientProject::solve(): Failed to solve central-difference step.");
                }

                Eigen::VectorXd v_new = (u_new - u_minus_1) / (2.0 * dt);
                Eigen::VectorXd a_new = (u_new - 2.0 * u + u_minus_1) / (dt * dt);

                if (!u_new.allFinite() || !v_new.allFinite() || !a_new.allFinite()) {
                    throw std::runtime_error("TransientProject::solve(): Central differences diverged.");
                }

                time.push_back(t_np1);
                U_hist.push_back(expandVector(u_new, free_dofs, ndof_full));
                V_hist.push_back(expandVector(v_new, free_dofs, ndof_full));
                A_hist.push_back(expandVector(a_new, free_dofs, ndof_full));

                u_minus_1 = u;
                u = u_new;
                v = v_new;
                a = a_new;
            }
            break;
        }

        default:
            throw std::runtime_error("TransientProject::solve(): Unknown transient solver.");
    }

    if (U_hist.empty()) {
        throw std::runtime_error("TransientProject::solve(): No displacement history was generated.");
    }

    mesh->U = U_hist.back();
    setState(ProjectState::Completed);
}

void TransientProject::printSummary(std::ostream& os) const {
    os << "Project summary\n";
    os << "Name         : " << projectName << '\n';
    os << "Type         : TransientAnalysis\n";

    os << "State        : ";
    switch (state) {
        case ProjectState::Created:            os << "Created"; break;
        case ProjectState::MaterialAssigned:   os << "MaterialAssigned"; break;
        case ProjectState::GeometryImported:   os << "GeometryImported"; break;
        case ProjectState::MeshReady:          os << "MeshReady"; break;
        case ProjectState::MatricesAssembled:  os << "MatricesAssembled"; break;
        case ProjectState::LoadsAssigned:      os << "LoadsAssigned"; break;
        case ProjectState::Ready:              os << "Ready"; break;
        case ProjectState::Completed:          os << "Completed"; break;
        default:                               os << "Unknown"; break;
    }
    os << '\n';

    os << "t0           : " << t0 << '\n';
    os << "tf           : " << tf << '\n';
    os << "dt           : " << dt << '\n';
    os << "steps        : " << time.size() << '\n';

    if (!U_hist.empty()) {
        os << "||U_final||  : " << U_hist.back().norm() << '\n';
    }
    if (!V_hist.empty()) {
        os << "||V_final||  : " << V_hist.back().norm() << '\n';
    }
    if (!A_hist.empty()) {
        os << "||A_final||  : " << A_hist.back().norm() << '\n';
    }
}

void TransientProject::setVariableLoad(const std::function<Eigen::VectorXd(double)>& f,const std::vector<double>& x) {
    if (x.size() != 2) {
        throw std::runtime_error("TransientProject::setVariableLoad(): coordinates must be 2D.");
    }

    loadFunctions.emplace_back(f, x);
}

void TransientProject::exportResultsToCSV(const std::string& folderPath,
                                          const std::string& fileName) const {
    namespace fs = std::filesystem;

    if (time.empty() || U_hist.empty() || V_hist.empty() || A_hist.empty()) {
        throw std::runtime_error("TransientProject::exportResultsToCSV(): no transient results available.");
    }

    if (U_hist.size() != time.size() || V_hist.size() != time.size() || A_hist.size() != time.size()) {
        throw std::runtime_error("TransientProject::exportResultsToCSV(): history sizes do not match time vector.");
    }

    if (folderPath.empty()) {
        throw std::runtime_error("TransientProject::exportResultsToCSV(): folder path is empty.");
    }

    fs::path folder(folderPath);
    fs::create_directories(folder);

    fs::path filePath = folder / fileName;
    const int ndof = static_cast<int>(U_hist[0].size());

    std::ofstream out(filePath);
    if (!out) {
        throw std::runtime_error("TransientProject::exportResultsToCSV(): cannot open file: " + filePath.string());
    }

    out << std::setprecision(16);

    out << "time";
    for (int i = 0; i < ndof; ++i) out << ",u_" << i;
    for (int i = 0; i < ndof; ++i) out << ",v_" << i;
    for (int i = 0; i < ndof; ++i) out << ",a_" << i;
    out << '\n';

    for (std::size_t k = 0; k < time.size(); ++k) {
        if (U_hist[k].size() != ndof || V_hist[k].size() != ndof || A_hist[k].size() != ndof) {
            throw std::runtime_error("TransientProject::exportResultsToCSV(): inconsistent DOF count in history.");
        }

        out << time[k];
        for (int i = 0; i < ndof; ++i) out << ',' << U_hist[k](i);
        for (int i = 0; i < ndof; ++i) out << ',' << V_hist[k](i);
        for (int i = 0; i < ndof; ++i) out << ',' << A_hist[k](i);
        out << '\n';
    }
}

void TransientProject::exportNodesToCSV(const std::string& folderPath,
                                        const std::string& fileName) const {
    namespace fs = std::filesystem;

    if (!mesh) {
        throw std::runtime_error("TransientProject::exportNodesToCSV(): no mesh assigned.");
    }

    if (folderPath.empty()) {
        throw std::runtime_error("TransientProject::exportNodesToCSV(): folder path is empty.");
    }

    fs::path folder(folderPath);
    fs::create_directories(folder);

    fs::path filePath = folder / fileName;

    std::ofstream out(filePath);
    if (!out) {
        throw std::runtime_error("TransientProject::exportNodesToCSV(): cannot open file: " + filePath.string());
    }

    out << std::setprecision(16);
    out << "node_id,x,y\n";

    for (std::size_t i = 0; i < mesh->nodes.size(); ++i) {
        const auto& node = mesh->nodes[i];
        if (!node || node->x.size() < 2) {
            throw std::runtime_error("TransientProject::exportNodesToCSV(): invalid node coordinates.");
        }

        out << node->ID << ',' << node->x[0] << ',' << node->x[1] << '\n';
    }
}

void TransientProject::exportConnectivityToCSV(const std::string& folderPath,
                                               const std::string& fileName) const {
    namespace fs = std::filesystem;

    if (!mesh) {
        throw std::runtime_error("TransientProject::exportConnectivityToCSV(): no mesh assigned.");
    }

    if (folderPath.empty()) {
        throw std::runtime_error("TransientProject::exportConnectivityToCSV(): folder path is empty.");
    }

    fs::path folder(folderPath);
    fs::create_directories(folder);

    fs::path filePath = folder / fileName;

    std::ofstream out(filePath);
    if (!out) {
        throw std::runtime_error("TransientProject::exportConnectivityToCSV(): cannot open file: " + filePath.string());
    }

    out << "element_id,n1,n2\n";

    std::size_t edge_id = 0;

    for (std::size_t e = 0; e < mesh->elements.size(); ++e) {
        const auto& elem = mesh->elements[e];
        const auto& conn = elem->connectivity();

        if (conn.size() == 2) {
            out << edge_id++ << ',' << conn[0] << ',' << conn[1] << '\n';
        }
        else if (conn.size() == 3) {
            out << edge_id++ << ',' << conn[0] << ',' << conn[1] << '\n';
            out << edge_id++ << ',' << conn[1] << ',' << conn[2] << '\n';
            out << edge_id++ << ',' << conn[2] << ',' << conn[0] << '\n';
        }
        else if (conn.size() == 4) {
            out << edge_id++ << ',' << conn[0] << ',' << conn[1] << '\n';
            out << edge_id++ << ',' << conn[1] << ',' << conn[2] << '\n';
            out << edge_id++ << ',' << conn[2] << ',' << conn[3] << '\n';
            out << edge_id++ << ',' << conn[3] << ',' << conn[0] << '\n';
        }
        else {
            throw std::runtime_error(
                "TransientProject::exportConnectivityToCSV(): unsupported element with " +
                std::to_string(conn.size()) + " nodes."
            );
        }
    }
}

void TransientProject::plotResults(const std::string& folderPath, int dof, const std::string& fileName, bool overwrite) const {
    namespace fs = std::filesystem;

    const fs::path pythonExe =
        R"(C:\Users\User\Desktop\8o_e3amhno\Analisis_of_machenical_structurs\SAMS\venv\Scripts\python.exe)";
    const fs::path scriptPath =
        R"(C:\Users\User\Desktop\8o_e3amhno\Analisis_of_machenical_structurs\SAMS\cmake-build-debug\plot_transient.py)";

    if (folderPath.empty()) {
        throw std::runtime_error("TransientProject::plotResults(): folder path is empty.");
    }
    if (dof < 0) {
        throw std::runtime_error("TransientProject::plotResults(): dof must be >= 0.");
    }
    if (!fs::exists(pythonExe)) {
        throw std::runtime_error("TransientProject::plotResults(): python executable not found: " + pythonExe.string());
    }
    if (!fs::exists(scriptPath)) {
        throw std::runtime_error("TransientProject::plotResults(): plot script not found: " + scriptPath.string());
    }

    fs::path folder(folderPath);
    fs::create_directories(folder);

    fs::path csvPath = folder / fileName;

    if (overwrite || !fs::exists(csvPath)) {
        exportResultsToCSV(folderPath, fileName);
    }

    std::ostringstream cmd;
    cmd << "cmd /c "
        << '\"'
        << '\"' << pythonExe.string() << '\"'
        << " "
        << '\"' << scriptPath.string() << '\"'
        << " "
        << '\"' << csvPath.string() << '\"'
        << " --dof " << dof
        << '\"';

    std::cout << "Command: " << cmd.str() << '\n';

    const int rc = std::system(cmd.str().c_str());
    if (rc != 0) {
        throw std::runtime_error(
            "TransientProject::plotResults(): failed to run python script. Return code = " +
            std::to_string(rc)
        );
    }
}

void TransientProject::animateSolution(const std::string& folderPath, double scale, int fps, const std::string& resultsFileName,
                                       const std::string& nodesFileName, const std::string& connectivityFileName,
                                       bool overwriteResults, bool overwriteGeometry, bool saveGif) const {
    namespace fs = std::filesystem;

    const fs::path pythonExe =
        R"(C:\Users\User\Desktop\8o_e3amhno\Analisis_of_machenical_structurs\SAMS\venv\Scripts\python.exe)";
    const fs::path scriptPath =
        R"(C:\Users\User\Desktop\8o_e3amhno\Analisis_of_machenical_structurs\SAMS\animate_transient.py)";

    if (folderPath.empty()) {
        throw std::runtime_error("TransientProject::animateSolution(): folder path is empty.");
    }
    if (scale <= 0.0) {
        throw std::runtime_error("TransientProject::animateSolution(): scale must be positive.");
    }
    if (fps <= 0) {
        throw std::runtime_error("TransientProject::animateSolution(): fps must be positive.");
    }
    if (!mesh) {
        throw std::runtime_error("TransientProject::animateSolution(): no mesh assigned.");
    }
    if (!fs::exists(pythonExe)) {
        throw std::runtime_error("TransientProject::animateSolution(): python executable not found: " + pythonExe.string());
    }
    if (!fs::exists(scriptPath)) {
        throw std::runtime_error("TransientProject::animateSolution(): animation script not found: " + scriptPath.string());
    }

    fs::path folder(folderPath);
    fs::create_directories(folder);

    fs::path resultsPath = folder / resultsFileName;
    fs::path nodesPath = folder / nodesFileName;
    fs::path connPath = folder / connectivityFileName;

    if (overwriteResults || !fs::exists(resultsPath)) {
        exportResultsToCSV(folderPath, resultsFileName);
    }

    if (overwriteGeometry || !fs::exists(nodesPath)) {
        exportNodesToCSV(folderPath, nodesFileName);
    }

    if (overwriteGeometry || !fs::exists(connPath)) {
        exportConnectivityToCSV(folderPath, connectivityFileName);
    }

    std::ostringstream cmd;
    cmd << "cmd /c "
        << '\"'
        << '\"' << pythonExe.string() << '\"'
        << " "
        << '\"' << scriptPath.string() << '\"'
        << " "
        << '\"' << folder.string() << '\"'
        << " --results " << '\"' << resultsFileName << '\"'
        << " --nodes " << '\"' << nodesFileName << '\"'
        << " --connectivity " << '\"' << connectivityFileName << '\"'
        << " --scale " << scale
        << " --fps " << fps;

    if (saveGif) {
        cmd << " --save-gif";
    }

    cmd << '\"';

    std::cout << "Command: " << cmd.str() << '\n';

    const int rc = std::system(cmd.str().c_str());
    if (rc != 0) {
        throw std::runtime_error(
            "TransientProject::animateSolution(): failed to run python script. Return code = " +
            std::to_string(rc)
        );
    }
}

// Static Small Deformations Protject ==================================================================================
void StaticProject::assemble() {
    if (!mesh) {
        throw std::runtime_error("StaticProject::assemble(): mesh not assigned.");
    }

    const int ndof = static_cast<int>(mesh->total_dofs);

    if (mesh->K.rows() != ndof || mesh->K.cols() != ndof || mesh->K.norm() == 0.0) {
        mesh->globalStiffnessMatrix();
    }

    if (mesh->K.rows() != ndof || mesh->K.cols() != ndof) {
        throw std::runtime_error("StaticProject::assemble(): global stiffness matrix has wrong dimensions.");
    }

    if (mesh->F.size() != ndof) {
        mesh->F = Eigen::VectorXd::Zero(ndof);
    }

    setState(ProjectState::MatricesAssembled);
    setState(ProjectState::Ready);
}

void StaticProject::setLoad(const std::vector<double>& F_ext, const std::vector<double>& x) {
    if (!mesh) {
        throw std::runtime_error("StaticProject::setLoad(): no mesh assigned.");
    }

    if (x.size() != 2) {
        throw std::runtime_error("StaticProject::setLoad(): coordinates must be 2D.");
    }

    const double tol = 1e-9;
    int nodeIndex = -1;

    for (std::size_t i = 0; i < mesh->nodes.size(); ++i) {
        const auto& node = mesh->nodes[i];
        if (!node || node->x.size() < 2) {
            continue;
        }

        if (std::abs(node->x[0] - x[0]) <= tol &&
            std::abs(node->x[1] - x[1]) <= tol) {
            nodeIndex = static_cast<int>(i);
            break;
        }
    }

    if (nodeIndex < 0) {
        throw std::runtime_error("StaticProject::setLoad(): no node found at given coordinates.");
    }

    const auto& node = mesh->nodes[nodeIndex];
    const std::vector<int>& dofs = mesh->node_dofs[nodeIndex];

    if (F_ext.size() != dofs.size()) {
        throw std::runtime_error("StaticProject::setLoad(): load size does not match node DOFs.");
    }

    if (mesh->F.size() == 0) {
        mesh->F = Eigen::VectorXd::Zero(mesh->total_dofs);
    }

    if (mesh->F.size() != mesh->total_dofs) {
        throw std::runtime_error("StaticProject::setLoad(): global load vector has wrong size.");
    }

    for (std::size_t i = 0; i < F_ext.size(); ++i) {
        mesh->F(dofs[i]) += F_ext[i];
    }

    std::cout << "StaticProject::setLoad()\n";
    std::cout << "  node index : " << nodeIndex << '\n';
    std::cout << "  node ID    : " << node->ID << '\n';
    std::cout << "  coords     : (" << node->x[0] << ", " << node->x[1] << ")\n";
    std::cout << "  dofs       : ";
    for (int dof : dofs) std::cout << dof << ' ';
    std::cout << "\n  load       : ";
    for (double v : F_ext) std::cout << v << ' ';
    std::cout << "\n  ||F||      : " << mesh->F.norm() << '\n';

    if (mesh->K.size()) {
        setState(ProjectState::Ready);
    } else {
        setState(ProjectState::LoadsAssigned);
    }
}

void StaticProject::solve() {
    if (!mesh) {
        throw std::runtime_error("StaticProject::solve(): no mesh assigned.");
    }

    if (state != ProjectState::Ready && state != ProjectState::LoadsAssigned) {
        throw std::runtime_error("StaticProject::solve(): project not ready.");
    }

    const int ndof_full = static_cast<int>(mesh->total_dofs);

    if (mesh->K.rows() != ndof_full || mesh->K.cols() != ndof_full) {
        throw std::runtime_error("StaticProject::solve(): stiffness matrix has wrong dimensions.");
    }

    if (mesh->F.size() != ndof_full) {
        throw std::runtime_error("StaticProject::solve(): global load vector has wrong dimensions.");
    }

    std::vector<int> fixed_dofs = mesh->getSupportDofs();
    std::vector<bool> is_fixed(ndof_full, false);

    for (int dof : fixed_dofs) {
        if (dof < 0 || dof >= ndof_full) {
            throw std::out_of_range("StaticProject::solve(): fixed DOF out of range.");
        }
        is_fixed[dof] = true;
    }

    std::vector<int> free_dofs;
    free_dofs.reserve(ndof_full);
    for (int i = 0; i < ndof_full; ++i) {
        if (!is_fixed[i]) {
            free_dofs.push_back(i);
        }
    }

    if (free_dofs.empty()) {
        throw std::runtime_error("StaticProject::solve(): all DOFs are constrained.");
    }

    Eigen::MatrixXd Kff = reduceMatrix(mesh->K, free_dofs);
    Eigen::VectorXd Ff  = reduceVector(mesh->F, free_dofs);

    if (Ff.norm() < 1e-12) {
        throw std::runtime_error(
            "StaticProject::solve(): reduced load vector is zero. "
            "No effective external load is applied on free DOFs."
        );
    }

    Eigen::FullPivLU<Eigen::MatrixXd> lu(Kff);
    const int rank = lu.rank();
    std::cout << "  rank(Kff)  = " << rank << " / " << Kff.rows() << '\n';

    if (!lu.isInvertible()) {
        throw std::runtime_error(
            "StaticProject::solve(): reduced stiffness matrix is singular or not invertible."
        );
    }

    Eigen::LDLT<Eigen::MatrixXd> solver(Kff);
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("StaticProject::solve(): failed to factorize reduced stiffness matrix.");
    }

    const double minAbsD = solver.vectorD().cwiseAbs().minCoeff();
    std::cout << "  min|D|     = " << minAbsD << '\n';

    if (minAbsD < 1e-14) {
        throw std::runtime_error(
            "StaticProject::solve(): reduced stiffness matrix is singular or nearly singular."
        );
    }

    Eigen::VectorXd Uf = solver.solve(Ff);
    if (solver.info() != Eigen::Success || !Uf.allFinite()) {
        throw std::runtime_error("StaticProject::solve(): failed to solve static system.");
    }

    U = expandVector(Uf, free_dofs, ndof_full);
    mesh->U = U;

    reactions = mesh->K * U - mesh->F;

    for (int dof : free_dofs) {
        reactions(dof) = 0.0;
    }

    if (!U.allFinite() || !reactions.allFinite()) {
        throw std::runtime_error("StaticProject::solve(): non-finite solution or reactions.");
    }

    setState(ProjectState::Completed);
}

void StaticProject::exportResultsToCSV(const std::string& folderPath, const std::string& fileName) const {
    namespace fs = std::filesystem;

    if (!mesh) {
        throw std::runtime_error("StaticProject::exportResultsToCSV(): no mesh assigned.");
    }

    if (U.size() == 0) {
        throw std::runtime_error("StaticProject::exportResultsToCSV(): no static results available.");
    }

    if (reactions.size() != U.size()) {
        throw std::runtime_error("StaticProject::exportResultsToCSV(): reaction vector size mismatch.");
    }

    if (folderPath.empty()) {
        throw std::runtime_error("StaticProject::exportResultsToCSV(): folder path is empty.");
    }

    fs::path folder(folderPath);
    fs::create_directories(folder);
    fs::path filePath = folder / fileName;

    std::ofstream out(filePath);
    if (!out) {
        throw std::runtime_error("StaticProject::exportResultsToCSV(): cannot open file: " + filePath.string());
    }

    out << std::setprecision(16);
    out << "node_id,x,y,ux,uy,umag,rx,ry,rmag\n";

    for (std::size_t i = 0; i < mesh->nodes.size(); ++i) {
        const auto& node = mesh->nodes[i];
        const auto& dofs = mesh->node_dofs[i];

        if (!node || node->x.size() < 2 || dofs.size() < 2) {
            throw std::runtime_error("StaticProject::exportResultsToCSV(): invalid node data.");
        }

        const double x = node->x[0];
        const double y = node->x[1];

        const double ux = U(dofs[0]);
        const double uy = U(dofs[1]);
        const double umag = std::sqrt(ux * ux + uy * uy);

        const double rx = reactions(dofs[0]);
        const double ry = reactions(dofs[1]);
        const double rmag = std::sqrt(rx * rx + ry * ry);

        out << node->ID << ','
            << x << ','
            << y << ','
            << ux << ','
            << uy << ','
            << umag << ','
            << rx << ','
            << ry << ','
            << rmag << '\n';
    }
}

void StaticProject::exportNodesToCSV(const std::string& folderPath, const std::string& fileName) const {
    namespace fs = std::filesystem;

    if (!mesh) {
        throw std::runtime_error("StaticProject::exportNodesToCSV(): no mesh assigned.");
    }

    if (folderPath.empty()) {
        throw std::runtime_error("StaticProject::exportNodesToCSV(): folder path is empty.");
    }

    fs::path folder(folderPath);
    fs::create_directories(folder);
    fs::path filePath = folder / fileName;

    std::ofstream out(filePath);
    if (!out) {
        throw std::runtime_error("StaticProject::exportNodesToCSV(): cannot open file: " + filePath.string());
    }

    out << std::setprecision(16);
    out << "node_id,x,y\n";

    for (std::size_t i = 0; i < mesh->nodes.size(); ++i) {
        const auto& node = mesh->nodes[i];

        if (!node || node->x.size() < 2) {
            throw std::runtime_error("StaticProject::exportNodesToCSV(): invalid node coordinates.");
        }

        out << node->ID << ',' << node->x[0] << ',' << node->x[1] << '\n';
    }
}

void StaticProject::exportConnectivityToCSV(const std::string& folderPath, const std::string& fileName) const {
    namespace fs = std::filesystem;

    if (!mesh) {
        throw std::runtime_error("StaticProject::exportConnectivityToCSV(): no mesh assigned.");
    }

    if (folderPath.empty()) {
        throw std::runtime_error("StaticProject::exportConnectivityToCSV(): folder path is empty.");
    }

    fs::path folder(folderPath);
    fs::create_directories(folder);
    fs::path filePath = folder / fileName;

    std::ofstream out(filePath);
    if (!out) {
        throw std::runtime_error("StaticProject::exportConnectivityToCSV(): cannot open file: " + filePath.string());
    }

    out << "element_id,type,n1,n2,n3,n4\n";

    for (std::size_t e = 0; e < mesh->elements.size(); ++e) {
        const auto& elem = mesh->elements[e];
        const auto& conn = elem->connectivity();

        if (conn.size() == 3) {
            out << e << ",Tri3,"
                << conn[0] << ','
                << conn[1] << ','
                << conn[2] << ",-1\n";
        }
        else if (conn.size() == 4) {
            out << e << ",Quad4,"
                << conn[0] << ','
                << conn[1] << ','
                << conn[2] << ','
                << conn[3] << '\n';
        }
        else {
            throw std::runtime_error(
                "StaticProject::exportConnectivityToCSV(): unsupported element with " +
                std::to_string(conn.size()) + " nodes."
            );
        }
    }
}

void StaticProject::plotResults(const std::string& folderPath,
                                double scale,
                                const std::string& resultsFileName,
                                const std::string& connectivityFileName,
                                const std::string& imageFileName,
                                bool overwriteResults,
                                bool overwriteConnectivity) const {
    namespace fs = std::filesystem;

    const fs::path pythonExe =
        R"(C:\Users\User\Desktop\8o_e3amhno\Analisis_of_machenical_structurs\SAMS\venv\Scripts\python.exe)";
    const fs::path scriptPath =
        R"(C:\Users\User\Desktop\8o_e3amhno\Analisis_of_machenical_structurs\SAMS\plot_static.py)";

    if (folderPath.empty()) {
        throw std::runtime_error("StaticProject::plotResults(): folder path is empty.");
    }

    if (scale <= 0.0) {
        throw std::runtime_error("StaticProject::plotResults(): scale must be positive.");
    }

    if (!mesh) {
        throw std::runtime_error("StaticProject::plotResults(): no mesh assigned.");
    }

    if (U.size() == 0) {
        throw std::runtime_error("StaticProject::plotResults(): no static results available.");
    }

    if (!fs::exists(pythonExe)) {
        throw std::runtime_error("StaticProject::plotResults(): python executable not found: " + pythonExe.string());
    }

    if (!fs::exists(scriptPath)) {
        throw std::runtime_error("StaticProject::plotResults(): plot script not found: " + scriptPath.string());
    }

    fs::path folder(folderPath);
    fs::create_directories(folder);

    fs::path resultsPath = folder / resultsFileName;
    fs::path connPath    = folder / connectivityFileName;
    fs::path imagePath   = folder / imageFileName;

    if (overwriteResults || !fs::exists(resultsPath)) {
        exportResultsToCSV(folderPath, resultsFileName);
    }

    if (overwriteConnectivity || !fs::exists(connPath)) {
        exportConnectivityToCSV(folderPath, connectivityFileName);
    }

    std::ostringstream cmd;
    cmd << "cmd /c "
        << '\"'
        << '\"' << pythonExe.string() << '\"'
        << " "
        << '\"' << scriptPath.string() << '\"'
        << " "
        << '\"' << resultsPath.string() << '\"'
        << " "
        << '\"' << connPath.string() << '\"'
        << " --scale " << scale
        << " --output " << '\"' << imagePath.string() << '\"'
        << " --no-show"
        << '\"';

    std::cout << "Command: " << cmd.str() << '\n';

    const int rc = std::system(cmd.str().c_str());
    if (rc != 0) {
        throw std::runtime_error(
            "StaticProject::plotResults(): failed to run python script. Return code = " +
            std::to_string(rc)
        );
    }
}

void StaticProject::printSummary(std::ostream& os) const {
    os << "Project summary\n";
    os << "Name         : " << projectName << '\n';
    os << "Type         : StaticAnalysis\n";

    os << "State        : ";
    switch (state) {
        case ProjectState::Created:           os << "Created"; break;
        case ProjectState::MaterialAssigned:  os << "MaterialAssigned"; break;
        case ProjectState::GeometryImported:  os << "GeometryImported"; break;
        case ProjectState::MeshReady:         os << "MeshReady"; break;
        case ProjectState::MatricesAssembled: os << "MatricesAssembled"; break;
        case ProjectState::LoadsAssigned:     os << "LoadsAssigned"; break;
        case ProjectState::Ready:             os << "Ready"; break;
        case ProjectState::Completed:         os << "Completed"; break;
        default:                              os << "Unknown"; break;
    }
    os << '\n';

    if (mesh && mesh->F.size() > 0) {
        os << "||F||        : " << mesh->F.norm() << '\n';
        os << "F.size()     : " << mesh->F.size() << '\n';
    }

    if (U.size() > 0) {
        os << "||U||        : " << U.norm() << '\n';
        os << "U.size()     : " << U.size() << '\n';
    }

    if (reactions.size() > 0) {
        os << "||R||        : " << reactions.norm() << '\n';
        os << "R.size()     : " << reactions.size() << '\n';
    }
}

// Static Large Deformations Protject ==================================================================================
void StaticLDProject::setMonitorDOF(std::vector<double> x, int component){
    if (!mesh) {
        throw std::runtime_error("StaticLDProject::setMonitorDOF(): no mesh assigned.");
    }

    if (x.size() != 2) {
        throw std::runtime_error("StaticLDProject::setMonitorDOF(): coordinates must be 2D.");
    }

    if (component < 0 || component > 1) {
        throw std::runtime_error("StaticLDProject::setMonitorDOF(): component must be 0 (ux) or 1 (uy).");
    }

    const double tol = 1e-9;
    int nodeIndex = -1;

    for (std::size_t i = 0; i < mesh->nodes.size(); ++i) {
        const auto& node = mesh->nodes[i];
        if (!node || node->x.size() < 2) {
            continue;
        }

        if (std::abs(node->x[0] - x[0]) <= tol &&
            std::abs(node->x[1] - x[1]) <= tol) {
            nodeIndex = static_cast<int>(i);
            break;
            }
    }

    if (nodeIndex < 0) {
        throw std::runtime_error("StaticLDProject::setMonitorDOF(): no node found at given coordinates.");
    }

    const auto& dofs = mesh->node_dofs[nodeIndex];

    if (component >= static_cast<int>(dofs.size())) {
        throw std::runtime_error("StaticLDProject::setMonitorDOF(): component out of range for node DOFs.");
    }

    monitorDof = dofs[component];

    std::cout << "StaticLDProject::setMonitorDOF()\n";
    std::cout << "  node index : " << nodeIndex << '\n';
    std::cout << "  node ID    : " << mesh->nodes[nodeIndex]->ID << '\n';
    std::cout << "  coords     : (" << mesh->nodes[nodeIndex]->x[0]
              << ", " << mesh->nodes[nodeIndex]->x[1] << ")\n";
    std::cout << "  component  : " << component << '\n';
    std::cout << "  monitorDof : " << monitorDof << '\n';
}

void StaticLDProject::assemble() {
    if (!mesh) {
        throw std::runtime_error("StaticLDProject::assemble(): mesh not assigned.");
    }

    const int ndof = static_cast<int>(mesh->total_dofs);

    if (mesh->F.size() != ndof) {
        mesh->F = Eigen::VectorXd::Zero(ndof);
    }

    if (Fref.size() != ndof) {
        Fref = Eigen::VectorXd::Zero(ndof);
    }

    U = Eigen::VectorXd::Zero(ndof);
    reactions = Eigen::VectorXd::Zero(ndof);
    lambda = 0.0;

    setState(ProjectState::MatricesAssembled);
    setState(ProjectState::Ready);
}

void StaticLDProject::setLoad(const std::vector<double>& F_ext, const std::vector<double>& x) {
    if (!mesh) {
        throw std::runtime_error("StaticLDProject::setLoad(): no mesh assigned.");
    }

    if (x.size() != 2) {
        throw std::runtime_error("StaticLDProject::setLoad(): coordinates must be 2D.");
    }

    const double tol = 1e-9;
    int nodeIndex = -1;

    for (std::size_t i = 0; i < mesh->nodes.size(); ++i) {
        const auto& node = mesh->nodes[i];
        if (!node || node->x.size() < 2) continue;

        if (std::abs(node->x[0] - x[0]) <= tol && std::abs(node->x[1] - x[1]) <= tol) {
            nodeIndex = static_cast<int>(i);
            break;
            }
    }

    if (nodeIndex < 0) {
        throw std::runtime_error("StaticLDProject::setLoad(): no node found at given coordinates.");
    }

    const auto& dofs = mesh->node_dofs[nodeIndex];

    if (F_ext.size() != dofs.size()) {
        throw std::runtime_error("StaticLDProject::setLoad(): load size does not match node DOFs.");
    }

    if (Fref.size() == 0) {
        Fref = Eigen::VectorXd::Zero(mesh->total_dofs);
    }

    for (std::size_t i = 0; i < F_ext.size(); ++i) {
        Fref(dofs[i]) += F_ext[i];
    }

    setState(ProjectState::LoadsAssigned);
}

void StaticLDProject::solve() {
    switch (method) {
        case Solver::LoadControl:
            solveLoadControl();
        break;
        case Solver::DisplacmentControl:
            solveDisplacementControl();
        break;
        case Solver::ArcLength:
            solveArcLength();
        break;
        default:
            throw std::runtime_error("StaticLDProject::solve(): unknown solver method.");
    }
}

void StaticLDProject::solveLoadControl(){
    if (!mesh) {
        throw std::runtime_error("StaticLDProject::solve(): no mesh assigned.");
    }

    if (state != ProjectState::Ready && state != ProjectState::LoadsAssigned) {
        throw std::runtime_error("StaticLDProject::solve(): project not ready.");
    }

    const int ndof = static_cast<int>(mesh->total_dofs);

    if (Fref.size() != ndof) {
        throw std::runtime_error("StaticLDProject::solve(): reference load vector has wrong size.");
    }

    // Assign supports
    std::vector<int> fixed_dofs = mesh->getSupportDofs();
    std::vector<bool> is_fixed(ndof, false);

    for (int dof : fixed_dofs) {
        if (dof < 0 || dof >= ndof) {
            throw std::out_of_range("StaticLDProject::solve(): fixed DOF out of range.");
        }
        is_fixed[dof] = true;
    }

    std::vector<int> free_dofs;
    free_dofs.reserve(ndof);
    for (int i = 0; i < ndof; ++i) {
        if (!is_fixed[i]) free_dofs.push_back(i);
    }

    if (free_dofs.empty()) {
        throw std::runtime_error("StaticLDProject::solve(): all DOFs are constrained.");
    }

    // Initializing
    U = Eigen::VectorXd::Zero(ndof);
    lambda = 0.0;
    loadHistory.clear();

    const double dLambda = 1.0 / static_cast<double>(N_Steps); // Set step

    for (int step = 1; step <= N_Steps; ++step) {
        lambda += dLambda; // Load update

        Eigen::VectorXd Ustep = U;
        bool converged = false;
        double finalResidual = 0.0;
        double magnitudeInternalLoad = 0.0;
        double magnitudeAppliedExternalLoad = 0.0;
        int finalIterations = 0;

        // Newton - Raphson iterations
        for (int iter = 0; iter <= maxIterations; ++iter) { // Step. 1

            // Step. 2
            Eigen::VectorXd Fint = mesh->globalInternalForce(Ustep);

            // Step. 3
            Eigen::VectorXd Fext = lambda * Fref;
            Eigen::MatrixXd KT   = mesh->globalTangentStiffness(Ustep);

            // Checking
            if (Fint.size() != ndof) {
                throw std::runtime_error("StaticLDProject::solve(): internal force vector has wrong size.");
            }

            if (KT.rows() != ndof || KT.cols() != ndof) {
                throw std::runtime_error("StaticLDProject::solve(): tangent stiffness matrix has wrong dimensions.");
            }

            // Step. 4
            Eigen::VectorXd R = Fext - Fint; // Residual j_R( j_u_hat(i))

            Eigen::VectorXd R_f = reduceVector(R, free_dofs); // Reduced residual
            Eigen::MatrixXd KT_ff = reduceMatrix(KT, free_dofs); // Reduced Stiffness Matrix



            finalResidual = R_f.norm();
            magnitudeInternalLoad = Fint.norm();
            magnitudeAppliedExternalLoad = Fext.norm();
            finalIterations = iter;
            // Step. 6
            if (finalResidual < tolerance) {
                converged = true;
                break;
            }
            std::cout << "[Step " << step << "/" << N_Steps
                      << ", Iter " << iter
                      << "] lambda = " << lambda
                      << ", ||F_int|| = " << magnitudeInternalLoad
                      << ", lambda*||F_ext|| = " << magnitudeAppliedExternalLoad<< '\n';

            // Set-up solver
            Eigen::FullPivLU<Eigen::MatrixXd> lu(KT_ff);
            if (!lu.isInvertible()) {
                throw std::runtime_error("StaticLDProject::solve(): reduced tangent stiffness is singular.");
            }

            Eigen::LDLT<Eigen::MatrixXd> solver(KT_ff);
            if (solver.info() != Eigen::Success) {
                throw std::runtime_error("StaticLDProject::solve(): failed to factorize reduced tangent stiffness.");
            }

            // Solve for Δu_hat(i)
            Eigen::VectorXd dUf = solver.solve(R_f);
            if (solver.info() != Eigen::Success || !dUf.allFinite()) {
                throw std::runtime_error("StaticLDProject::solve(): failed to solve for displacement increment.");
            }

            // Expand of the resulted vector
            Eigen::VectorXd dU = expandVector(dUf, free_dofs, ndof);

            // Step. 5
            Ustep += dU;

            // Checking
            if (!Ustep.allFinite()) {
                throw std::runtime_error("StaticLDProject::solve(): non-finite displacement state.");
            }
        }

        if (!converged) {
            throw std::runtime_error(
                "StaticLDProject::solve(): Newton-Raphson did not converge at load step " +
                std::to_string(step)
            );
        }

        U = Ustep;

        // Save load-displacement history if monitor DOF is defined
        if (monitorDof >= 0) {
            if (monitorDof >= ndof) {
                throw std::runtime_error("StaticLDProject::solve(): monitorDof out of range.");
            }

            const double monitoredDisp = U(monitorDof);

            double monitoredLoad = 0.0;
            if (loadDof >= 0) {
                if (loadDof >= ndof) {
                    throw std::runtime_error("StaticLDProject::solve(): loadDof out of range.");
                }
                monitoredLoad = lambda * Fref(loadDof);
            } else {
                monitoredLoad = lambda * Fref.norm();
            }

            loadHistory.push_back({
                step,
                lambda,
                std::fabs(monitoredDisp),
                monitoredLoad,
                finalResidual,
                finalIterations
            });
        }
    }

    // Save Results
    mesh->U = U;
    Eigen::VectorXd Fext_final = lambda * Fref;
    Eigen::VectorXd Fint_final = mesh->globalInternalForce(U);
    reactions = Fint_final - Fext_final;

    for (int dof : free_dofs) {
        reactions(dof) = 0.0;
    }

    // Export load-displacement history if monitor DOF is defined
    if (monitorDof >= 0) {
        if (historyFolderPath.empty()) {
            throw std::runtime_error(
                "StaticLDProject::solve(): monitorDof is set, but no history output folder was defined."
            );
        }

        exportLoadHistoryToCSV(historyFolderPath, historyFileName);
    }

    setState(ProjectState::Completed);
}

void StaticLDProject::setControlDOF(const std::vector<double> &x, int component) {
    if (!mesh) {
        throw std::runtime_error("StaticLDProject::setMonitorDOF(): no mesh assigned.");
    }

    if (x.size() != 2) {
        throw std::runtime_error("StaticLDProject::setMonitorDOF(): coordinates must be 2D.");
    }

    if (component < 0 || component > 1) {
        throw std::runtime_error("StaticLDProject::setMonitorDOF(): component must be 0 (ux) or 1 (uy).");
    }

    const double tol = 1e-9;
    int nodeIndex = -1;

    for (std::size_t i = 0; i < mesh->nodes.size(); ++i) {
        const auto& node = mesh->nodes[i];
        if (!node || node->x.size() < 2) {
            continue;
        }

        if (std::abs(node->x[0] - x[0]) <= tol &&
            std::abs(node->x[1] - x[1]) <= tol) {
            nodeIndex = static_cast<int>(i);
            break;
            }
    }

    if (nodeIndex < 0) {
        throw std::runtime_error("StaticLDProject::setMonitorDOF(): no node found at given coordinates.");
    }

    const auto& dofs = mesh->node_dofs[nodeIndex];

    if (component >= static_cast<int>(dofs.size())) {
        throw std::runtime_error("StaticLDProject::setMonitorDOF(): component out of range for node DOFs.");
    }

    controlDof = dofs[component];
}
void StaticLDProject::solveDisplacementControl() {
    if (!mesh) {
        throw std::runtime_error("StaticLDProject::solveDisplacementControl(): no mesh assigned.");
    }

    if (state != ProjectState::Ready && state != ProjectState::LoadsAssigned) {
        throw std::runtime_error("StaticLDProject::solveDisplacementControl(): project not ready.");
    }

    const int ndof = static_cast<int>(mesh->total_dofs);

    if (Fref.size() != ndof) {
        throw std::runtime_error("StaticLDProject::solveDisplacementControl(): reference load vector has wrong size.");
    }

    if (controlDof < 0 || controlDof >= ndof) {
        throw std::runtime_error("StaticLDProject::solveDisplacementControl(): control DOF is not set or out of range.");
    }

    // Assign supports
    std::vector<int> fixed_dofs = mesh->getSupportDofs();
    std::vector<bool> is_fixed(ndof, false);

    for (int dof : fixed_dofs) {
        if (dof < 0 || dof >= ndof) {
            throw std::out_of_range("StaticLDProject::solveDisplacementControl(): fixed DOF out of range.");
        }
        is_fixed[dof] = true;
    }

    if (is_fixed[controlDof]) {
        throw std::runtime_error("StaticLDProject::solveDisplacementControl(): control DOF is constrained.");
    }

    std::vector<int> free_dofs;
    free_dofs.reserve(ndof);
    for (int i = 0; i < ndof; ++i) {
        if (!is_fixed[i]) free_dofs.push_back(i);
    }

    if (free_dofs.empty()) {
        throw std::runtime_error("StaticLDProject::solveDisplacementControl(): all DOFs are constrained.");
    }

    int controlIndexFree = -1;
    for (std::size_t i = 0; i < free_dofs.size(); ++i) {
        if (free_dofs[i] == controlDof) {
            controlIndexFree = static_cast<int>(i);
            break;
        }
    }

    if (controlIndexFree < 0) {
        throw std::runtime_error("StaticLDProject::solveDisplacementControl(): control DOF is not in free DOFs.");
    }

    // Initializing
    U = Eigen::VectorXd::Zero(ndof);
    lambda = 0.0;
    loadHistory.clear();

    const double dUbar = umax / static_cast<double>(N_Steps); // Displacement step
    const int nfree = static_cast<int>(free_dofs.size());
    double lambdaStep = lambda;
    for (int step = 1; step <= N_Steps; ++step) {
        const double uTarget = step * dUbar;

        Eigen::VectorXd Ustep = U;
        Ustep(controlDof) += dUbar; // predictor

        //lambdaStep = lambda;
        bool converged = false;
        double finalResidual = 0.0;
        int finalIterations = 0;

        // Newton - Raphson iterations
        for (int iter = 1; iter <= maxIterations; ++iter) { // Step. 1

            // Step. 2
            Eigen::VectorXd Fint = mesh->globalInternalForce(Ustep);
            Eigen::MatrixXd KT   = mesh->globalTangentStiffness(Ustep);

            // Checking
            if (Fint.size() != ndof) {
                throw std::runtime_error("StaticLDProject::solveDisplacementControl(): internal force vector has wrong size.");
            }

            if (KT.rows() != ndof || KT.cols() != ndof) {
                throw std::runtime_error("StaticLDProject::solveDisplacementControl(): tangent stiffness matrix has wrong dimensions.");
            }

            Eigen::VectorXd R = lambdaStep * Fref - Fint;

            Eigen::VectorXd R_f = reduceVector(R, free_dofs);
            Eigen::VectorXd Fref_f = reduceVector(Fref, free_dofs);
            Eigen::MatrixXd KT_ff = reduceMatrix(KT, free_dofs);

            Eigen::VectorXd U_f = reduceVector(Ustep, free_dofs);

            Eigen::VectorXd T_f = Eigen::VectorXd::Zero(nfree);
            T_f(controlIndexFree) = 1.0;

            // Step. 3
            Eigen::MatrixXd A = Eigen::MatrixXd::Zero(nfree + 1, nfree + 1);
            A.block(0, 0, nfree, nfree) = KT_ff;
            A.block(0, nfree, nfree, 1) = -Fref_f;
            A.block(nfree, 0, 1, nfree) = T_f.transpose();

            Eigen::VectorXd b = Eigen::VectorXd::Zero(nfree + 1);
            b.head(nfree) = R_f;
            b(nfree) = uTarget - T_f.dot(U_f);

            finalResidual = R_f.norm();
            finalIterations = iter;

            std::cout << "[Step " << step << "/" << N_Steps
                      << ", Iter " << iter
                      << "] lambda = " << lambdaStep
                      << ", ||R_f|| = " << finalResidual
                      << ", control = " << T_f.dot(U_f)
                      << ", target = " << uTarget << '\n';

            if (finalResidual < tolerance && std::abs(uTarget - T_f.dot(U_f)) < tolerance) {
                converged = true;
                break;
            }

            Eigen::FullPivLU<Eigen::MatrixXd> lu(A);

            if (!lu.isInvertible()) {
                throw std::runtime_error("StaticLDProject::solveDisplacementControl(): augmented reduced system is singular.");
            }

            Eigen::VectorXd sol = lu.solve(b);
            if (!sol.allFinite()) {
                throw std::runtime_error("StaticLDProject::solveDisplacementControl(): failed to solve augmented system.");
            }

            Eigen::VectorXd dUf = sol.head(nfree);
            double dLambda = sol(nfree);

            Eigen::VectorXd dU = expandVector(dUf, free_dofs, ndof);

            // Step. 4
            Ustep += dU;
            lambdaStep += dLambda;

            // Checking
            if (!Ustep.allFinite() || !std::isfinite(lambdaStep)) {
                throw std::runtime_error("StaticLDProject::solveDisplacementControl(): non-finite displacement state or load factor.");
            }
        }

        if (!converged) {
            throw std::runtime_error(
                "StaticLDProject::solveDisplacementControl(): Newton-Raphson did not converge at displacement step " +
                std::to_string(step)
            );
        }

        U = Ustep;
        lambda = lambdaStep;

        // Save load-displacement history if monitor DOF is defined
        if (monitorDof >= 0) {
            if (monitorDof >= ndof) {
                throw std::runtime_error("StaticLDProject::solveDisplacementControl(): monitorDof out of range.");
            }

            const double monitoredDisp = U(monitorDof);

            double monitoredLoad = 0.0;
            if (loadDof >= 0) {
                if (loadDof >= ndof) {
                    throw std::runtime_error("StaticLDProject::solveDisplacementControl(): loadDof out of range.");
                }
                monitoredLoad = lambda * Fref(loadDof);
            } else {
                monitoredLoad = lambda * Fref.norm();
            }

            loadHistory.push_back({
                step,
                lambda,
                std::fabs(monitoredDisp),
                monitoredLoad,
                finalResidual,
                finalIterations
            });
        }
    }

    // Save Results
    mesh->U = U;
    Eigen::VectorXd Fext_final = lambda * Fref;
    Eigen::VectorXd Fint_final = mesh->globalInternalForce(U);
    reactions = Fint_final - Fext_final;

    for (int dof : free_dofs) {
        reactions(dof) = 0.0;
    }

    // Export load-displacement history if monitor DOF is defined
    if (monitorDof >= 0) {
        if (historyFolderPath.empty()) {
            throw std::runtime_error(
                "StaticLDProject::solveDisplacementControl(): monitorDof is set, but no history output folder was defined."
            );
        }

        exportLoadHistoryToCSV(historyFolderPath, historyFileName);
    }

    setState(ProjectState::Completed);
}

void StaticLDProject::solveArcLength() {
    if (!mesh) {
        throw std::runtime_error("StaticLDProject::solveArcLength(): no mesh assigned.");
    }

    if (state != ProjectState::Ready && state != ProjectState::LoadsAssigned) {
        throw std::runtime_error("StaticLDProject::solveArcLength(): project not ready.");
    }

    const int ndof = static_cast<int>(mesh->total_dofs);

    if (Fref.size() != ndof) {
        throw std::runtime_error("StaticLDProject::solveArcLength(): reference load vector has wrong size.");
    }

    // Assign supports
    std::vector<int> fixed_dofs = mesh->getSupportDofs();
    std::vector<bool> is_fixed(ndof, false);

    for (int dof : fixed_dofs) {
        if (dof < 0 || dof >= ndof) {
            throw std::out_of_range("StaticLDProject::solveArcLength(): fixed DOF out of range.");
        }
        is_fixed[dof] = true;
    }

    std::vector<int> free_dofs;
    free_dofs.reserve(ndof);
    for (int i = 0; i < ndof; ++i) {
        if (!is_fixed[i]) free_dofs.push_back(i);
    }

    if (free_dofs.empty()) {
        throw std::runtime_error("StaticLDProject::solveArcLength(): all DOFs are constrained.");
    }

    const int nfree = static_cast<int>(free_dofs.size());

    // Arc-length parameters (Crisfield 1981)
    const double ds = arcLengthRadius;         // Δs (arc-length radius)
    const double psi = arcLengthAlpha;         // ψ (load scaling parameter)
    const double psi2 = psi * psi;

    if (ds <= 0.0) {
        throw std::runtime_error("StaticLDProject::solveArcLength(): arc-length radius must be positive.");
    }

    // Initializing
    U = Eigen::VectorXd::Zero(ndof);
    lambda = 0.0;
    loadHistory.clear();

    for (int step = 1; step <= N_Steps; ++step) {

        // Previous equilibrium point (u0, lambda0) — the center of the sphere
        Eigen::VectorXd U0 = U;
        double lambda0 = lambda;

        // Step. 0: Predictor (elastic predictor using tangent stiffness)
        Eigen::MatrixXd KT0 = mesh->globalTangentStiffness(U0);  // FIXED: MatrixXd, not VectorXd

        if (KT0.rows() != ndof || KT0.cols() != ndof) {
            throw std::runtime_error("StaticLDProject::solveArcLength(): tangent stiffness matrix has wrong dimensions in predictor.");
        }

        Eigen::MatrixXd KT0_ff = reduceMatrix(KT0, free_dofs);
        Eigen::VectorXd Fref_f = reduceVector(Fref, free_dofs);

        Eigen::FullPivLU<Eigen::MatrixXd> lu0(KT0_ff);
        if (!lu0.isInvertible()) {
            throw std::runtime_error("StaticLDProject::solveArcLength(): predictor tangent is singular.");
        }

        Eigen::VectorXd du_dir_f = lu0.solve(Fref_f);

        const double denom2 = du_dir_f.squaredNorm() + psi2;
        if (denom2 <= 0.0) {
            throw std::runtime_error("StaticLDProject::solveArcLength(): invalid predictor denominator.");
        }

        double dLambdaPred = ds / std::sqrt(denom2);
        Eigen::VectorXd dU_pred_f = dLambdaPred * du_dir_f;
        Eigen::VectorXd dU_pred = expandVector(dU_pred_f, free_dofs, ndof);

        Eigen::VectorXd Ustep = U0 + dU_pred;
        double lambdaStep = lambda0 + dLambdaPred;

        bool converged = false;
        double finalResidual = 0.0;
        int finalIterations = 0;

        // Newton - Raphson iterations
        for (int iter = 1; iter <= maxIterations; ++iter) { // Step. 1

            // Step. 2
            Eigen::VectorXd Fint = mesh->globalInternalForce(Ustep);
            Eigen::MatrixXd KT   = mesh->globalTangentStiffness(Ustep);

            // Checking
            if (Fint.size() != ndof) {
                throw std::runtime_error("StaticLDProject::solveArcLength(): internal force vector has wrong size.");
            }

            if (KT.rows() != ndof || KT.cols() != ndof) {
                throw std::runtime_error("StaticLDProject::solveArcLength(): tangent stiffness matrix has wrong dimensions.");
            }

            // Step. 3
            Eigen::VectorXd R = lambdaStep * Fref - Fint;

            Eigen::VectorXd R_f = reduceVector(R, free_dofs);
            Eigen::MatrixXd KT_ff = reduceMatrix(KT, free_dofs);

            Eigen::VectorXd Ustep_f = reduceVector(Ustep, free_dofs);
            Eigen::VectorXd U0_f    = reduceVector(U0, free_dofs);

            // Step. 4
            // Crisfield spherical arc-length constraint:
            // g = sqrt((u-u0)^T (u-u0) + psi^2 (lambda-lambda0)^2) - ds
            Eigen::VectorXd dUarc_f = Ustep_f - U0_f;
            double dLambdaArc = lambdaStep - lambda0;

            double normArc = std::sqrt(
                dUarc_f.squaredNorm() + psi2 * dLambdaArc * dLambdaArc
            );

            double g = normArc - ds;

            Eigen::VectorXd h_f = Eigen::VectorXd::Zero(nfree);
            double s = 0.0;

            if (normArc > 1e-14) {
                h_f = dUarc_f / normArc;    // h = ∂g/∂u = (u - u0) / normArc
                s   = psi2 * dLambdaArc / normArc; // s = ∂g/∂lambda = psi^2 (lambda - lambda0) / normArc
            } else {
                throw std::runtime_error("StaticLDProject::solveArcLength(): zero arc-length norm in constraint evaluation.");
            }

            finalResidual = R_f.norm();
            finalIterations = iter;

            std::cout << "[Step " << step << "/" << N_Steps
                      << ", Iter " << iter
                      << "] lambda = " << lambdaStep
                      << ", ||R_f|| = " << finalResidual
                      << ", g = " << g << '\n';

            // Step. 6
            // Check convergence of both equilibrium and arc-length constraint
            if (finalResidual < tolerance && std::abs(g) < tolerance) {
                converged = true;
                break;
            }

            // Step. 5
            // Solve the augmented arc-length system (Crisfield 1981):
            // [ KT   -Fref ] [δu     ] = [ R ]
            // [ h^T   s    ] [δlambda]   [ -g ]
            Eigen::MatrixXd A = Eigen::MatrixXd::Zero(nfree + 1, nfree + 1);
            A.block(0, 0, nfree, nfree)   = KT_ff;
            A.block(0, nfree, nfree, 1)   = -Fref_f;
            A.block(nfree, 0, 1, nfree)   = h_f.transpose();
            A(nfree, nfree)               = s;

            Eigen::VectorXd b = Eigen::VectorXd::Zero(nfree + 1);
            b.head(nfree) = R_f;
            b(nfree) = -g;

            Eigen::FullPivLU<Eigen::MatrixXd> lu(A);
            if (!lu.isInvertible()) {
                throw std::runtime_error("StaticLDProject::solveArcLength(): augmented system is singular.");
            }

            Eigen::VectorXd sol = lu.solve(b);
            if (!sol.allFinite()) {
                throw std::runtime_error("StaticLDProject::solveArcLength(): failed to solve augmented system.");
            }

            Eigen::VectorXd dUf = sol.head(nfree);
            double dLambda = sol(nfree);

            Eigen::VectorXd dU = expandVector(dUf, free_dofs, ndof);

            // Step. 5
            Ustep += dU;
            lambdaStep += dLambda;

            // Checking
            if (!Ustep.allFinite() || !std::isfinite(lambdaStep)) {
                throw std::runtime_error("StaticLDProject::solveArcLength(): non-finite displacement state or load factor.");
            }
        }

        if (!converged) {
            throw std::runtime_error(
                "StaticLDProject::solveArcLength(): Newton-Raphson did not converge at arc-length step " +
                std::to_string(step)
            );
        }

        U = Ustep;
        lambda = lambdaStep;

        // Save load-displacement history if monitor DOF is defined
        if (monitorDof >= 0) {
            if (monitorDof >= ndof) {
                throw std::runtime_error("StaticLDProject::solveArcLength(): monitorDof out of range.");
            }

            const double monitoredDisp = U(monitorDof);

            double monitoredLoad = 0.0;
            if (loadDof >= 0) {
                if (loadDof >= ndof) {
                    throw std::runtime_error("StaticLDProject::solveArcLength(): loadDof out of range.");
                }
                monitoredLoad = lambda * Fref(loadDof);
            } else {
                monitoredLoad = lambda * Fref.norm();
            }

            loadHistory.push_back({
                step,
                lambda,
                std::fabs(monitoredDisp),
                monitoredLoad,
                finalResidual,
                finalIterations
            });
        }
    }

    // Save Results
    mesh->U = U;
    Eigen::VectorXd Fext_final = lambda * Fref;
    Eigen::VectorXd Fint_final = mesh->globalInternalForce(U);
    reactions = Fint_final - Fext_final;

    for (int dof : free_dofs) {
        reactions(dof) = 0.0;
    }

    // Export load-displacement history if monitor DOF is defined
    if (monitorDof >= 0) {
        if (historyFolderPath.empty()) {
            throw std::runtime_error(
                "StaticLDProject::solveArcLength(): monitorDof is set, but no history output folder was defined."
            );
        }

        exportLoadHistoryToCSV(historyFolderPath, historyFileName);
    }

    setState(ProjectState::Completed);
}

void StaticLDProject::exportResultsToCSV(const std::string& folderPath, const std::string& fileName) const {
    namespace fs = std::filesystem;

    if (!mesh) {
        throw std::runtime_error("StaticLDProject::exportResultsToCSV(): no mesh assigned.");
    }

    if (U.size() == 0) {
        throw std::runtime_error("StaticLDProject::exportResultsToCSV(): no static LD results available.");
    }

    if (reactions.size() != U.size()) {
        throw std::runtime_error("StaticLDProject::exportResultsToCSV(): reaction vector size mismatch.");
    }

    if (folderPath.empty()) {
        throw std::runtime_error("StaticLDProject::exportResultsToCSV(): folder path is empty.");
    }

    fs::path folder(folderPath);
    fs::create_directories(folder);
    fs::path filePath = folder / fileName;

    std::ofstream out(filePath);
    if (!out) {
        throw std::runtime_error("StaticLDProject::exportResultsToCSV(): cannot open file: " + filePath.string());
    }

    out << std::setprecision(16);
    out << "node_id,x,y,ux,uy,umag,rx,ry,rmag\n";

    for (std::size_t i = 0; i < mesh->nodes.size(); ++i) {
        const auto& node = mesh->nodes[i];
        const auto& dofs = mesh->node_dofs[i];

        if (!node || node->x.size() < 2 || dofs.size() < 2) {
            throw std::runtime_error("StaticLDProject::exportResultsToCSV(): invalid node data.");
        }

        const double x = node->x[0];
        const double y = node->x[1];

        const double ux = U(dofs[0]);
        const double uy = U(dofs[1]);
        const double umag = std::sqrt(ux * ux + uy * uy);

        const double rx = reactions(dofs[0]);
        const double ry = reactions(dofs[1]);
        const double rmag = std::sqrt(rx * rx + ry * ry);

        out << node->ID << ','
            << x << ','
            << y << ','
            << ux << ','
            << uy << ','
            << umag << ','
            << rx << ','
            << ry << ','
            << rmag << '\n';
    }
}

void StaticLDProject::exportNodesToCSV(const std::string& folderPath, const std::string& fileName) const {
    namespace fs = std::filesystem;

    if (!mesh) {
        throw std::runtime_error("StaticProject::exportNodesToCSV(): no mesh assigned.");
    }

    if (folderPath.empty()) {
        throw std::runtime_error("StaticProject::exportNodesToCSV(): folder path is empty.");
    }

    fs::path folder(folderPath);
    fs::create_directories(folder);
    fs::path filePath = folder / fileName;

    std::ofstream out(filePath);
    if (!out) {
        throw std::runtime_error("StaticProject::exportNodesToCSV(): cannot open file: " + filePath.string());
    }

    out << std::setprecision(16);
    out << "node_id,x,y\n";

    for (std::size_t i = 0; i < mesh->nodes.size(); ++i) {
        const auto& node = mesh->nodes[i];

        if (!node || node->x.size() < 2) {
            throw std::runtime_error("StaticProject::exportNodesToCSV(): invalid node coordinates.");
        }

        out << node->ID << ',' << node->x[0] << ',' << node->x[1] << '\n';
    }
}

void StaticLDProject::exportConnectivityToCSV(const std::string& folderPath, const std::string& fileName) const {
    namespace fs = std::filesystem;

    if (!mesh) {
        throw std::runtime_error("StaticProject::exportConnectivityToCSV(): no mesh assigned.");
    }

    if (folderPath.empty()) {
        throw std::runtime_error("StaticProject::exportConnectivityToCSV(): folder path is empty.");
    }

    fs::path folder(folderPath);
    fs::create_directories(folder);
    fs::path filePath = folder / fileName;

    std::ofstream out(filePath);
    if (!out) {
        throw std::runtime_error("StaticProject::exportConnectivityToCSV(): cannot open file: " + filePath.string());
    }

    out << "element_id,type,n1,n2,n3,n4\n";

    for (std::size_t e = 0; e < mesh->elements.size(); ++e) {
        const auto& elem = mesh->elements[e];
        const auto& conn = elem->connectivity();

        if (conn.size() == 3) {
            out << e << ",Tri3,"
                << conn[0] << ','
                << conn[1] << ','
                << conn[2] << ",-1\n";
        }
        else if (conn.size() == 4) {
            out << e << ",Quad4,"
                << conn[0] << ','
                << conn[1] << ','
                << conn[2] << ','
                << conn[3] << '\n';
        }
        else {
            throw std::runtime_error(
                "StaticProject::exportConnectivityToCSV(): unsupported element with " +
                std::to_string(conn.size()) + " nodes."
            );
        }
    }
}

void StaticLDProject::exportLoadHistoryToCSV(const std::string& folderPath, const std::string& fileName) const {
    namespace fs = std::filesystem;

    if (folderPath.empty()) {
        throw std::runtime_error("StaticLDProject::exportLoadHistoryToCSV(): folder path is empty.");
    }

    fs::path folder(folderPath);
    fs::create_directories(folder);
    fs::path filePath = folder / fileName;

    std::ofstream out(filePath);
    if (!out) {
        throw std::runtime_error("StaticLDProject::exportLoadHistoryToCSV(): cannot open file: " + filePath.string());
    }

    out << std::setprecision(16);
    out << "step,lambda,displacement,load,residual,iterations\n";

    for (const auto& rec : loadHistory) {
        out << rec.step << ','
            << rec.lambda << ','
            << rec.displacement << ','
            << rec.load << ','
            << rec.residual << ','
            << rec.iterations << '\n';
    }
}

void StaticLDProject::plotResults(const std::string& folderPath,
                                  double scale,
                                  const std::string& resultsFileName,
                                  const std::string& connectivityFileName,
                                  const std::string& imageFileName,
                                  bool overwriteResults,
                                  bool overwriteConnectivity) const {
    namespace fs = std::filesystem;

    const fs::path pythonExe =
        R"(C:\Users\User\Desktop\8o_e3amhno\Analisis_of_machenical_structurs\SAMS\venv\Scripts\python.exe)";
    const fs::path scriptPath =
        R"(C:\Users\User\Desktop\8o_e3amhno\Analisis_of_machenical_structurs\SAMS\plot_staticLD.py)";

    if (folderPath.empty()) {
        throw std::runtime_error("StaticLDProject::plotResults(): folder path is empty.");
    }

    if (scale <= 0.0) {
        throw std::runtime_error("StaticLDProject::plotResults(): scale must be positive.");
    }

    if (!mesh) {
        throw std::runtime_error("StaticLDProject::plotResults(): no mesh assigned.");
    }

    if (U.size() == 0) {
        throw std::runtime_error("StaticLDProject::plotResults(): no static LD results available.");
    }

    if (!fs::exists(pythonExe)) {
        throw std::runtime_error(
            "StaticLDProject::plotResults(): python executable not found: " + pythonExe.string()
        );
    }

    if (!fs::exists(scriptPath)) {
        throw std::runtime_error(
            "StaticLDProject::plotResults(): plot script not found: " + scriptPath.string()
        );
    }

    fs::path folder(folderPath);
    fs::create_directories(folder);

    fs::path resultsPath = folder / resultsFileName;
    fs::path connPath    = folder / connectivityFileName;
    fs::path imagePath   = folder / imageFileName;
    fs::path historyPath = folder / historyFileName;

    if (overwriteResults || !fs::exists(resultsPath)) {
        exportResultsToCSV(folderPath, resultsFileName);
    }

    if (overwriteConnectivity || !fs::exists(connPath)) {
        exportConnectivityToCSV(folderPath, connectivityFileName);
    }

    const bool hasHistory = (monitorDof >= 0) && fs::exists(historyPath);

    std::ostringstream cmd;
    cmd << "cmd /c "
        << '\"'
        << '\"' << pythonExe.string() << '\"'
        << " "
        << '\"' << scriptPath.string() << '\"'
        << " "
        << '\"' << resultsPath.string() << '\"'
        << " "
        << '\"' << connPath.string() << '\"';

    if (hasHistory) {
        cmd << " --history " << '\"' << historyPath.string() << '\"';
    }

    cmd << " --scale " << scale
        << " --output " << '\"' << imagePath.string() << '\"'
        << " --no-show"
        << '\"';

    std::cout << "Command: " << cmd.str() << '\n';

    const int rc = std::system(cmd.str().c_str());
    if (rc != 0) {
        throw std::runtime_error(
            "StaticLDProject::plotResults(): failed to run python script. Return code = " +
            std::to_string(rc)
        );
    }
}

void StaticLDProject::printSummary(std::ostream& os) const {
    os << "Project summary\n";
    os << "Name         : " << projectName << '\n';
    os << "Type         : StaticLDAnalysis\n";

    os << "State        : ";
    switch (state) {
        case ProjectState::Created:           os << "Created"; break;
        case ProjectState::MaterialAssigned:  os << "MaterialAssigned"; break;
        case ProjectState::GeometryImported:  os << "GeometryImported"; break;
        case ProjectState::MeshReady:         os << "MeshReady"; break;
        case ProjectState::MatricesAssembled: os << "MatricesAssembled"; break;
        case ProjectState::LoadsAssigned:     os << "LoadsAssigned"; break;
        case ProjectState::Ready:             os << "Ready"; break;
        case ProjectState::Completed:         os << "Completed"; break;
        default:                              os << "Unknown"; break;
    }
    os << '\n';

    os << "N_Steps       : " << N_Steps << '\n';
    os << "lambda       : " << lambda << '\n';
    os << "maxIterations: " << maxIterations << '\n';
    os << "tolerance    : " << tolerance << '\n';

    if (Fref.size() > 0) {
        os << "||Fref||     : " << Fref.norm() << '\n';
        os << "Fref.size()  : " << Fref.size() << '\n';
    }

    if (U.size() > 0) {
        os << "||U||        : " << U.norm() << '\n';
        os << "U.size()     : " << U.size() << '\n';
    }

    if (reactions.size() > 0) {
        os << "||R||        : " << reactions.norm() << '\n';
        os << "R.size()     : " << reactions.size() << '\n';
    }
}