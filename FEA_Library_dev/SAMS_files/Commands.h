//
// Created by User on 1/5/2026.
//

#ifndef COMMANDS_H
#define COMMANDS_H

#pragma once

#include <cmath>
#include <string>
#include <vector>
#include <memory>
#include <fstream>
#include <stdexcept>
#include <iomanip>
#include <ostream>
#include <functional>
#include "ExternalLibraries/json.hpp"
#include <filesystem>
#include <sstream>
#include <cstdlib>
#include <iostream>
using json = nlohmann::json;

#include "Physics.h"
#include "Import_Handler.h"

// Globals =============================================================================================================
inline constexpr const char* Mat_Library = "Materials.json";

// Material library ====================================================================================================
std::shared_ptr<Material> loadMaterial(const std::string& materialName,const std::string& filePath = Mat_Library);

void addMaterialToLibrary(const Material& material,const std::string& filePath = Mat_Library,bool overwrite = false);

std::shared_ptr<Material> createMaterialFromUserInput();

// Project Library =====================================================================================================
enum class ProjectType {
    EigenFrequencies,
    StaticAnalysis,
    TransientAnalysis,
    StaticLDAnalysis
};

enum class ProjectState {
    Created,
    MaterialAssigned,
    GeometryImported,
    MeshReady,
    MatricesAssembled,
    LoadsAssigned,
    Ready,
    Completed
};

enum class Solver{
    Neumark,
    CentralDifs,
    LoadControl,
    DisplacmentControl,
    ArcLength
};

class Project {
    protected:
        std::string projectName;
        ProjectType projectType;
        ProjectState state = ProjectState::Created;

        std::shared_ptr<Material> defaultMaterial;
        std::string geometryFile;


        void setState(ProjectState newState) { state = newState; }

    public:
        std::unique_ptr<Mesh> mesh;

        Project(std::string name, ProjectType type)
            : projectName(std::move(name)), projectType(type) {}

        virtual ~Project() = default;

        [[nodiscard]] const std::string& name() const { return projectName; }
        [[nodiscard]] ProjectType type() const { return projectType; }
        [[nodiscard]] ProjectState getState() const { return state; }
        [[nodiscard]] bool hasMesh() const { return static_cast<bool>(mesh); }

        [[nodiscard]] bool isReady() const {
            return state == ProjectState::Ready;
        }

        [[nodiscard]] bool isCompleted() const {
            return state == ProjectState::Completed;
        }

        virtual void assignMaterial(const std::shared_ptr<Material>& material);
        virtual void importGeometry(const std::string& mshFile,const std::shared_ptr<Material>& material,double thickness,unsigned int dofsPerNode,double meshDensity);

        void setupElementAssumptions(const std::vector<std::string>& assumptions) const;
        void removeElementAssumptions(AssumptionType assumptions) const;

        void addNodalLoad(unsigned int node_id, const std::vector<double>& values) const;
        void setNodalLoad(unsigned int node_id, const std::vector<double>& values) const;

        virtual void assemble() = 0;
        virtual void solve() = 0;
        virtual void printSummary(std::ostream& os) const;
};

class EigenFrequencyProject : public Project {
    protected:
        Eigen::VectorXd omega;
        Eigen::VectorXd frequenciesHz;
        Eigen::MatrixXd modes;
        int numModes = 12;

    public:
        explicit EigenFrequencyProject(const std::string& name): Project(name, ProjectType::EigenFrequencies) {}

        void setNumModes(int n);
        void assemble() override;
        void solve() override;
        [[nodiscard]] std::vector<double> getEigenFrequencies() const;

        void exportModesToCSV(const std::string& folderPath, int nModes) const;
        void exportModeToCSV(const std::string& filePath, int modeIndex) const;
        void exportConnectivityToCSV(const std::string& filePath) const;
        static void plotModes(const std::string& folderPath, int nModes);
        void printSummary(std::ostream& os) const override;
};

class TransientProject : public Project {
    protected:
        double t0 = 0.0;
        double tf = 1.0;
        double dt = 0.01;

        double beta  = 0.25;
        double gamma = 0.5;

        Eigen::VectorXd U0;
        Eigen::VectorXd V0;
        Eigen::VectorXd A0;

        std::vector<double> time;
        std::vector<Eigen::VectorXd> U_hist;
        std::vector<Eigen::VectorXd> V_hist;
        std::vector<Eigen::VectorXd> A_hist;

        std::vector<std::pair<std::function<Eigen::VectorXd(double)>, std::vector<double>>> loadFunctions;
        Solver method = Solver::Neumark;

    public:
        explicit TransientProject(const std::string& name): Project(name, ProjectType::TransientAnalysis) {}

        void setTimeSpan(double t_start, double t_final, double time_step);
        void setInitialConditions(const Eigen::VectorXd& u0, const Eigen::VectorXd& v0);

        void assemble() override;
        void setSolver(Solver solver) { method = solver; }

        void solve() override;
        void printSummary(std::ostream& os) const override;

        void setVariableLoad(const std::function<Eigen::VectorXd(double)>& f, const std::vector<double>& x);

        void exportResultsToCSV(const std::string& folderPath, const std::string& fileName = "transient_results.csv") const;
        void exportNodesToCSV(const std::string& folderPath, const std::string& fileName = "nodes.csv") const;
        void exportConnectivityToCSV(const std::string& folderPath, const std::string& fileName = "connectivity.csv") const;
        void plotResults(const std::string& folderPath, int dof, const std::string& fileName = "transient_results.csv", bool overwrite = false) const;
        void animateSolution(const std::string& folderPath, double scale = 1.0, int fps = 30, const std::string& resultsFileName = "transient_results.csv", const std::string& nodesFileName = "nodes.csv", const std::string& connectivityFileName = "connectivity.csv", bool overwriteResults = false, bool overwriteGeometry = false, bool saveGif = false) const;
};

class StaticProject : public Project {
    protected:
        Eigen::VectorXd U;
        Eigen::VectorXd reactions;

    public:
        explicit StaticProject(const std::string& name)
            : Project(name, ProjectType::StaticAnalysis) {}

        void assemble() override;
        void solve() override;
        void printSummary(std::ostream& os) const override;

        void setLoad(const std::vector<double>& F_ext, const std::vector<double>& x);

        void exportResultsToCSV(const std::string& folderPath,const std::string& fileName = "static_results.csv") const;
        void exportNodesToCSV(const std::string& folderPath, const std::string& fileName = "nodes.csv") const;
        void exportConnectivityToCSV(const std::string& folderPath, const std::string& fileName = "connectivity.csv") const;
        void plotResults(const std::string& folderPath, double scale = 10, const std::string& resultsFileName = "static_results.csv",
                 const std::string& connectivityFileName = "connectivity.csv", const std::string& imageFileName = "static_plot.png",
                 bool overwriteResults = false,bool overwriteConnectivity = false) const;

        [[nodiscard]] const Eigen::VectorXd& displacements() const { return U; }
};


class StaticLDProject : public Project {
    protected:
        Eigen::VectorXd U;
        Eigen::VectorXd reactions;

        Eigen::VectorXd Fref;
        int N_Steps = 100;
        int maxIterations = 100;
        double tolerance = 1e-4;
        double lambda = 0.0;

        // Load - Deformation Diagram
        struct LoadStepRecord {
            int step;
            double lambda;
            double displacement;
            double load;
            double residual;
            int iterations;
        };

        std::vector<LoadStepRecord> loadHistory;
        int monitorDof = -1;
        int loadDof = -1;
        std::string historyFolderPath;
        std::string historyFileName = "load_history.csv";

        Solver method = Solver::LoadControl;

        // Displacement Control
        int controlDof = -1;
        double umax = 0.0;   // maximum prescribed displacement

        // Arc-length
        double arcLengthRadius = 1e-2;   // Δs
        double arcLengthPsi = 0.0;       // 1.0 = spherical, 0.0 = cylindrical
        double arcLengthAlpha = 0.001;
    public:
        explicit StaticLDProject(const std::string& name): Project(name, ProjectType::StaticLDAnalysis) {}

        void assemble() override;
        void setLoad(const std::vector<double>& F_ext, const std::vector<double>& x);

        void setSolver(Solver solver) { method = solver; }
        void solveLoadControl();
        void solveDisplacementControl();
        void solveArcLength();
        void solve() override;

        // Load Control
        void setNumberOfSteps(int steps) {
            if (steps <= 0) throw std::runtime_error("StaticLDProject::setNumberOfSteps(): steps must be positive.");
            N_Steps = steps;
        }
        void setTolerance(double tol) {
            if (tol <= 0.0) throw std::runtime_error("StaticLDProject::setTolerance(): tolerance must be positive.");
            tolerance = tol;
        }
        void setMonitorDOF(std::vector<double> x, int component);// component = 0 -> x, component = 1 -> y
        void setMaxIterations(int maxIts) {
            if (maxIts <= 0) throw std::runtime_error("StaticLDProject::setMaxIterations(): maxIterations must be positive.");
            maxIterations = maxIts;
        }

        // Displacement Control
        void setControlDOF(const std::vector<double>& x, int component);// component = 0 -> x, component = 1 -> y
        void setMaximumDisplacement(double uMax) {
            if (uMax == 0.0) {
                throw std::runtime_error("StaticLDProject::setMaximumDisplacement(): maximum displacement cannot be zero.");
            }
            umax = uMax;
        }

        // Arc Length
        void setArcLengthRadius(double ds) {
                if (ds <= 0.0) {
                    throw std::runtime_error("StaticLDProject::setArcLengthRadius(): radius must be positive.");
                }
                arcLengthRadius = ds;
            }
        void setArcLengthPsi(double psi) {
                if (psi < 0.0 || psi > 1.0) {
                    throw std::runtime_error("StaticLDProject::setArcLengthPsi(): psi must be in [0,1].");
                }
                arcLengthPsi = psi;
            }

        // Results
        [[nodiscard]] const Eigen::VectorXd& displacements() const { return U; }

        void exportResultsToCSV(const std::string& folderPath,const std::string& fileName = "static_results.csv") const;
        void exportNodesToCSV(const std::string& folderPath, const std::string& fileName = "nodes.csv") const;
        void exportConnectivityToCSV(const std::string& folderPath, const std::string& fileName = "connectivity.csv") const;
        void exportLoadHistoryToCSV(const std::string& folderPath, const std::string& fileName) const;

        void plotResults(const std::string& folderPath,
                 double scale = 10,
                 const std::string& resultsFileName = "static_ld_results.csv",
                 const std::string& connectivityFileName = "connectivity.csv",
                 const std::string& imageFileName = "static_ld_plot.png",
                 bool overwriteResults = false,
                 bool overwriteConnectivity = false) const;

        void printSummary(std::ostream& os) const override;

        void setHistoryOutput(const std::string& folderPath, const std::string& fileName) {
            if (folderPath.empty()) {
                throw std::runtime_error("StaticLDProject::setHistoryOutput(): folder path is empty.");
            }

            historyFolderPath = folderPath;
            historyFileName = fileName.empty() ? "load_history.csv" : fileName;
        }
};
#endif // COMMANDS_H