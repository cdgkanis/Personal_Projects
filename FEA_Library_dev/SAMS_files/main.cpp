//Project libraries
#include "Commands.h"
#include "Import_Handler.h"
#include "Physics.h"
#include <filesystem>
#include <iostream>
#include <vector>
#include <cmath>
#include <memory>
using namespace std;
namespace fs = std::filesystem;

// Definitions
#define Geometry_Path "C:\\Users\\User\\Desktop\\8o_e3amhno\\Analisis_of_machenical_structurs\\SAMS\\Geometry.msh"
#define TRresults R"(C:\Users\User\Desktop\TRresults)"
#define EFresults R"(C:\Users\User\Desktop\EFresults)"

#define Geometry_Path2 "C:\\Users\\User\\Desktop\\8o_e3amhno\\Analisis_of_machenical_structurs\\SAMS\\Geometry2.msh"
#define SLresults R"(C:\Users\User\Desktop\SLresults)"
#define SNLresults R"(C:\Users\User\Desktop\SNLresults)"

std::function<Eigen::VectorXd(double)> F = [](double t) {
    Eigen::VectorXd F_(2);
    F_.setZero();

    if (t >= 0.0 && t <= 0.5) {
        F_(0) = 0.0;
        F_(1) = -100.0 * std::sin(20.0 * t) * 1e3;
    }

    return F_;
};

int main() {
    try {

        // Erotima 1
        auto mater = loadMaterial("Steel");

        const std::vector<std::string> modal_assumptions = {
            "PlaneStress",
            "Consistent"
        };
        /*
        //------ Beam's Eigen Frequencies Study
        EigenFrequencyProject project("Plate eigen-frequency test");
        project.setNumModes(10);
        project.assignMaterial(mater);
        project.importGeometry(Geometry_Path, mater, 0.2, 2, 4);

        project.setupElementAssumptions(modal_assumptions);

        // Supports
        project.mesh->addSupportOnEdge('L', SupportType::Fixed);

        project.assemble();
        project.solve();
        project.printSummary(std::cout);
        project.exportModesToCSV(EFresults, 4);
        project.exportConnectivityToCSV(EFresults);
        EigenFrequencyProject::plotModes(EFresults, 4);

        std::vector<double> res = project.getEigenFrequencies();

        //------ Beam's Transient Study
        if (res.size() >= 2) {
            project.mesh->globalDumpingMatrix_RayleighMethode(res.at(0), res.at(1), 0.05);
        }


        TransientProject TRproject = TransientProject("Transient study");
        TRproject.mesh = std::move(project.mesh);

        TRproject.setTimeSpan(0.0, 1.0, 0.1e-4);
        std::vector<double> F_place = {4, 0};

        TRproject.setVariableLoad(F, F_place);
        TRproject.assemble();
        TRproject.setSolver(Solver::Neumark);
        TRproject.solve();

        TRproject.plotResults(TRresults, 40, "transient_results.csv", true);

        TRproject.animateSolution(
            TRresults,
            20.0,
            20,
            "transient_results.csv",
            "nodes.csv",
            "connectivity.csv",
            true,
            true,
            true
        );

        TRproject.printSummary(std::cout);

        */
        // Erotima 2
        /*
        //------ Beam's Static Study
        StaticProject SLProject("Static Linear Analysis");
        SLProject.assignMaterial(mater);
        SLProject.importGeometry(Geometry_Path2, mater, 0.2, 2, 5);

        SLProject.setupElementAssumptions(modal_assumptions);
        SLProject.mesh->addSupportOnEdge('L', SupportType::Fixed);

        SLProject.setLoad({0.0, -20e6}, {6.0, 0.0});
        SLProject.assemble();

        SLProject.solve();
        SLProject.printSummary(std::cout);

        SLProject.plotResults(
            SLresults,
            0.20,
            "static_results.csv",
            "connectivity.csv",
            "static_plot.png",
            true,
            true
        );
        */
        //------ Beam's Static Large Deformation Study
        StaticLDProject SLDProject("Static Linear Analysis");
        SLDProject.assignMaterial(mater);
        SLDProject.importGeometry(Geometry_Path2, mater, 0.2, 2, 6);

        SLDProject.setupElementAssumptions(modal_assumptions);
        SLDProject.mesh->addSupportOnEdge('L', SupportType::Fixed);

        SLDProject.setLoad({0.0, -20e6}, {6.0, 0.0});
        SLDProject.assemble();

        SLDProject.setControlDOF({6.0, 0.0}, 1);
        SLDProject.setMonitorDOF({6.0, 0.0}, 1);
        SLDProject.setMaximumDisplacement(-2.3);
        SLDProject.setSolver(Solver::LoadControl);
        SLDProject.setHistoryOutput(SNLresults, "load_history.csv");
        SLDProject.solve();
        SLDProject.printSummary(std::cout);

        SLDProject.plotResults(
            SNLresults,
            1.0,
            "static_results.csv",
            "connectivity.csv",
            "static_plotLD.png",
            true,
            true
        );

    }catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }


    return 0;
}
