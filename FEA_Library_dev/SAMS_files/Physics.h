//
// Created by User on 30/4/2026.
//

#ifndef Physics_H
#define Physics_H

#pragma once
//Standard libraries
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <ExternalLibraries/Eigen/Dense>



//======================================================================================================================
// Definition and assigning Material
enum class MaterialType {
    Isotropic,
    Orthotropic
};


// General definition of material
class Material {
    protected:
        MaterialType Type;
        std::string Name;
        double rho_;

    public:
        Material(MaterialType type, std::string name, double rho): Type(type), Name(std::move(name)), rho_(rho) {}
        virtual ~Material() = default;

        [[nodiscard]] MaterialType type() const{return Type;}
        [[nodiscard]] const std::string& name() const{return Name;}
        [[nodiscard]] double rho() const{return rho_;}

        //Polymorfism
        virtual std::ostream& print(std::ostream& os) const = 0;
};

// Isotropic Material
class IsotropicMaterial : public Material {
    protected:
        double E_, nu_;

    public:
        IsotropicMaterial(std::string name, double rho, double E, double nu): Material(MaterialType::Isotropic, std::move(name), rho), E_(E), nu_(nu) {}

        // Basic accessing
        [[nodiscard]] double E() const { return E_; }
        [[nodiscard]] double nu() const { return nu_; }

        // Calculations
        [[nodiscard]] double G() const{return this->E_/(2*(1 + this->nu_));};

        std::ostream& print(std::ostream& os) const override;
};

// Orthotropic Material
class OrthotropicMaterial : public Material {
    protected:
        double Ex_, Ey_, Ez_;
        double nuxy_, nuyz_, nuxz_;
        double Gxy_, Gyz_, Gxz_;
    public:
        OrthotropicMaterial(std::string name, double rho, double Ex, double Ey, double Ez,
                                            double nuxy, double nuyz, double nuxz,
                                             double Gxy, double Gyz, double Gxz)
        : Material(MaterialType::Orthotropic, std::move(name), rho), Ex_(Ex), Ey_(Ey), Ez_(Ez),
                                                                    nuxy_(nuxy), nuyz_(nuyz), nuxz_(nuxz),
                                                                    Gxy_(Gxy), Gyz_(Gyz), Gxz_(Gxz){}
        // Basic accessing
        [[nodiscard]] double Ex()   const { return Ex_; }
        [[nodiscard]] double Ey()   const { return Ey_; }
        [[nodiscard]] double Ez()   const { return Ez_; }
        [[nodiscard]] double nuxy() const { return nuxy_; }
        [[nodiscard]] double nuyz() const { return nuyz_; }
        [[nodiscard]] double nuxz() const { return nuxz_; }
        [[nodiscard]] double Gxy()  const { return Gxy_; }
        [[nodiscard]] double Gyz()  const { return Gyz_; }
        [[nodiscard]] double Gxz()  const { return Gxz_; }

        std::ostream& print(std::ostream& os) const override;
};

inline std::ostream& operator<<(std::ostream& os, const Material& mat){return mat.print(os);}

//======================================================================================================================
// Definition of n dimensional Node
typedef struct Node_nDim{
    unsigned int ID;
    std::vector<double> x;
}Node;

// General definition of an element
enum class AssumptionType {
    PlaneStress,
    PlaneStrain,
    Lumped,
    Consistent,
    LargeDeformation
};



class Element {
    protected:
        unsigned int ID;
        std::vector<unsigned int> conn;
        std::vector<std::shared_ptr<Node>> nodes;

        std::shared_ptr<Material> material;
        std::vector<AssumptionType> assumptions;

    public:
        Element(unsigned int El_ID,
                std::vector<unsigned int> connectivity,
                std::vector<std::shared_ptr<Node>> Nodes,
                std::shared_ptr<Material> mat)
            : ID(El_ID),
              conn(std::move(connectivity)),
              nodes(std::move(Nodes)),
              material(std::move(mat)) {}

        virtual ~Element() = default;

        [[nodiscard]] const std::vector<unsigned int>& connectivity() const { return conn; }
        [[nodiscard]] const std::vector<AssumptionType>& getAssumptions() const { return assumptions; }

        void AddAssumption(AssumptionType assumption);
        void AddAssumption(const std::string& assumption);
        void RemoveAssumption(AssumptionType assumption);
        [[nodiscard]] bool HasAssumption(AssumptionType assumption) const;

        [[nodiscard]] virtual int numNodes() const = 0;
        [[nodiscard]] virtual int dimension() const = 0;
        [[nodiscard]] virtual std::string name() const = 0;

        [[nodiscard]] virtual double ShapeFunction(int i, double r, double s) const = 0;
        [[nodiscard]] virtual double dShapeFunction_dr(int i, double r, double s) const = 0;
        [[nodiscard]] virtual double dShapeFunction_ds(int i, double r, double s) const = 0;
        [[nodiscard]] virtual Eigen::MatrixXd J(double r, double s) const = 0;

        [[nodiscard]] virtual Eigen::MatrixXd stiffnessMatrix() const = 0;
        [[nodiscard]] virtual Eigen::MatrixXd massMatrix() const = 0;
        virtual std::ostream& print(std::ostream& os) const = 0;

        // Large Deformation
        [[nodiscard]] virtual Eigen::MatrixXd tangentStiffness(const Eigen::VectorXd& ue) const = 0;
        [[nodiscard]] virtual Eigen::VectorXd internalForce(const Eigen::VectorXd& ue) const = 0;
};

class Tri3 : public Element {
    public:
        double t;

    Tri3(unsigned int id, std::vector<unsigned int> connectivity, std::vector<std::shared_ptr<Node>> Nodes, std::shared_ptr<Material> mat, double Thickness): Element(id, std::move(connectivity), std::move(Nodes) ,std::move(mat)){
        if (conn.size() != 3) {throw std::runtime_error("Tri3(): Node number should be 3.");}

        if (Thickness > 0.0) {
            t = Thickness;
        } else {
            throw std::runtime_error("Tri3(): Thickness should be positive.");
        }
    }

        [[nodiscard]] int numNodes() const override { return 3; }
        [[nodiscard]] int dimension() const override { return 2; }
        [[nodiscard]] std::string name() const override { return "Tri3"; }

        [[nodiscard]] double ShapeFunction(int i, double r, double s) const override;
        [[nodiscard]] double dShapeFunction_dr(int i, double r, double s) const override;
        [[nodiscard]] double dShapeFunction_ds(int i, double r, double s) const override;
        [[nodiscard]] Eigen::MatrixXd J(double r, double s) const override;
        [[nodiscard]] Eigen::MatrixXd stiffnessMatrix() const override;
        [[nodiscard]] Eigen::MatrixXd massMatrix() const override;

        std::ostream& print(std::ostream& os) const override;
        // Large Deformation
        [[nodiscard]] Eigen::VectorXd internalForce(const Eigen::VectorXd& ue) const override {
            return stiffnessMatrix() * ue;
        }

        [[nodiscard]] Eigen::MatrixXd tangentStiffness(const Eigen::VectorXd& ue) const override {
            (void)ue;
            return stiffnessMatrix();
        }
};

class Quad4 : public Element {

    public:
        double t;

        Quad4(unsigned int id, std::vector<unsigned int> connectivity, std::vector<std::shared_ptr<Node>> Nodes, std::shared_ptr<Material> mat, double Thickness): Element(id, std::move(connectivity), std::move(Nodes) ,std::move(mat)){
            if (conn.size() != 4) {throw std::runtime_error("Quad4(): Node number should be 4.");}

            if (Thickness > 0.0) {
                t = Thickness;
            } else {
                throw std::runtime_error("Quad4(): Thickness should be positive.");
            }
        }

        [[nodiscard]] int numNodes() const override { return 4; }
        [[nodiscard]] int dimension() const override { return 2; }
        [[nodiscard]] std::string name() const override { return "Quad4"; }

        [[nodiscard]] double ShapeFunction(int i, double r, double s) const override;
        [[nodiscard]] double dShapeFunction_dr(int i, double r, double s) const override;
        [[nodiscard]] double dShapeFunction_ds(int i, double r, double s) const override;
        [[nodiscard]] Eigen::MatrixXd J( double r, double s) const override;

        [[nodiscard]] Eigen::MatrixXd stiffnessMatrix() const override;
        [[nodiscard]] Eigen::MatrixXd massMatrix() const override;

        std::ostream& print(std::ostream& os) const override;

        // Large Deformation
        [[nodiscard]] Eigen::Matrix2d displacementGradient(const Eigen::VectorXd& U) const;
        [[nodiscard]] Eigen::Matrix2d deformationGradient(const Eigen::VectorXd& U) const;
        [[nodiscard]] Eigen::VectorXd internalForce(const Eigen::VectorXd& U) const override;
        [[nodiscard]] Eigen::MatrixXd tangentStiffness(const Eigen::VectorXd& U) const override;
};

inline std::ostream& operator<<(std::ostream& os, const Element& el) {return el.print(os);}

//======================================================================================================================
class Geometry {
    public:
        std::vector<std::shared_ptr<Node>> nodes;
        std::vector<std::unique_ptr<Element>> elements;


        Geometry(std::vector<std::shared_ptr<Node>> Node_list,std::vector<std::unique_ptr<Element>> Elements_list): nodes(std::move(Node_list)), elements(std::move(Elements_list)) {}
        ~Geometry() = default;

        [[nodiscard]] const std::shared_ptr<Node>& node(std::size_t i) const { return nodes.at(i); }
        [[nodiscard]] const Element& element(std::size_t i) const { return *elements.at(i); }
        [[nodiscard]] std::size_t numNodes() const { return nodes.size(); }
        [[nodiscard]] std::size_t numElements() const { return elements.size(); }

        Geometry(const Geometry&) = delete;
        Geometry& operator=(const Geometry&) = delete;

        Geometry(Geometry&&) noexcept = default;
        Geometry& operator=(Geometry&&) noexcept = default;
};

enum class SupportType {
    Fixed,
    SlidingX,
    SlidingY
};

class Mesh {
    public:
        std::vector<std::shared_ptr<Node>> nodes;
        std::vector<std::pair<unsigned int, SupportType>> supports;
        std::vector<std::unique_ptr<Element>> elements;

        std::vector<std::vector<int>> node_dofs;
        std::vector<std::vector<int>> element_dofs;

        Eigen::VectorXd U;
        Eigen::VectorXd F;
        Eigen::MatrixXd M;
        Eigen::MatrixXd K;
        Eigen::MatrixXd C;

        unsigned int dofs_per_node = 0;
        unsigned int total_dofs = 0;

        Mesh(std::vector<std::shared_ptr<Node>> node_list,std::vector<std::unique_ptr<Element>> element_list,unsigned int ndpn)
            : nodes(std::move(node_list)),
              elements(std::move(element_list)),
              dofs_per_node(ndpn)
        {
            total_dofs = static_cast<unsigned int>(nodes.size()) * ndpn;
        }

        Mesh(Geometry&& geometry, int ndpn);

        Mesh(const Mesh&) = delete;
        Mesh& operator=(const Mesh&) = delete;
        Mesh(Mesh&&) noexcept = default;
        Mesh& operator=(Mesh&&) noexcept = default;

        static Mesh buildMesh(Geometry&& geometry,int ndpn,double meshDensity,const std::shared_ptr<Material>& material,double thickness);
        void buildElementDofs(unsigned int dofsPerNode);
        ~Mesh() = default;

        void setNodalLoad(unsigned int node_id, const std::vector<double>& Vector_value);
        void addNodalLoad(unsigned int node_id, const std::vector<double>& Vector_value);
        [[nodiscard]] std::vector<unsigned int> findNodesOnEdge(char edge, double tol) const;

        void addSupport(unsigned int node_id, SupportType type);
        void addSupportOnEdge(char edge, SupportType type, double tol = 1e-9);

        [[nodiscard]] std::vector<int> getSupportDofs() const;
        [[nodiscard]] int findNodeByCoordinates(double x, double y, double tol) const;
        Eigen::MatrixXd globalStiffnessMatrix();
        Eigen::MatrixXd globalDumpingMatrix_RayleighMethode(double f1, double f2, double zeta);
        Eigen::MatrixXd globalMassMatrix();

        // Large Deformationd
        [[nodiscard]] Eigen::VectorXd globalInternalForce(const Eigen::VectorXd& Uglob) const;
        [[nodiscard]] Eigen::MatrixXd globalTangentStiffness(const Eigen::VectorXd& Uglob) const;
};

#endif //Physics_H
