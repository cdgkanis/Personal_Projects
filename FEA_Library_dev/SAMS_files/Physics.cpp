//
// Created by User on 30/4/2026.
//

#include "Physics.h"

#include <vector>

//======================================================================================================================
static Eigen::MatrixXd removeRowsCols(const Eigen::MatrixXd& A,const std::vector<int>& fixed_dofs){
    const int n = static_cast<int>(A.rows());
    std::vector<bool> isFixed(n, false);

    for (int dof : fixed_dofs) {
        if (dof < 0 || dof >= n) {
            throw std::out_of_range("removeRowsCols(): DOF index out of range.");
        }
        isFixed[dof] = true;
    }

    std::vector<int> free_dofs;
    free_dofs.reserve(n - static_cast<int>(fixed_dofs.size()));

    for (int i = 0; i < n; ++i) {
        if (!isFixed[i]) {
            free_dofs.push_back(i);
        }
    }

    Eigen::MatrixXd R(free_dofs.size(), free_dofs.size());

    for (int i = 0; i < static_cast<int>(free_dofs.size()); ++i) {
        for (int j = 0; j < static_cast<int>(free_dofs.size()); ++j) {
            R(i, j) = A(free_dofs[i], free_dofs[j]);
        }
    }

    return R;
}

static std::string toLowerCopy(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) {
            return static_cast<char>(std::tolower(c));
        });
    return s;
}

//======================================================================================================================
// Isotropic material functions
std::ostream& IsotropicMaterial::print(std::ostream& os) const {
    os << "Material properties\n";
    os << "Name : " << Name << '\n';
    os << "Type : Isotropic\n";
    os << "rho  : " << rho_ << '\n';
    os << "E    : " << E_ << '\n';
    os << "nu   : " << nu_ << '\n';
    return os;
}

// Orthotropic Material functions
std::ostream& OrthotropicMaterial::print(std::ostream& os) const {
    os << "Material properties\n";
    os << "Name : " << Name << '\n';
    os << "Type : Orthotropic\n";
    os << "rho  : " << rho_ << '\n';
    os << "Ex   : " << Ex_ << '\n';
    os << "Ey   : " << Ey_ << '\n';
    os << "Ez   : " << Ez_ << '\n';
    os << "nuxy : " << nuxy_ << '\n';
    os << "nuyz : " << nuyz_ << '\n';
    os << "nuxz : " << nuxz_ << '\n';
    os << "Gxy  : " << Gxy_ << '\n';
    os << "Gyz  : " << Gyz_ << '\n';
    os << "Gxz  : " << Gxz_ << '\n';
    return os;
}

//======================================================================================================================
bool Element::HasAssumption(AssumptionType assumption) const {
    return std::find(assumptions.begin(), assumptions.end(), assumption) != assumptions.end();
}

void Element::AddAssumption(AssumptionType assumption) {
    if (!HasAssumption(assumption)) {
        assumptions.push_back(assumption);
    }
}

void Element::RemoveAssumption(AssumptionType assumption) {
    if (HasAssumption(assumption)) {
        auto it = std::find(assumptions.begin(), assumptions.end(), assumption);
        assumptions.erase(it);
    }
}

void Element::AddAssumption(const std::string& assumption) {
    if (assumption.empty()) {
        throw std::invalid_argument("Element::AddAssumption(): assumption is empty.");
    }

    const std::string a = toLowerCopy(assumption);

    if (a == "planestress") {
        AddAssumption(AssumptionType::PlaneStress);
    }
    else if (a == "planestrain") {
        AddAssumption(AssumptionType::PlaneStrain);
    }
    else if (a == "lumped") {
        AddAssumption(AssumptionType::Lumped);
    }
    else if (a == "consistent") {
        AddAssumption(AssumptionType::Consistent);
    }
    else {
        throw std::invalid_argument("Element::AddAssumption(): invalid assumption.");
    }
}

// Tri3 functions
std::ostream & Tri3::print(std::ostream &os) const  {
    os << "Element ID : " << ID << '\n'
       << "Type       : Tri3\n"
       << "Material   : " << material->name() << '\n'
       << "Dimension  : 2D\n"
       << "Thickness  : " << t << " [m]\n"
       << "Nodes      : ";
    for (unsigned int n : conn) os << n << ' ';
    os << '\n';
    return os;
}

double Tri3::ShapeFunction(int i, double r, double s) const {
    if (i == 0) return 1.0 - r - s;
    if (i == 1) return r;
    if (i == 2) return s;
    throw std::out_of_range("Tri3::ShapeFunction(): invalid i");
}

double Tri3::dShapeFunction_dr(int i, double r, double s) const {
    (void)r;
    (void)s;

    switch (i) {
        case 0: return -1.0;
        case 1: return  1.0;
        case 2: return  0.0;
        default:
            throw std::out_of_range("Tri3::dShapeFunction_dr(): i must be 0, 1 or 2.");
    }
}

double Tri3::dShapeFunction_ds(int i, double r, double s) const {
    (void)r;
    (void)s;

    switch (i) {
        case 0: return -1.0;
        case 1: return  0.0;
        case 2: return  1.0;
        default:
            throw std::out_of_range("Tri3::dShapeFunction_ds(): i must be 0, 1 or 2.");
    }
}

Eigen::MatrixXd Tri3::J(double r, double s) const {
    (void)r;
    (void)s;

    if (nodes.size() != 3) {
        throw std::runtime_error("Tri3::J(): Tri3 must have exactly 3 nodes.");
    }

    const double x1 = nodes[0]->x[0];
    const double y1 = nodes[0]->x[1];
    const double x2 = nodes[1]->x[0];
    const double y2 = nodes[1]->x[1];
    const double x3 = nodes[2]->x[0];
    const double y3 = nodes[2]->x[1];

    Eigen::Matrix2d Jmat;
    Jmat(0,0) = x2 - x1;   // dx/dr
    Jmat(0,1) = x3 - x1;   // dx/ds
    Jmat(1,0) = y2 - y1;   // dy/dr
    Jmat(1,1) = y3 - y1;   // dy/ds

    return Jmat;
}

Eigen::MatrixXd Tri3::stiffnessMatrix() const {
    if (nodes.size() != 3) {
        throw std::runtime_error("Tri3::stiffnessMatrix(): Tri3 must have exactly 3 nodes.");
    }

    if (material->type() != MaterialType::Isotropic) {
        throw std::runtime_error("Tri3::stiffnessMatrix(): Material has to be isotropic.");
    }

    auto iso = std::dynamic_pointer_cast<IsotropicMaterial>(material);
    if (!iso) {
        throw std::runtime_error("Tri3::stiffnessMatrix(): Failed to cast material to IsotropicMaterial.");
    }

    const double E  = iso->E();
    const double nu = iso->nu();

    Eigen::Matrix2d Jmat = J(0.0, 0.0);
    const double detJ = Jmat.determinant();

    if (detJ <= 0.0) {
        throw std::runtime_error("Tri3::stiffnessMatrix(): det(J) <= 0. Check node ordering; Tri3 should be counter-clockwise.");
    }

    const double A = 0.5 * detJ;

    Eigen::Matrix<double, 2, 3> dNdrs;
    dNdrs <<
        -1.0,  1.0,  0.0,
        -1.0,  0.0,  1.0;

    Eigen::Matrix2d invJ = Jmat.inverse();
    Eigen::Matrix<double, 2, 3> dNdxy = invJ * dNdrs;

    Eigen::Matrix<double, 3, 6> B;
    B.setZero();

    for (int i = 0; i < 3; ++i) {
        const double dNdx = dNdxy(0, i);
        const double dNdy = dNdxy(1, i);

        B(0, 2*i    ) = dNdx;
        B(1, 2*i + 1) = dNdy;
        B(2, 2*i    ) = dNdy;
        B(2, 2*i + 1) = dNdx;
    }

    Eigen::Matrix3d D;
    D.setZero();

    if (this->HasAssumption(AssumptionType::PlaneStrain)) {
        const double c = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
        D <<
            c * (1.0 - nu), c * nu,           0.0,
            c * nu,         c * (1.0 - nu),   0.0,
            0.0,            0.0,              c * (1.0 - 2.0 * nu) / 2.0;
    }
    else if(this->HasAssumption(AssumptionType::PlaneStress)){ // default: PlaneStress
        const double c = E / (1.0 - nu * nu);
        D <<
            c,         c * nu,    0.0,
            c * nu,    c,         0.0,
            0.0,       0.0,       E / (2.0 * (1.0 + nu));
    }else {
        throw std::runtime_error("Tri3::stiffnessMatrix(): Assumption type is not supported.");
    }

    Eigen::Matrix<double, 6, 6> Ke = t * A * (B.transpose() * D * B);
    return Ke;
}

Eigen::MatrixXd Tri3::massMatrix() const {
    if (nodes.size() != 3) {
        throw std::runtime_error("Tri3::massMatrix(): Tri3 must have exactly 3 nodes.");
    }

    const double rho = material->rho();

    Eigen::Matrix<double, 6, 6> Me;
    Me.setZero();

    // 1-point Gauss rule for linear triangle: (1/3, 1/3)
    const double r = 1.0 / 3.0;
    const double s = 1.0 / 3.0;

    Eigen::Matrix2d Jmat = J(r, s);
    const double detJ = Jmat.determinant();

    if (detJ <= 0.0) {
        throw std::runtime_error("Tri3::massMatrix(): det(J) <= 0. Check node ordering.");
    }

    const double N1 = ShapeFunction(0, r, s);  // 1 - r - s
    const double N2 = ShapeFunction(1, r, s);  // r
    const double N3 = ShapeFunction(2, r, s);  // s

    Eigen::Matrix<double, 2, 6> N;
    N.setZero();

    N(0,0) = N1;  N(0,2) = N2;  N(0,4) = N3;
    N(1,1) = N1;  N(1,3) = N2;  N(1,5) = N3;

    Me = rho * t * (N.transpose() * N) * detJ;

    if (HasAssumption(AssumptionType::Lumped)) {
        Eigen::Matrix<double, 6, 6> Ml;
        Ml.setZero();

        for (int i = 0; i < 6; ++i) {
            Ml(i, i) = Me.row(i).sum();
        }

        return Ml;
    }

    if (HasAssumption(AssumptionType::Consistent) ||
        (!HasAssumption(AssumptionType::Lumped) && !HasAssumption(AssumptionType::Consistent))) {
        return Me;
        }

    throw std::runtime_error("Tri3::massMatrix(): invalid mass assumption.");
}


// Quad4 functions
std::ostream & Quad4::print(std::ostream &os) const  {
    os << "Element ID : " << ID << '\n'
       << "Type       : Quad4\n"
       << "Material   : " << material->name() << '\n'
       << "Dimension  : 2D\n"
       << "Thickness  : " << t << " [m]\n"
       << "Nodes      : ";
    for (unsigned int n : conn) os << n << ' ';
    os << '\n';
    return os;
}

double Quad4::ShapeFunction(int i, double r, double s) const{
    if (i < 0 || i > 3) {
        throw std::invalid_argument("ShapeFunction(): Number of shape function is out of bounds.");
    }

    const int r_i[4] = {-1,  1,  1, -1};
    const int s_i[4] = {-1, -1,  1,  1};

    return 0.25 * (1.0 + r_i[i] * r) * (1.0 + s_i[i] * s);
}

double Quad4::dShapeFunction_dr(int i, double r, double s) const{
    if (i < 0 || i > 3) {
        throw std::invalid_argument("dShapeFunction_dr(): Number of shape function is out of bounds.");
    }

    const int r_i[4] = {-1,  1,  1, -1};
    const int s_i[4] = {-1, -1,  1,  1};

    return 0.25 * r_i[i] * (1.0 + s_i[i] * s);
}

double Quad4::dShapeFunction_ds(int i, double r, double s) const{
    if (i < 0 || i > 3) {
        throw std::invalid_argument("dShapeFunction_ds(): Number of shape function is out of bounds.");
    }

    const int r_i[4] = {-1,  1,  1, -1};
    const int s_i[4] = {-1, -1,  1,  1};

    return 0.25 * s_i[i] * (1.0 + r_i[i] * r);
}

Eigen::MatrixXd Quad4::J(double r, double s) const {
    Eigen::Matrix2d Jac = Eigen::Matrix2d::Zero();

    for (int a = 0; a < 4; ++a) {
        const Node& node = *this->nodes.at(a);

        const double x = node.x.at(0);
        const double y = node.x.at(1);

        const double dNdr = dShapeFunction_dr(a, r, s);
        const double dNds = dShapeFunction_ds(a, r, s);

        Jac(0,0) += dNdr * x;
        Jac(0,1) += dNds * x;
        Jac(1,0) += dNdr * y;
        Jac(1,1) += dNds * y;
    }

    return Jac;
}

Eigen::MatrixXd Quad4::stiffnessMatrix() const {
    if (nodes.size() != 4) {
        throw std::runtime_error("Quad4::stiffnessMatrix(): Quad4 must have exactly 4 nodes.");
    }

    if (material->type() != MaterialType::Isotropic) {
        throw std::runtime_error("Quad4::stiffnessMatrix(): Material has to be isotropic.");
    }

    auto iso = std::dynamic_pointer_cast<IsotropicMaterial>(material);
    if (!iso) {
        throw std::runtime_error("Quad4::stiffnessMatrix(): Failed to cast material to IsotropicMaterial.");
    }

    const double E  = iso->E();
    const double nu = iso->nu();

    Eigen::Matrix3d D;
    D.setZero();

    if (this->HasAssumption(AssumptionType::PlaneStrain)) {
        const double c = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
        D <<
            c * (1.0 - nu), c * nu,           0.0,
            c * nu,         c * (1.0 - nu),   0.0,
            0.0,            0.0,              c * (1.0 - 2.0 * nu) / 2.0;
    } else if (this->HasAssumption(AssumptionType::PlaneStress)){
        const double c = E / (1.0 - nu * nu);
        D <<
            c,         c * nu,    0.0,
            c * nu,    c,         0.0,
            0.0,       0.0,       E / (2.0 * (1.0 + nu));
    }else {
        throw std::runtime_error("Quad4::stiffnessMatrix(): Assumption does not match with the problem.");
    }

    Eigen::Matrix<double, 8, 8> Ke;
    Ke.setZero();

    //Gauss Nodes #2
    const double g = 1.0 / std::sqrt(3.0);
    const double gp[] = {-g, g};

    for (double r : gp) {
        for (double s : gp) {
            Eigen::Matrix2d Jmat = J(r, s);
            const double detJ = Jmat.determinant();

            if (detJ <= 0.0) {
                throw std::runtime_error("Quad4::stiffnessMatrix(): det(J) <= 0. Check node ordering/distortion.");
            }

            Eigen::Matrix2d invJ = Jmat.inverse();

            Eigen::Matrix<double, 2, 4> dNdrs;
            for (int i = 0; i < 4; ++i) {
                dNdrs(0, i) = dShapeFunction_dr(i, r, s);
                dNdrs(1, i) = dShapeFunction_ds(i, r, s);
            }

            Eigen::Matrix<double, 2, 4> dNdxy = invJ * dNdrs;

            Eigen::Matrix<double, 3, 8> B;
            B.setZero();

            for (int i = 0; i < 4; ++i) {
                const double dNdx = dNdxy(0, i);
                const double dNdy = dNdxy(1, i);

                B(0, 2 * i    ) = dNdx;
                B(1, 2 * i + 1) = dNdy;
                B(2, 2 * i    ) = dNdy;
                B(2, 2 * i + 1) = dNdx;
            }

            Ke += t * (B.transpose() * D * B) * detJ;
        }
    }



    return Ke;
}

Eigen::MatrixXd Quad4::massMatrix() const {
    if (nodes.size() != 4) {
        throw std::runtime_error("Quad4::massMatrix(): Quad4 must have exactly 4 nodes.");
    }

    const double rho = material->rho();

    Eigen::Matrix<double, 8, 8> Me;
    Me.setZero();

    // 2x2 Gauss integration
    const double g = 1.0 / std::sqrt(3.0);
    const double gp[] = { -g, g };

    for (double r : gp) {
        for (double s : gp) {
            Eigen::Matrix2d Jmat = J(r, s);
            const double detJ = Jmat.determinant();

            if (detJ <= 0.0) {
                throw std::runtime_error("Quad4::massMatrix(): det(J) <= 0. Check node ordering/distortion.");
            }

            const double N1 = ShapeFunction(0, r, s);
            const double N2 = ShapeFunction(1, r, s);
            const double N3 = ShapeFunction(2, r, s);
            const double N4 = ShapeFunction(3, r, s);

            Eigen::Matrix<double, 2, 8> N;
            N.setZero();

            N(0,0) = N1;  N(0,2) = N2;  N(0,4) = N3;  N(0,6) = N4;
            N(1,1) = N1;  N(1,3) = N2;  N(1,5) = N3;  N(1,7) = N4;

            Me += rho * t * (N.transpose() * N) * detJ;
        }
    }

    if (HasAssumption(AssumptionType::Lumped)) {
        Eigen::Matrix<double, 8, 8> Ml;
        Ml.setZero();

        for (int i = 0; i < 8; ++i) {
            Ml(i, i) = Me.row(i).sum();
        }

        return Ml;
    }

    if (HasAssumption(AssumptionType::Consistent)) {
        return Me;
    }

    throw std::runtime_error("Quad4::massMatrix(): invalid mass assumption.");
}

// Large Deformation
Eigen::Matrix2d Quad4::displacementGradient(const Eigen::VectorXd& U) const {
    if (this->nodes.size() != 4) {
        throw std::runtime_error("Quad4::displacementGradient(): Quad4 must have exactly 4 nodes.");
    }

    if (U.size() != 8) {
        throw std::runtime_error("Quad4::displacementGradient(): displacement vector must have size 8.");
    }

    Eigen::Matrix<double, 4, 2> dNdrs;
    dNdrs.setZero();

    for (int i = 0; i < 4; ++i) {
        dNdrs(i, 0) = dShapeFunction_dr(i, 0.0, 0.0);
        dNdrs(i, 1) = dShapeFunction_ds(i, 0.0, 0.0);
    }

    Eigen::Matrix2d Jmat = J(0.0, 0.0);
    const double detJ = Jmat.determinant();

    if (detJ <= 0.0) {
        throw std::runtime_error("Quad4::displacementGradient(): det(J) <= 0.");
    }

    Eigen::Matrix2d invJ = Jmat.inverse();
    Eigen::Matrix<double, 4, 2> dNdxy = dNdrs * invJ.transpose();

    double uX = 0.0, uY = 0.0, vX = 0.0, vY = 0.0;

    for (int a = 0; a < 4; ++a) {
        const double ua = U(2 * a);
        const double va = U(2 * a + 1);

        const double dNdx = dNdxy(a, 0);
        const double dNdy = dNdxy(a, 1);

        uX += dNdx * ua;
        uY += dNdy * ua;
        vX += dNdx * va;
        vY += dNdy * va;
    }

    Eigen::Matrix2d gradU;
    gradU << uX, uY,
             vX, vY;

    return gradU;
}

Eigen::Matrix2d Quad4::deformationGradient(const Eigen::VectorXd &U) const {
    return Eigen::Matrix2d::Identity() + displacementGradient(U);
}

Eigen::VectorXd Quad4::internalForce(const Eigen::VectorXd& U) const {
    if (U.size() != 8) {
        throw std::runtime_error("Quad4::internalForce(): displacement vector must have size 8.");
    }

    if (material->type() != MaterialType::Isotropic) {
        throw std::runtime_error("Quad4::internalForce(): Material has to be isotropic.");
    }

    auto iso = std::dynamic_pointer_cast<IsotropicMaterial>(material);
    if (!iso) {
        throw std::runtime_error("Quad4::internalForce(): Failed to cast material to IsotropicMaterial.");
    }

    const double E  = iso->E();
    const double nu = iso->nu();

    Eigen::Matrix3d D;
    D.setZero();

    if (this->HasAssumption(AssumptionType::PlaneStrain)) {
        const double c = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
        D <<
            c * (1.0 - nu), c * nu,           0.0,
            c * nu,         c * (1.0 - nu),   0.0,
            0.0,            0.0,              c * (1.0 - 2.0 * nu) / 2.0;
    } else if (this->HasAssumption(AssumptionType::PlaneStress)) {
        const double c = E / (1.0 - nu * nu);
        D <<
            c,         c * nu,    0.0,
            c * nu,    c,         0.0,
            0.0,       0.0,       E / (2.0 * (1.0 + nu));
    } else {
        throw std::runtime_error("Quad4::internalForce(): Assumption does not match with the problem.");
    }

    Eigen::Matrix<double, 8, 1> fint;
    fint.setZero();

    const double g = 1.0 / std::sqrt(3.0);
    const double gp[] = { -g, g };

    for (double r : gp) {
        for (double s : gp) {
            Eigen::Matrix2d Jmat = J(r, s);
            const double detJ = Jmat.determinant();

            if (detJ <= 0.0) {
                throw std::runtime_error("Quad4::internalForce(): det(J) <= 0. Check node ordering/distortion.");
            }

            Eigen::Matrix2d invJ = Jmat.inverse();

            Eigen::Matrix<double, 2, 4> dNdrs;
            for (int i = 0; i < 4; ++i) {
                dNdrs(0, i) = dShapeFunction_dr(i, r, s);
                dNdrs(1, i) = dShapeFunction_ds(i, r, s);
            }

            Eigen::Matrix<double, 2, 4> dNdxy = invJ * dNdrs;


            Eigen::Matrix<double, 4, 8> B_nabla;
            B_nabla.setZero();

            Eigen::Matrix<double, 3, 8> B;
            B.setZero();

            Eigen::Matrix<double, 3, 4> A_L;
            A_L.setZero();

            for (int i = 0; i < 4; ++i) {
                const double ui = U(2 * i);
                const double vi = U(2 * i + 1);

                const double dNdx = dNdxy(0, i);
                const double dNdy = dNdxy(1, i);


                B_nabla(0, 2 * i    ) = dNdx;
                B_nabla(1, 2 * i + 1) = dNdy;
                B_nabla(2, 2 * i    ) = dNdy;
                B_nabla(3, 2 * i + 1) = dNdx;

                B(0, 2 * i    ) = dNdx;
                B(1, 2 * i + 1) = dNdy;
                B(2, 2 * i    ) = dNdy;
                B(2, 2 * i + 1) = dNdx;
            }

            Eigen::Matrix2d F = deformationGradient(U);
            double det_F = F.determinant();
            if (det_F == 0.0) {
                throw std::runtime_error("Quad4::internalForce(): F determinant() was zero.");
            }else if(det_F < 0.0) {
                throw std::runtime_error("Quas4::internalForce(): F determinant() was negative.");
            }

            A_L(0, 0) = F(0, 0);
            A_L(1, 1) = F(1, 1);
            A_L(2, 0) = F(0, 1);
            A_L(2, 1) = F(1, 0);

            A_L(1, 2) = F(0, 1);
            A_L(2, 2) = F(0, 0);
            A_L(0, 3) = F(1, 0);
            A_L(2, 3) = F(1, 1);


            Eigen::Matrix<double, 3, 8> B_L = A_L * B_nabla;

            const Eigen::Vector3d strain = B * U;
            const Eigen::Vector3d stress = D * strain;

            Eigen::Matrix2d sigma;
            sigma << stress(0), stress(2),
                     stress(2), stress(1);

            const Eigen::Matrix2d F_inv = F.inverse();

            const Eigen::Matrix2d P_1st_tensor = det_F * sigma * F_inv.transpose() ;
            const Eigen::Matrix2d S_2nd_tensor = F_inv * P_1st_tensor;

            Eigen::Vector3d S_2nd;
            S_2nd << S_2nd_tensor( 0, 0), S_2nd_tensor( 1, 1), S_2nd_tensor( 0, 1);

            fint += t * (B_L.transpose() * S_2nd)* detJ;
        }
    }

    return fint;
}

Eigen::MatrixXd Quad4::tangentStiffness(const Eigen::VectorXd &U) const {
    if (U.size() != 8) {
        throw std::runtime_error("Quad4::internalForce(): displacement vector must have size 8.");
    }

    if (material->type() != MaterialType::Isotropic) {
        throw std::runtime_error("Quad4::internalForce(): Material has to be isotropic.");
    }

    auto iso = std::dynamic_pointer_cast<IsotropicMaterial>(material);
    if (!iso) {
        throw std::runtime_error("Quad4::internalForce(): Failed to cast material to IsotropicMaterial.");
    }

    const double E  = iso->E();
    const double nu = iso->nu();

    Eigen::Matrix3d D;
    D.setZero();

    if (this->HasAssumption(AssumptionType::PlaneStrain)) {
        const double c = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
        D <<
            c * (1.0 - nu), c * nu,           0.0,
            c * nu,         c * (1.0 - nu),   0.0,
            0.0,            0.0,              c * (1.0 - 2.0 * nu) / 2.0;
    } else if (this->HasAssumption(AssumptionType::PlaneStress)) {
        const double c = E / (1.0 - nu * nu);
        D <<
            c,         c * nu,    0.0,
            c * nu,    c,         0.0,
            0.0,       0.0,       E / (2.0 * (1.0 + nu));
    } else {
        throw std::runtime_error("Quad4::internalForce(): Assumption does not match with the problem.");
    }

    const double g = 1.0 / std::sqrt(3.0);
    const double gp[] = { -g, g };

    Eigen::MatrixXd K_L  = Eigen::MatrixXd::Zero(8, 8);
    Eigen::MatrixXd K_NL = Eigen::MatrixXd::Zero(8, 8);

    for (double r : gp) {
        for (double s : gp) {
            Eigen::Matrix2d Jmat = J(r, s);
            const double detJ = Jmat.determinant();

            if (detJ <= 0.0) {
                throw std::runtime_error("Quad4::internalForce(): det(J) <= 0. Check node ordering/distortion.");
            }

            Eigen::Matrix2d invJ = Jmat.inverse();

            Eigen::Matrix<double, 2, 4> dNdrs;
            for (int i = 0; i < 4; ++i) {
                dNdrs(0, i) = dShapeFunction_dr(i, r, s);
                dNdrs(1, i) = dShapeFunction_ds(i, r, s);
            }

            Eigen::Matrix<double, 2, 4> dNdxy = invJ * dNdrs;


            Eigen::Matrix<double, 4, 8> B_nabla;
            B_nabla.setZero();

            Eigen::Matrix<double, 3, 8> B;
            B.setZero();

            Eigen::Matrix<double, 3, 4> A_L;
            A_L.setZero();

            for (int i = 0; i < 4; ++i) {
                const double dNdx = dNdxy(0, i);
                const double dNdy = dNdxy(1, i);


                B_nabla(0, 2 * i    ) = dNdx;
                B_nabla(1, 2 * i + 1) = dNdy;
                B_nabla(2, 2 * i    ) = dNdy;
                B_nabla(3, 2 * i + 1) = dNdx;

                B(0, 2 * i    ) = dNdx;
                B(1, 2 * i + 1) = dNdy;
                B(2, 2 * i    ) = dNdy;
                B(2, 2 * i + 1) = dNdx;
            }

            Eigen::Matrix2d F = deformationGradient(U);
            double det_F = F.determinant();
            if (det_F == 0.0) {
                throw std::runtime_error("Quad4::internalForce(): F determinant() was zero.");
            }else if(det_F < 0.0) {
                throw std::runtime_error("Quas4::internalForce(): F determinant() was negative.");
            }

            A_L(0, 0) = F(0, 0);
            A_L(1, 1) = F(1, 1);
            A_L(2, 0) = F(0, 1);
            A_L(2, 1) = F(1, 0);

            A_L(1, 2) = F(0, 1);
            A_L(2, 2) = F(0, 0);
            A_L(0, 3) = F(1, 0);
            A_L(2, 3) = F(1, 1);


            Eigen::Matrix<double, 3, 8> B_L = A_L * B_nabla;

            const Eigen::Vector3d strain = B * U;
            const Eigen::Vector3d stress = D * strain;

            Eigen::Matrix2d sigma;
            sigma << stress(0), stress(2),
                     stress(2), stress(1);

            const Eigen::Matrix2d F_inv = F.inverse();

            const Eigen::Matrix2d P_1st_tensor = det_F * sigma * F_inv.transpose() ;
            const Eigen::Matrix2d S_2nd_tensor = F_inv * P_1st_tensor;

            Eigen::Vector3d S_2nd;
            S_2nd << S_2nd_tensor( 0, 0), S_2nd_tensor( 1, 1), S_2nd_tensor( 0, 1);

            Eigen::Matrix<double, 4, 4> S_2nd_bar;
            S_2nd_bar.setZero();

            S_2nd_bar(0, 0) = S_2nd(0);
            S_2nd_bar(2, 0) = S_2nd(2);
            S_2nd_bar(1, 1) = S_2nd(1);
            S_2nd_bar(3, 1) = S_2nd(2);

            S_2nd_bar(0, 2) = S_2nd(2);
            S_2nd_bar(2, 2) = S_2nd(1);
            S_2nd_bar(1, 3) = S_2nd(2);
            S_2nd_bar(3, 3) = S_2nd(0);

            K_L  += t * (B_L.transpose() * D * B_L) * detJ;
            K_NL += t * (B_nabla.transpose() * S_2nd_bar * B_nabla) * detJ;
        }
    }

    return K_L + K_NL;
}

//======================================================================================================================
Mesh::Mesh(Geometry&& geometry, int ndpn)
    : nodes(std::move(geometry.nodes)),
      elements(std::move(geometry.elements)),
      dofs_per_node(ndpn){
    total_dofs = static_cast<unsigned int>(nodes.size()) * dofs_per_node;

    node_dofs.resize(nodes.size());
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        node_dofs[i].resize(dofs_per_node);
        for (unsigned int d = 0; d < dofs_per_node; ++d) {
            node_dofs[i][d] = static_cast<int>(i * dofs_per_node + d);
        }
    }

    element_dofs.resize(elements.size());
    for (std::size_t e = 0; e < elements.size(); ++e) {
        const auto& conn = elements[e]->connectivity();
        for (unsigned int node_id : conn) {
            for (unsigned int d = 0; d < dofs_per_node; ++d) {
                element_dofs[e].push_back(node_dofs[node_id][d]);
            }
        }
    }

    U = Eigen::VectorXd::Zero(total_dofs);
    F = Eigen::VectorXd::Zero(total_dofs);
}

Mesh Mesh::buildMesh(Geometry&& geometry,int ndpn,double meshDensity,const std::shared_ptr<Material>& material, double thickness){
    if (geometry.numNodes() == 0 || geometry.numElements() == 0) {
        throw std::invalid_argument("Mesh::buildMesh(): empty geometry.");
    }

    if (ndpn <= 0) {
        throw std::invalid_argument("Mesh::buildMesh(): dofs_per_node must be positive.");
    }

    if (meshDensity <= 0.0) {
        throw std::invalid_argument("Mesh::buildMesh(): meshDensity must be positive.");
    }

    double x_min = +std::numeric_limits<double>::max();
    double x_max = -std::numeric_limits<double>::max();
    double y_min = +std::numeric_limits<double>::max();
    double y_max = -std::numeric_limits<double>::max();

    for (const auto& node : geometry.nodes) {
        if (node->x.size() < 2) {
            throw std::runtime_error("Mesh::buildMesh(): node coordinate vector too small.");
        }

        x_min = std::min(x_min, node->x[0]);
        x_max = std::max(x_max, node->x[0]);
        y_min = std::min(y_min, node->x[1]);
        y_max = std::max(y_max, node->x[1]);
    }

    const double Lx = x_max - x_min;
    const double Ly = y_max - y_min;

    if (Lx <= 0.0 || Ly <= 0.0) {
        throw std::runtime_error("Mesh::buildMesh(): domain has zero size.");
    }

    const int nx = std::max(1, static_cast<int>(std::round(meshDensity * Lx)));
    const int ny = std::max(1, static_cast<int>(std::round(meshDensity * Ly)));

    const int nNodesX = nx + 1;
    const int nNodesY = ny + 1;
    const int nNodes  = nNodesX * nNodesY;
    const int nElems  = nx * ny;

    std::vector<std::shared_ptr<Node>> new_nodes;
    std::vector<std::unique_ptr<Element>> new_elements;

    new_nodes.reserve(nNodes);
    new_elements.reserve(nElems);

    const double dx = Lx / static_cast<double>(nx);
    const double dy = Ly / static_cast<double>(ny);

    for (int iy = 0; iy < nNodesY; ++iy) {
        for (int ix = 0; ix < nNodesX; ++ix) {
            auto node = std::make_shared<Node>();
            node->ID = static_cast<unsigned int>(new_nodes.size());
            node->x = {x_min + ix * dx, y_min + iy * dy};
            new_nodes.push_back(std::move(node));
        }
    }
    unsigned int elemId = 0;

    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            const unsigned int n0 = static_cast<unsigned int>(j * nNodesX + i);
            const unsigned int n1 = n0 + 1;
            const unsigned int n2 = n0 + nNodesX + 1;
            const unsigned int n3 = n0 + nNodesX;

            std::vector<unsigned int> conn = {n0, n1, n2, n3};
            std::vector<std::shared_ptr<Node>> quadNodes = {
                new_nodes[n0], new_nodes[n1], new_nodes[n2], new_nodes[n3]
            };

            new_elements.push_back(
                std::make_unique<Quad4>(
                    elemId++,
                    conn,
                    std::move(quadNodes),
                    material,
                    thickness
                )
            );
        }
    }

    Mesh mesh(std::move(new_nodes), std::move(new_elements), static_cast<unsigned int>(ndpn));

    mesh.node_dofs.resize(mesh.nodes.size());
    for (std::size_t i = 0; i < mesh.nodes.size(); ++i) {
        mesh.node_dofs[i].resize(mesh.dofs_per_node);
        for (unsigned int d = 0; d < mesh.dofs_per_node; ++d) {
            mesh.node_dofs[i][d] = static_cast<int>(i * mesh.dofs_per_node + d);
        }
    }

    mesh.element_dofs.resize(mesh.elements.size());
    for (std::size_t e = 0; e < mesh.elements.size(); ++e) {
        const auto& conn = mesh.elements[e]->connectivity();
        for (unsigned int node_id : conn) {
            for (unsigned int d = 0; d < mesh.dofs_per_node; ++d) {
                mesh.element_dofs[e].push_back(mesh.node_dofs[node_id][d]);
            }
        }
    }

    mesh.total_dofs = static_cast<unsigned int>(mesh.nodes.size()) * mesh.dofs_per_node;

    mesh.U = Eigen::VectorXd::Zero(mesh.total_dofs);
    mesh.F = Eigen::VectorXd::Zero(mesh.total_dofs);
    mesh.M = Eigen::MatrixXd::Zero(mesh.total_dofs, mesh.total_dofs);
    mesh.K = Eigen::MatrixXd::Zero(mesh.total_dofs, mesh.total_dofs);
    mesh.C = Eigen::MatrixXd::Zero(mesh.total_dofs, mesh.total_dofs);

    return mesh;
}

void Mesh::addNodalLoad(unsigned int node_id, const std::vector<double>& Vector_value) {
    if (node_id >= this->node_dofs.size()) {
        throw std::out_of_range("addNodalLoad(): Node ID out of range.");
    }

    if (Vector_value.size() != this->dofs_per_node) {
        throw std::out_of_range("addNodalLoad(): Vector and problem have different dimensions.");
    }

    const std::vector<int>& global_dof = this->node_dofs[node_id];

    for (unsigned int i = 0; i < global_dof.size(); ++i) {
        int dof = global_dof[i];

        if (dof < 0 || dof >= F.size()) {
            throw std::runtime_error("addNodalLoad(): Invalid global DOF index.");
        }

        F(dof) += Vector_value[i];
    }
}

void Mesh::setNodalLoad(unsigned int node_id, const std::vector<double>& Vector_value) {
    if (node_id >= this->node_dofs.size()) {
        throw std::out_of_range("addNodalLoad(): Node ID out of range.");
    }

    if (Vector_value.size() != this->dofs_per_node) {
        throw std::out_of_range("addNodalLoad(): Vector and problem have different dimensions.");
    }

    const std::vector<int>& global_dof = this->node_dofs[node_id];

    //std::cout << global_dof[0] << ' ' << global_dof[1] << '\n';
    for (unsigned int i = 0; i < global_dof.size(); ++i) {
        int dof = global_dof[i];

        if (dof < 0 || dof >= F.size()) {
            throw std::runtime_error("addNodalLoad(): Invalid global DOF index.");
        }

        F(dof) = Vector_value[i];
    }
}

std::vector<int> Mesh::getSupportDofs() const {
    std::vector<int> fixed_dofs;

    for (const auto& s : supports) {
        unsigned int nid = s.first;
        SupportType type = s.second;

        if (nid >= node_dofs.size()) {
            throw std::out_of_range("Mesh::getSupportDofs(): node_id out of range in supports.");
        }

        switch (type) {
            case SupportType::Fixed:
                fixed_dofs.push_back(node_dofs[nid][0]);
            fixed_dofs.push_back(node_dofs[nid][1]);
            break;

            case SupportType::SlidingX:
                fixed_dofs.push_back(node_dofs[nid][1]); // Uy fixed
            break;

            case SupportType::SlidingY:
                fixed_dofs.push_back(node_dofs[nid][0]); // Ux fixed
            break;
        }
    }

    std::sort(fixed_dofs.begin(), fixed_dofs.end());
    fixed_dofs.erase(std::unique(fixed_dofs.begin(), fixed_dofs.end()), fixed_dofs.end());

    return fixed_dofs;
}

void Mesh::addSupportOnEdge(char edge, SupportType type, double tol) {
    std::vector<unsigned int> edge_nodes = findNodesOnEdge(edge, tol);

    for (unsigned int nid : edge_nodes) {
        addSupport(nid, type);
    }
}

void Mesh::addSupport(unsigned int node_id, SupportType type) {
    if (node_id >= nodes.size()) {
        throw std::out_of_range("Mesh::addSupport(): node_id out of range.");
    }

    for (auto& s : supports) {
        if (s.first == node_id) {
            s.second = type;
            return;
        }
    }

    supports.emplace_back(node_id, type);
}

Eigen::MatrixXd Mesh::globalStiffnessMatrix() {

    Eigen::MatrixXd K_ = Eigen::MatrixXd::Zero(total_dofs, total_dofs);
    Eigen::MatrixXd Ke;

    for (std::size_t e = 0; e < elements.size(); ++e) {
        if (elements[e]->HasAssumption(AssumptionType::LargeDeformation)) {
           // Ke = elements[e]->stiffnessMatrixLD();
        }else{
            Ke = elements[e]->stiffnessMatrix();
        }
        const std::vector<int>& map = element_dofs[e];

        for (int i = 0; i < static_cast<int>(map.size()); ++i) {
            for (int j = 0; j < static_cast<int>(map.size()); ++j) {
                K_(map[i], map[j]) += Ke(i, j);
            }
        }
    }
    this->K = K_;
    return K_;
}

Eigen::MatrixXd Mesh::globalDumpingMatrix_RayleighMethode(double f1, double f2, double zeta){
    Eigen::MatrixXd c = Eigen::MatrixXd::Zero(total_dofs, total_dofs);
    if (f1 <= 0.0 || f2 <= 0.0) {
        throw std::runtime_error("globalDampingMatrix(): frequencies must be positive.");
    }
    if (zeta < 0.0) {
        throw std::runtime_error("globalDampingMatrix(): damping ratio must be non-negative.");
    }

    const double w1 = 2.0 * M_PI * f1;
    const double w2 = 2.0 * M_PI * f2;

    const double b = 2.0 * zeta / (w1 + w2);
    const double a = b * w1 * w2;

    std::cout << "Rayleigh estimation: " << std::endl;
    std::cout << "For: F1 = " << f1 << "[Hz], F2 = " << f2 << "[Hz] and Z = " << zeta << "[-]"<< std::endl;
    std::cout << "Rayleigh coefficients: " << " a = " << a << "[-], b = "<< b << "[-]"<<std::endl;
    std::cout << std::endl;

    c = a * this->M + b * this->K;

    this->C = c;
    return c;
}

Eigen::MatrixXd Mesh::globalMassMatrix(){
    Eigen::MatrixXd M_ = Eigen::MatrixXd::Zero(total_dofs, total_dofs);

    for (std::size_t e = 0; e < elements.size(); ++e) {
        Eigen::MatrixXd Me = elements[e]->massMatrix();
        const std::vector<int>& map = element_dofs[e];

        for (int i = 0; i < static_cast<int>(map.size()); ++i) {
            for (int j = 0; j < static_cast<int>(map.size()); ++j) {
                M_(map[i], map[j]) += Me(i, j);
            }
        }
    }

    this->M = M_;
    return M_;
}

std::vector<unsigned int> Mesh::findNodesOnEdge(char edge, double tol = 1e-9) const {
    std::vector<unsigned int> edgeNodes;

    double xmin =  std::numeric_limits<double>::max();
    double xmax = -std::numeric_limits<double>::max();
    double ymin =  std::numeric_limits<double>::max();
    double ymax = -std::numeric_limits<double>::max();

    for (const auto& n : nodes) {
        xmin = std::min(xmin, n->x[0]);
        xmax = std::max(xmax, n->x[0]);
        ymin = std::min(ymin, n->x[1]);
        ymax = std::max(ymax, n->x[1]);
    }

    for (const auto& n : nodes) {
        double x = n->x[0];
        double y = n->x[1];

        bool onEdge = false;

        if (edge == 'L') onEdge = std::abs(x - xmin) < tol;
        if (edge == 'R') onEdge = std::abs(x - xmax) < tol;
        if (edge == 'B') onEdge = std::abs(y - ymin) < tol;
        if (edge == 'T') onEdge = std::abs(y - ymax) < tol;

        if (onEdge) {
            edgeNodes.push_back(n->ID);
        }
    }

    return edgeNodes;
}

int Mesh::findNodeByCoordinates(double x, double y, double tol) const {
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        const double xi = nodes[i]->x[0];
        const double yi = nodes[i]->x[1];

        if (std::abs(xi - x) < tol && std::abs(yi - y) < tol) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

// Large Deformation
Eigen::VectorXd Mesh::globalInternalForce(const Eigen::VectorXd& Uglob) const {
    if (Uglob.size() != static_cast<int>(total_dofs)) {
        throw std::runtime_error("Mesh::globalInternalForce(): global displacement vector has wrong size.");
    }

    Eigen::VectorXd Fint = Eigen::VectorXd::Zero(total_dofs);

    for (std::size_t e = 0; e < elements.size(); ++e) {
        const auto& edofs = element_dofs[e];
        Eigen::VectorXd ue(edofs.size());

        for (std::size_t i = 0; i < edofs.size(); ++i) {
            ue(i) = Uglob(edofs[i]);
        }

        Eigen::VectorXd fe = elements[e]->internalForce(ue);

        for (std::size_t i = 0; i < edofs.size(); ++i) {
            Fint(edofs[i]) += fe(i);
        }
    }

    return Fint;
}

Eigen::MatrixXd Mesh::globalTangentStiffness(const Eigen::VectorXd& Uglob) const {
    Eigen::MatrixXd KT = Eigen::MatrixXd::Zero(total_dofs, total_dofs);

    for (std::size_t e = 0; e < elements.size(); ++e) {
        const auto& edofs = element_dofs[e];
        Eigen::VectorXd ue(edofs.size());

        for (std::size_t i = 0; i < edofs.size(); ++i) {
            ue(i) = Uglob(edofs[i]);
        }

        Eigen::MatrixXd Ke = elements[e]->tangentStiffness(ue);

        if (Ke.rows() != static_cast<int>(edofs.size()) || Ke.cols() != static_cast<int>(edofs.size())) {
            throw std::runtime_error("Mesh::globalTangentStiffness(): element tangent matrix has wrong size.");
        }

        for (std::size_t i = 0; i < edofs.size(); ++i) {
            for (std::size_t j = 0; j < edofs.size(); ++j) {
                KT(edofs[i], edofs[j]) += Ke(i, j);
            }
        }
    }

    return KT;
}