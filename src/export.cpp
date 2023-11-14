#include "orlqp/symbolic/export.hpp"
#include <sstream>
#include <fstream>

namespace orlqp
{

    ExportSymbolicQPProblem::ExportSymbolicQPProblem(const SymbolicQPProblem::Ptr sym_qp,
                                                     const std::string &class_name,
                                                     const std::string &file_path)
        : ClassName(class_name), QP(sym_qp)
    {
        std::ofstream file_out;
        file_out.open(file_path);
        file_out.clear();

        file_out << generateHeader();
        file_out << generateHessian();
        file_out << generateGradient();
        file_out << generateLinearConstraint();
        file_out << generateLowerBound();
        file_out << generateUpperBound();
        file_out << generateEnder();

        file_out.close();
    }

    std::string ExportSymbolicQPProblem::generateHeader()
    {
        std::string CLASSNAME = ClassName;
        std::transform(CLASSNAME.begin(), CLASSNAME.end(), CLASSNAME.begin(), ::toupper);
        std::ostringstream oss;
        oss
        << "#ifndef ORLQP_EXPORTED_QP_" << CLASSNAME << "_HPP_\n"
        << "#define ORLQP_EXPORTED_QP_" << CLASSNAME << "_HPP_\n"
        << "\n"
        << "#include \"orlqp/orlqp.hpp\"\n"
        << "\n"
        << "namespace orlqp\n"
        << "{\n"
        << "class " << ClassName << "\n"
        << "{\n"
        << "public:\n"
        << "using Ptr = std::shared_ptr<" << ClassName << ">;\n"
        << "\n"
        << ClassName << "() {}\n"
        << "\n";
    }

    std::string ExportSymbolicQPProblem::generateHessian()
    {
        /* TODO */
    }

    std::string ExportSymbolicQPProblem::generateGradient()
    {
        /* TODO */
    }

    std::string ExportSymbolicQPProblem::generateLinearConstraint()
    {
        /* TODO */
    }

    std::string ExportSymbolicQPProblem::generateLowerBound()
    {
        /* TODO */
    }

    std::string ExportSymbolicQPProblem::generateUpperBound()
    {
        /* TODO */
    }

    std::string ExportSymbolicQPProblem::generateEnder()
    {
        /* TODO */
    }

}