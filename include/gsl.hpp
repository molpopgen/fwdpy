#ifndef FWDPY_GSL_HPP_
#define FWDPY_GSL_HPP_

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <memory>
namespace fwdpy
{
    namespace gsl
    {
        // Enable use of smart pointers for GSL types
        struct gsl_matrix_deleter
        {
            void
            operator()(gsl_matrix *l) noexcept
            {
                gsl_matrix_free(l);
            }
        };

        struct gsl_vector_deleter
        {
            void
            operator()(gsl_vector *l) noexcept
            {
                gsl_vector_free(l);
            }
        };

        using gsl_vector_ptr_t
            = std::unique_ptr<gsl_vector, gsl_vector_deleter>;
        using gsl_matrix_ptr_t
            = std::unique_ptr<gsl_matrix, gsl_matrix_deleter>;
    }
}

#endif
