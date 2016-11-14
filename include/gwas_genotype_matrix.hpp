#ifndef FWDPY_GWAS_GENOTYPE_MATRIX_HPP
#define FWDPY_GWAS_GENOTYPE_MATRIX_HPP

#include "types.hpp"
#include <vector>
#include <cstdint>
namespace fwdpy
{
    namespace gwas
    {
        struct genotype_matrix
        /*!
          Intended to be coerced to 2D numpy matrix
          */
        {
            std::vector<double> G, E, neutral, causative, npos, cpos,
                q_neutral, q_causative;
            std::size_t N, n_neutral, n_causative;
            genotype_matrix()
                : G{}, E{}, neutral{}, causative{}, npos{}, cpos{},
                  q_neutral{}, q_causative{}, N{ 0 }, n_neutral{ 0 },
                  n_causative{ 0 }
            {
            }
            genotype_matrix(const genotype_matrix &) = default;
        };

        genotype_matrix make_geno_matrix(const singlepop_t *pop,
                                         std::vector<std::size_t> &individuals,
                                         bool maf = true);
        genotype_matrix make_geno_matrix(const singlepop_t *pop,
                                         bool maf = true);
    }
}

#endif
