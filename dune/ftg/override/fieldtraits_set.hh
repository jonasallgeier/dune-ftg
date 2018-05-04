// -*- tab-width: 2; indent-tabs-mode: nil -*-
#ifndef DUNE_RANDOMFIELD_FIELDTRAITS_HH
#define	DUNE_RANDOMFIELD_FIELDTRAITS_HH

#include<array>
#include<vector>

#include <fftw3.h>
#include <fftw3-mpi.h>

#include<dune/common/parametertreeparser.hh>

namespace Dune {
  namespace RandomField {

    /**
     * @brief Default load balance strategy, taken from dune-grid to avoid hard dependency
     */
    template<long unsigned int dim>
      class DefaultLoadBalance
      {
        public:

          /** @brief Distribute a structured grid across a set of processors
           *
           * @param [in] size Number of elements in each coordinate direction, for the entire grid
           * @param [in] P    Number of processors
           */
          void loadbalance (const std::array<int, dim>& size, int P, std::array<int,dim>& dims) const
          {
            double opt = 1e100;
            std::array<int,dim> trydims;

            optimize_dims(dim-1,size,P,dims,trydims,opt);
          }

        private:

          void optimize_dims (
              int i,
              const std::array<int,dim>& size,
              int P, std::array<int,dim>& dims,
              std::array<int,dim>& trydims,
              double& opt
              ) const
          {
            if (i > 0) // test all subdivisions recursively
            {
              for (unsigned int k = 1; k <= (unsigned int)P; k++)
                if (P%k == 0)
                {
                  // P divisible by k
                  trydims[i] = k;
                  optimize_dims(i-1,size,P/k,dims,trydims,opt);
                }
            }
            else
            {
              // found a possible combination
              trydims[0] = P;

              // check for optimality
              double m = -1.;

              for (unsigned int k = 0; k < dim; k++)
              {
                double mm = ((double)size[k])/((double)trydims[k]);
                if (fmod((double)size[k],(double)trydims[k]) > 0.0001)
                  mm *= 3;
                if ( mm > m )
                  m = mm;
              }
              if (m < opt)
              {
                opt = m;
                dims = trydims;
              }
            }
          }
      };

    template<typename Traits> class TrendPart;
    template<typename Traits> class TrendComponent;
    template<typename Traits> class StochasticPart;
    template<typename Traits> class RandomFieldMatrix;
    template<typename GridTraits, bool storeInvMat, bool storeInvRoot> class RandomField;

    /**
     * @brief Traits for the RandomField class
     */
    template<typename GridTraits, bool storeInvMat, bool storeInvRoot>
      class RandomFieldTraits
      {
        typedef RandomFieldTraits<GridTraits,storeInvMat,storeInvRoot> ThisType;

        public:

        enum {dim = GridTraits::dim};

        typedef typename GridTraits::RangeField  RF;
        typedef typename GridTraits::DomainField DomainField;
        typedef typename GridTraits::Domain      DomainType;
        typedef typename GridTraits::Scalar      RangeType;

#if HAVE_DUNE_PDELAB
        // allows treating a RandomField as a PDELab function
        typedef typename Dune::YaspGrid<dim>::LeafGridView GridViewType;
        enum {dimRange  = 1};
        enum {dimDomain = dim};
#endif // HAVE_DUNE_PDELAB

        private:

        friend class TrendPart<ThisType>;
        friend class TrendComponent<ThisType>;
        friend class StochasticPart<ThisType>;
        friend class RandomFieldMatrix<ThisType>;
        friend class RandomField<GridTraits,storeInvMat,storeInvRoot>;

        // MPI constants
        int rank, commSize;

        std::array<int,dim> procPerDim;

        const Dune::ParameterTree& config;
        const MPI_Comm             comm;

        const std::array<RF,dim> extensions;
        unsigned int             level;
        std::array<RF,dim>       meshsize;
        RF                       cellVolume;

        const RF           variance;
        const std::string  covariance;
        const bool         periodic;
        const bool         verbose;
        const unsigned int cgIterations;

        ptrdiff_t allocLocal, localN0, local0Start;

        // factor used in domain embedding
        unsigned int embeddingFactor;

        // properties of random field
        std::array<unsigned int,dim> cells;
        unsigned int                 domainSize;
        std::array<unsigned int,dim> localCells;
        std::array<unsigned int,dim> localOffset;
        unsigned int                 localDomainSize;

        //std::set<RF>               set_x;
        //std::set<RF>               set_y;
        //std::set<RF>               set_z;
        std::array<std::set<RF>,dim> gridCoords;

        // properties on extended domain
        std::array<unsigned int,dim> extendedCells;
        unsigned int                 extendedDomainSize;
        std::array<unsigned int,dim> localExtendedCells;
        std::array<unsigned int,dim> localExtendedOffset;
        unsigned int                 localExtendedDomainSize;

        mutable std::array<unsigned int,dim> globalIndices;
        mutable std::array<unsigned int,dim> localIndices;

        public:

        template<typename LoadBalance>
          RandomFieldTraits(
              const Dune::ParameterTree& config_,
              const LoadBalance& loadBalance,
              const MPI_Comm comm_
              )
          : config(config_), comm(comm_),
          extensions     (config.get<std::array<RF,dim> >          ("grid.extensions")),
          variance       (config.get<RF>                           ("stochastic.variance")),
          covariance     (config.get<std::string>                  ("stochastic.covariance")),
          periodic       (config.get<bool>                         ("randomField.periodic",false)),
          verbose        (config.get<bool>                         ("randomField.verbose",false)),
          cgIterations   (config.get<unsigned int>                 ("randomField.cgIterations",100)),
          embeddingFactor(config.get<unsigned int>                 ("randomField.embeddingFactor",2)),
          cells          (config.get<std::array<unsigned int,dim> >("grid.cells"))
        {

          std::vector<RF> test1 = config.get<std::vector<RF> >("grid.vector_x");
          std::vector<RF> test2 = config.get<std::vector<RF> >("grid.vector_y");
          std::vector<RF> test3 = config.get<std::vector<RF> >("grid.vector_z");

          std::set<RF> interm1(test1.begin(), test1.end());
          std::set<RF> interm2(test2.begin(), test2.end());
          std::set<RF> interm3(test3.begin(), test3.end());
          gridCoords[0] = interm1;
          gridCoords[1] = interm2;
          gridCoords[2] = interm3;
          //std::vector<RF> vector_y = config.get<std::vector<RF> >("grid.vector_y");
          //std::vector<RF> vector_z = config.get<std::vector<RF> >("grid.vector_z");
          

          MPI_Comm_rank(comm,&rank);
          MPI_Comm_size(comm,&commSize);

          // dune-grid load balancers want int as data type
          std::array<int,dim> intCells;
          for (unsigned int i = 0; i < dim; i++)
            intCells[i] = cells[i];
          loadBalance.loadbalance(intCells,commSize,procPerDim);

          level = 0;

          if (periodic && embeddingFactor != 1)
          {
            if (verbose && rank == 0)
              std::cout << "periodic boundary conditions are synonymous with embeddingFactor == 1,"
                << " enforcing consistency" << std::endl;
            embeddingFactor = 1;
          }

          fftw_mpi_init();
          update();
        }

        /**
         * @brief Compute constants after construction or refinement
         */
        void update()
        {
          // ensures that FFTW can divide data equally between processes
          if (cells[dim-1] % commSize != 0)
            DUNE_THROW(Dune::Exception,"number of cells in last dimension has to be multiple of numProc");
          if (dim == 1 && cells[0] % (commSize*commSize) != 0)
            DUNE_THROW(Dune::Exception,"in 1D, number of cells has to be multiple of numProc^2");

          for (unsigned int i = 0; i < dim; i++)
          {
            meshsize[i]      = extensions[i] / cells[i];
            extendedCells[i] = embeddingFactor*cells[i];
          }

          getFFTData(allocLocal, localN0, local0Start);

          for (unsigned int i = 0; i < dim - 1; i++)
          {
            localExtendedCells [i] = extendedCells[i];
            localExtendedOffset[i] = 0;
            localCells [i] = cells[i];
            localOffset[i] = 0;
          }
          localExtendedCells [dim-1] = localN0;
          localExtendedOffset[dim-1] = local0Start;
          localCells [dim-1] = localN0/embeddingFactor;
          localOffset[dim-1] = local0Start/embeddingFactor;

          domainSize              = 1;
          extendedDomainSize      = 1;
          localDomainSize         = 1;
          localExtendedDomainSize = 1;
          cellVolume              = 1.;
          for (unsigned int i = 0; i < dim; i++)
          {
            domainSize              *= cells[i];
            extendedDomainSize      *= extendedCells[i];
            localDomainSize         *= localCells[i];
            localExtendedDomainSize *= localExtendedCells[i];
            cellVolume              *= meshsize[i];
          }

          if (verbose && rank == 0)
          {
            std::cout << "RandomField size:        " << localDomainSize << std::endl;
            std::cout << "RandomField cells:       ";
            for (unsigned int i = 0; i < dim; i++)
            {
              std::cout << cells[i] << " ";
            }
            std::cout << std::endl;
            std::cout << "RandomField local cells: ";
            for (unsigned int i = 0; i < dim; i++)
            {
              std::cout << localCells[i] << " ";
            }
            std::cout << std::endl;
            std::cout << "RandomField cell volume: " << cellVolume << std::endl;
          }
        }

        /**
         * @brief Request global refinement of the data structure
         */
        void refine()
        {
          for (unsigned int i = 0; i < dim; ++i)
            cells[i] *= 2;

          level++;

          update();
        }

        /**
         * @brief Get the domain decomposition data of the Fourier transform
         */
        template<typename T>
          void getFFTData(T& allocLocal, T& localN0, T& local0Start) const
          {
            if (dim == 3)
            {
              ptrdiff_t n[] = {(ptrdiff_t)extendedCells[0],
                (ptrdiff_t)extendedCells[1],(ptrdiff_t)extendedCells[2]};
              allocLocal = fftw_mpi_local_size_3d(n[2] , n[1], n[0], comm, &localN0, &local0Start);
            }
            else if (dim == 2)
            {
              ptrdiff_t n[] = {(ptrdiff_t)extendedCells[0],(ptrdiff_t)extendedCells[1]};
              allocLocal = fftw_mpi_local_size_2d(n[1], n[0], comm, &localN0, &local0Start);
            }
            else if (dim == 1)
            {
              ptrdiff_t n[] = {(ptrdiff_t)extendedCells[0]};
              ptrdiff_t localN02, local0Start2;
              allocLocal = fftw_mpi_local_size_1d(n[0], comm, FFTW_FORWARD, FFTW_ESTIMATE,
                  &localN0, &local0Start, &localN02, &local0Start2);
              if (localN0 != localN02 || local0Start != local0Start2)
                DUNE_THROW(Dune::Exception,"1d size / offset results don't match");
            }
            else
              DUNE_THROW(Dune::Exception,"dimension of field has to be 1, 2 or 3");
          }

        /**
         * @brief Convert an index tuple into a one dimensional encoding
         */
        unsigned int indicesToIndex(
            const std::array<unsigned int,dim>& indices,
            const std::array<unsigned int,dim>& bound
            ) const
        {
          if (dim == 3)
          {
            return indices[0] + bound[0] * (indices[1] + bound[1]*indices[2]);
          }
          else if (dim == 2)
          {
            return indices[1] * bound[0] + indices[0];
          }
          else if (dim == 1)
          {
            return indices[0];
          }
          else
            DUNE_THROW(Dune::Exception,"dimension of field has to be 1, 2 or 3");
        }

        /**
         * @brief Convert a one dimensional encoding into the original index tuple
         */
        void indexToIndices(
            const unsigned int index,
            std::array<unsigned int,dim>& indices,
            const std::array<unsigned int,dim>& bound
            ) const
        {
          if (dim == 3)
          {
            indices[0] = index % bound[0];
            indices[1] = (index / bound[0]) % bound[1];
            indices[2] = (index / bound[0]) / bound[1];
          }
          else if (dim == 2)
          {
            indices[0] = index % bound[0];
            indices[1] = index / bound[0];
          }
          else if (dim == 1)
          {
            indices[0] = index;
          }
          else
            DUNE_THROW(Dune::Exception,"dimension of field has to be 1, 2 or 3");
        }

        /**
         * @brief Convert spatial coordinates into the corresponding integer indices
         */
        void coordsToIndices(
            const DomainType& location,
            std::array<unsigned int,dim>& localIndices,
            const std::array<unsigned int,dim>& offset
            ) const
        {
          for (unsigned int i = 0; i < dim; i++)
          {
            /*if (i==0)
            {
              for (unsigned int j = 1; j!=vector_x.size();j++)
              {
                if ((location[i] < vector_x[j]) & (location[i] >=vector_x[j-1]) ) 
                {
                  globalIndices[i] = j-1;
                  break;
                };
                globalIndices[i] = vector_x.size()-2;
              }
            }
            else if (i==1)
            {
              for (unsigned int j = 1; j!=vector_y.size();j++)
              {
                if ((location[i] < vector_y[j]) & (location[i] >=vector_y[j-1]) ) 
                {
                  globalIndices[i] = j-1;
                  break;
                };
                globalIndices[i] = vector_y.size()-2;
              }
            }
            else if (i==2)
            {
              for (unsigned int j = 1; j!=vector_z.size();j++)
              {
                if ((location[i] < vector_z[j]) & (location[i] >=vector_z[j-1]) ) 
                {
                  globalIndices[i] = j-1;
                  break;
                };
                globalIndices[i] = vector_z.size()-2; // subtract 2, because: there is one a cell fewer than mesh vector indices and at the very boundary we use the value of the last cell!
              }
            } */
            
            
            std::set<RF> coords = gridCoords[i];
            typename std::set<RF>::iterator it;
            
            it = std::upper_bound (coords.begin(), coords.end(), location[i]);
            if (it != coords.end())
            {
              //std::cout << *std::prev(it);
              //std::cout << " "<< location[i] << " ";
              //std::cout << *it  << "-> cell: ";
              //std::cout << std::distance(coords.begin(), it)-1 << " || ";
              globalIndices[i] = std::distance(coords.begin(), it)-1;
            }
            else
            {
              //std::cout << *std::prev(it);
              //std::cout << " "<< location[i] << " ";
              //std::cout << *(coords.rbegin())  << "-> cell: ";
              //std::cout << coords.size()-2 << "(!) ";
              //std::cout << std::endl;
              globalIndices[i] = coords.size()-2;
            }
          //globalIndices[i] = (unsigned int) (location[i] * cells[i] / extensions[i]);
            localIndices[i]  = globalIndices[i] - offset[i];
          }
        }

        /**
         * @brief Convert integer indices into corresponding spatial coordinates
         */
        void indicesToCoords(
            const std::array<unsigned int,dim>& localIndices,
            const std::array<unsigned int,dim>& offset,
            DomainType& location
            ) const
        {
          for (unsigned int i = 0; i < dim; i++)
          {
            globalIndices[i] = localIndices[i] + offset[i];
            location[i]      = (globalIndices[i] * extensions[i] + 0.5) / cells[i];
          }
        }

      };

  }
}

#endif // DUNE_RANDOMFIELD_FIELDTRAITS_HH
