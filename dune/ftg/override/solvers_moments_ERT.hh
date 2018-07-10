// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MODELLING_SOLVERS_MOMENTS_ERT_HH
#define DUNE_MODELLING_SOLVERS_MOMENTS_ERT_HH

#include<dune/pdelab/finiteelement/localbasiscache.hh>
#include<dune/ftg/override/linearproblem.hh>
#include<dune/ftg/override/onestepparameter.hh>
#include<dune/pdelab/gridoperator/onestep.hh>
#include<dune/modelling/declarations.hh>

namespace Dune {
  namespace Modelling {

    /**
     * @brief Solver for stationary linear PDEs
     */
    template<typename Traits, typename ModelType, typename DirectionType>
      class StationarySolverMomentsERT_CG_AMG_SSOR_reuse_matrix
      {
        private:

          using RF = typename Traits::GridTraits::RangeField;
          using DF = typename Traits::GridTraits::DomainField;

          using DiscType   = typename EquationTraits<Traits,ModelType,DirectionType>::DiscretizationType;
          using GFS        = typename EquationTraits<Traits,ModelType,DirectionType>::GridFunctionSpace;
          using GridVector = typename EquationTraits<Traits,ModelType,DirectionType>::GridVector;
          using C          = typename GFS::template ConstraintsContainer<RF>::Type;

          using LOP = SpatialOperator<Traits,ModelType,DiscType,DirectionType>;

          using MBE = Dune::PDELab::istl::BCRSMatrixBackend<>;
          using GO  = Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,DF,RF,RF,C,C>;
          using M   = typename GO::template MatrixContainer<RF>::Type;

          using P0FEM = Dune::PDELab::P0LocalFiniteElementMap<DF,RF,Traits::GridTraits::dim>;
          using P0VBE = Dune::PDELab::istl::VectorBackend<Dune::PDELab::istl::Blocking::fixed,1>;
          using P0CON = Dune::PDELab::P0ParallelConstraints;
          using P0GFS = Dune::PDELab::GridFunctionSpace<typename Traits::GridTraits::GridView,P0FEM,P0CON,P0VBE>;
          using P0C   = typename P0GFS::template ConstraintsContainer<RF>::Type;
          using LS    = Dune::PDELab::ISTLBackend_CG_AMG_SSOR<GO>; // ISTLBackend_BCGS_AMG_ILU0/ISTLBackend_CG_AMG_SSOR

          const Traits& traits;
          const ModelParameters<Traits,ModelType>&  parameters;
          const EquationTraits<Traits,ModelType,DirectionType>& equationTraits;
          LOP   lop;

          MBE   mbe;
          GO    go;
          M     m;

          P0FEM p0fem;
          P0GFS p0gfs;
          P0C   p0cg;
          //LS    ls;

          std::shared_ptr<LS> ls_ptr;
          using ERTMatrixContainer = typename Traits::ERTMatrixContainer;
          std::shared_ptr<ERTMatrixContainer> ertMatrixContainer;

        public:

          /**
           * @brief Constructor
           */
          StationarySolverMomentsERT_CG_AMG_SSOR_reuse_matrix(
              const Traits& traits_,
              const EquationTraits<Traits,ModelType,DirectionType>& equationTraits_,
              const C& cg,
              const ModelParameters<Traits,ModelType>& parameters_,
              RF Tend,
              std::shared_ptr<typename Traits::ERTMatrixContainer>& ertMatrixContainer_
              )
            : traits(traits_),
            parameters(parameters_),
            equationTraits(equationTraits_), lop(traits_,parameters_), mbe(9),
            go(equationTraits.gfs(),cg,equationTraits.gfs(),cg,lop,mbe), m(go),
            p0fem(Dune::GeometryType(Dune::GeometryType::cube,Traits::GridTraits::dim)),
            p0gfs(equationTraits.gfs().gridView(),p0fem),
            ertMatrixContainer(ertMatrixContainer_)
        {
          if (parameters.model_number==0 && parameters.k == 0)
          {
            std::shared_ptr<LS> ls_ptr(new LS(equationTraits.gfs(),5000,1,true,true)); // max_iter, verbose, reuse, superLU
            (*ertMatrixContainer).set_ls_moments_ERT(ls_ptr);
          }
        }

          /**
           * @brief Compute solution (since problem is stationary)
           */
          unsigned int step(RF time, RF timestep, GridVector& oldSolution, GridVector& solution)
          {
            GridVector residual(equationTraits.gfs(),0.);
            lop.setTime(time+timestep);

            go.residual(solution, residual);
            // only re-evaluate matrix for the first ERT model
            //if (parameters.model_number==0)
            //{
            //  m = 0.;
            //  go.jacobian(solution,m);
            //  (*ertMatrixContainer).set_matrix_moments_ERT(m);
            //} else {
            m = (*ertMatrixContainer).read_matrix_ERT();
            //}
            GridVector z(equationTraits.gfs(),0.);
            (*(*ertMatrixContainer).read_ls_moments_ERT()).apply(m,z,residual,1e-12);   // check if pointer is valid!
            solution -= z;

            return 1;
          }
      };

    /**
     * @brief Dummy implementation for transient solver
     */
    template<typename Traits, typename ModelType, typename DirectionType>
      class TransientSolverMomentsERT
      {
        private:

          using GridTraits = typename Traits::GridTraits;
          using RF     = typename GridTraits::RangeField;
          using GFS        = typename EquationTraits<Traits,ModelType,DirectionType>::GridFunctionSpace;
          using GridVector = typename EquationTraits<Traits,ModelType,DirectionType>::GridVector;
          using C          = typename GFS::template ConstraintsContainer<RF>::Type;

          using ERTMatrixContainer = typename Traits::ERTMatrixContainer;
          std::shared_ptr<ERTMatrixContainer> ertMatrixContainer;
        public:

          /**
           * @brief Constructor
           */
          TransientSolverMomentsERT(
              const Traits& traits,
              const EquationTraits<Traits,ModelType,DirectionType>& equationTraits_,
              const C& cg,
              ModelParameters<Traits,ModelType>& parameters_,
              RF Tend,
              std::shared_ptr<typename Traits::ERTMatrixContainer>& ertMatrixContainer_)
        {}

          /**
           * @brief Compute a time step
           */
          unsigned int step(RF timeFrom, RF timestep, GridVector& oldSolution, GridVector& newSolution)
          {
            DUNE_THROW(Dune::Exception,"transient solver for ERT moments model not implemented!");
          }
      };
  }
}

#endif // DUNE_MODELLING_SOLVERS_MOMENTS_ERT_HH
