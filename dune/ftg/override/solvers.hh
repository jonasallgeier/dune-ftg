// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MODELLING_SOLVERS_HH
#define DUNE_MODELLING_SOLVERS_HH

#include<dune/pdelab/finiteelement/localbasiscache.hh>
#include<dune/pdelab/gridoperator/onestep.hh>
#include<dune/pdelab/stationary/linearproblem.hh>
#include<dune/pdelab/instationary/onestepparameter.hh>

#include<dune/modelling/declarations.hh>

namespace Dune {
  namespace Modelling {

    /**
     * @brief Solver for stationary linear PDEs
     */
    template<typename Traits, typename ModelType, typename DirectionType>
      class StationaryLinearSolver_CG_AMG_SSOR_reuse_matrix
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
          StationaryLinearSolver_CG_AMG_SSOR_reuse_matrix(
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
          if (parameters.model_number==0)
          {
            //ls(equationTraits.gfs(),5000,1,false,true); // max_iter, verbose, reuse, superLU
            
            std::shared_ptr<LS> ls_ptr(new LS(equationTraits.gfs(),5000,1,true,true));
            (*ertMatrixContainer).set_ls(ls_ptr);
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
            if (parameters.model_number==0)
            {
              m = 0.;
              go.jacobian(solution,m);
              (*ertMatrixContainer).set_matrix(m);
            } else {
              m = (*ertMatrixContainer).read_matrix();
            }
            GridVector z(equationTraits.gfs(),0.);
            (*(*ertMatrixContainer).read_ls()).apply(m,z,residual,1e-6);   // check if pointer is valid!
            solution -= z;

            return 1;
          }
      };

    /**
     * @brief Solver for stationary linear PDEs
     */
    template<typename Traits, typename ModelType, typename DirectionType>
      class StationaryLinearSolver_CG_AMG_SSOR
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

          const EquationTraits<Traits,ModelType,DirectionType>& equationTraits;
          LOP   lop;

          MBE   mbe;
          GO    go;
          M     m;

          P0FEM p0fem;
          P0GFS p0gfs;
          P0C   p0cg;
          LS    ls;
          
          using ERTMatrixContainer = typename Traits::ERTMatrixContainer;
          std::shared_ptr<ERTMatrixContainer> ertMatrixContainer;

        public:

          /**
           * @brief Constructor
           */
          StationaryLinearSolver_CG_AMG_SSOR(
              const Traits& traits,
              const EquationTraits<Traits,ModelType,DirectionType>& equationTraits_,
              const C& cg,
              const ModelParameters<Traits,ModelType>& parameters,
              RF Tend,
              std::shared_ptr<typename Traits::ERTMatrixContainer>& ertMatrixContainer_
              )
            : equationTraits(equationTraits_), lop(traits,parameters), mbe(9),
            go(equationTraits.gfs(),cg,equationTraits.gfs(),cg,lop,mbe), m(go),
            p0fem(Dune::GeometryType(Dune::GeometryType::cube,Traits::GridTraits::dim)),
            p0gfs(equationTraits.gfs().gridView(),p0fem),
            ls(equationTraits.gfs(),5000,1,false,true), // max_iter, verbose, reuse, superLU
            ertMatrixContainer(ertMatrixContainer_)
        {}

          /**
           * @brief Compute solution (since problem is stationary)
           */
          unsigned int step(RF time, RF timestep, GridVector& oldSolution, GridVector& solution)
          {
            GridVector residual(equationTraits.gfs(),0.);
            lop.setTime(time+timestep);

            go.residual(solution, residual);
            m = 0.;
            go.jacobian(solution,m);
            GridVector z(equationTraits.gfs(),0.);
            ls.apply(m,z,residual,1e-10);
            solution -= z;

            return 1;
          }
      };

    /**
     * @brief Solver for stationary linear PDEs
     */
    template<typename Traits, typename ModelType, typename DirectionType>
      class StationaryLinearSolver_BCGS_AMG_ILU0
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
          using LS    = Dune::PDELab::ISTLBackend_BCGS_AMG_ILU0<GO>; // ISTLBackend_BCGS_AMG_ILU0/ISTLBackend_CG_AMG_SSOR

          const EquationTraits<Traits,ModelType,DirectionType>& equationTraits;
          LOP   lop;

          MBE   mbe;
          GO    go;
          M     m;

          P0FEM p0fem;
          P0GFS p0gfs;
          P0C   p0cg;
          LS    ls;

          using ERTMatrixContainer = typename Traits::ERTMatrixContainer;
          std::shared_ptr<ERTMatrixContainer> ertMatrixContainer;

        public:

          /**
           * @brief Constructor
           */
          StationaryLinearSolver_BCGS_AMG_ILU0(
              const Traits& traits,
              const EquationTraits<Traits,ModelType,DirectionType>& equationTraits_,
              const C& cg,
              const ModelParameters<Traits,ModelType>& parameters,
              RF Tend,
              std::shared_ptr<typename Traits::ERTMatrixContainer>& ertMatrixContainer_)
            : equationTraits(equationTraits_), lop(traits,parameters), mbe(9),
            go(equationTraits.gfs(),cg,equationTraits.gfs(),cg,lop,mbe), m(go),
            p0fem(Dune::GeometryType(Dune::GeometryType::cube,Traits::GridTraits::dim)),
            p0gfs(equationTraits.gfs().gridView(),p0fem),
            ls(equationTraits.gfs(),5000,1,false,true), // max_iter, verbose, reuse, superLU
            ertMatrixContainer(ertMatrixContainer_)
        {}

          /**
           * @brief Compute solution (since problem is stationary)
           */
          unsigned int step(RF time, RF timestep, GridVector& oldSolution, GridVector& solution)
          {
            GridVector residual(equationTraits.gfs(),0.);
            lop.setTime(time+timestep);

            go.residual(solution, residual);
            m = 0.;
            go.jacobian(solution,m);
            GridVector z(equationTraits.gfs(),0.);
            ls.apply(m,z,residual,1e-10);
            solution -= z;

            return 1;
          }
      };

    /**
     * @brief Solver for linear PDEs using implicit timestepping scheme
     */
    template<typename Traits, typename ModelType, typename DirectionType>
      class ImplicitLinearSolver
      {
        private:

          using GridTraits = typename Traits::GridTraits;

          using GV     = typename GridTraits::GridView;
          using RF     = typename GridTraits::RangeField;
          using Scalar = typename GridTraits::Scalar;
          using DF     = typename GridTraits::DomainField;
          using Domain = typename GridTraits::Domain;

          using DiscType   = typename EquationTraits<Traits,ModelType,DirectionType>::DiscretizationType;
          using GFS        = typename EquationTraits<Traits,ModelType,DirectionType>::GridFunctionSpace;
          using GridVector = typename EquationTraits<Traits,ModelType,DirectionType>::GridVector;
          using C          = typename GFS::template ConstraintsContainer<RF>::Type;

          using LOP  = SpatialOperator <Traits,ModelType,DiscType,DirectionType>;
          using TLOP = TemporalOperator<Traits,ModelType,DiscType,DirectionType>;

          using Method = typename EquationTraits<Traits,ModelType,DirectionType>::OneStepScheme;

          using MBE = Dune::PDELab::istl::BCRSMatrixBackend<>;
          using GO0 = Dune::PDELab::GridOperator<GFS,GFS, LOP,MBE,DF,RF,RF,C,C>;
          using GO1 = Dune::PDELab::GridOperator<GFS,GFS,TLOP,MBE,DF,RF,RF,C,C>;
          using IGO = Dune::PDELab::OneStepGridOperator<GO0,GO1>;

          using P0FEM = Dune::PDELab::P0LocalFiniteElementMap<DF,RF,Traits::GridTraits::dim>;
          using P0VBE = Dune::PDELab::istl::VectorBackend<Dune::PDELab::istl::Blocking::fixed,1>;
          using P0CON = Dune::PDELab::P0ParallelConstraints;
          using P0GFS = Dune::PDELab::GridFunctionSpace<typename Traits::GridTraits::GridView,P0FEM,P0CON,P0VBE>;
          using P0C   = typename P0GFS::template ConstraintsContainer<RF>::Type;
          using LS    = Dune::PDELab::ISTLBackend_CG_AMG_SSOR<IGO>;

          using PDESolver = Dune::PDELab::StationaryLinearProblemSolver<IGO,LS,GridVector>;
          using OSM       = Dune::PDELab::OneStepMethod<RF,IGO,PDESolver,GridVector,GridVector>;

          const Dune::ParameterTree&                            config;
          const EquationTraits<Traits,ModelType,DirectionType>& equationTraits;
          ModelParameters<Traits,ModelType>&                    parameters;

          const RF minStep;
          RF smallStep;

          LOP lop;
          TLOP tlop;
          const Method method;
          const Dune::PDELab::ImplicitEulerParameter<RF> eulerMethod;

          MBE       mbe;
          GO0       go0;
          GO1       go1;
          IGO       igo;

          P0FEM     p0fem;
          P0GFS     p0gfs;
          P0C       p0cg;
          LS        ls;
          PDESolver solver;
          OSM       osm;

          using ERTMatrixContainer = typename Traits::ERTMatrixContainer;
          std::shared_ptr<ERTMatrixContainer> ertMatrixContainer;
        public:

          /**
           * @brief Constructor
           */
          ImplicitLinearSolver(
              const Traits& traits,
              const EquationTraits<Traits,ModelType,DirectionType>& equationTraits_,
              const C& cg,
              ModelParameters<Traits,ModelType>& parameters_,
              RF Tend,
              std::shared_ptr<typename Traits::ERTMatrixContainer>& ertMatrixContainer_)
            : config(traits.config()), equationTraits(equationTraits_), parameters(parameters_),
            minStep(config.get<RF>("time.step_transport")),
            smallStep(0.),
            lop(traits,parameters), tlop(traits,parameters), mbe(9),
            go0(equationTraits.gfs(),cg,equationTraits.gfs(),cg,lop,mbe),
            go1(equationTraits.gfs(),cg,equationTraits.gfs(),cg,tlop,mbe), igo(go0,go1),
            p0fem(Dune::GeometryType(Dune::GeometryType::cube,Traits::GridTraits::dim)),
            p0gfs(equationTraits.gfs().gridView(),p0fem),
            //ls(igo,cg,p0gfs,p0cg,config.sub("amg")), solver(igo,ls,1e-6),
            ls(equationTraits.gfs()), solver(igo,ls,1e-6),
            osm(method,igo,solver),
            ertMatrixContainer(ertMatrixContainer_)
        {
          osm.setVerbosityLevel(config.get<int>("solver.verbosity",3));
        }

          /**
           * @brief Compute a time step
           */
          unsigned int step(RF timeFrom, RF timestep, GridVector& oldSolution, GridVector& newSolution)
          {
            osm.apply(timeFrom,timestep,oldSolution,newSolution);
            return 1;
          }

      };

   /**
     * @brief Solver for linear PDEs using implicit timestepping scheme
     */
    template<typename Traits, typename ModelType, typename DirectionType>
      class CrankNicolsonSolver
      {
        private:

          using GridTraits = typename Traits::GridTraits;

          using GV     = typename GridTraits::GridView;
          using RF     = typename GridTraits::RangeField;
          using Scalar = typename GridTraits::Scalar;
          using DF     = typename GridTraits::DomainField;
          using Domain = typename GridTraits::Domain;

          using DiscType   = typename EquationTraits<Traits,ModelType,DirectionType>::DiscretizationType;
          using GFS        = typename EquationTraits<Traits,ModelType,DirectionType>::GridFunctionSpace;
          using GridVector = typename EquationTraits<Traits,ModelType,DirectionType>::GridVector;
          using C          = typename GFS::template ConstraintsContainer<RF>::Type;

          using LOP  = SpatialOperator <Traits,ModelType,DiscType,DirectionType>;
          using TLOP = TemporalOperator<Traits,ModelType,DiscType,DirectionType>;

          using Method = Dune::PDELab::OneStepThetaParameter<RF>;

          using MBE = Dune::PDELab::istl::BCRSMatrixBackend<>;
          using GO0 = Dune::PDELab::GridOperator<GFS,GFS, LOP,MBE,DF,RF,RF,C,C>;
          using GO1 = Dune::PDELab::GridOperator<GFS,GFS,TLOP,MBE,DF,RF,RF,C,C>;
          using IGO = Dune::PDELab::OneStepGridOperator<GO0,GO1>;

          using P0FEM = Dune::PDELab::P0LocalFiniteElementMap<DF,RF,Traits::GridTraits::dim>;
          using P0VBE = Dune::PDELab::istl::VectorBackend<Dune::PDELab::istl::Blocking::fixed,1>;
          using P0CON = Dune::PDELab::P0ParallelConstraints;
          using P0GFS = Dune::PDELab::GridFunctionSpace<typename Traits::GridTraits::GridView,P0FEM,P0CON,P0VBE>;
          using P0C   = typename P0GFS::template ConstraintsContainer<RF>::Type;
          using LS    = Dune::PDELab::ISTLBackend_CG_AMG_SSOR<IGO>;

          using PDESolver = Dune::PDELab::StationaryLinearProblemSolver<IGO,LS,GridVector>;
          using OSM       = Dune::PDELab::OneStepMethod<RF,IGO,PDESolver,GridVector,GridVector>;

          const Dune::ParameterTree&                            config;
          const EquationTraits<Traits,ModelType,DirectionType>& equationTraits;
          ModelParameters<Traits,ModelType>&                    parameters;

          const RF minStep;
          RF smallStep;

          const Method method;
          LOP lop;
          TLOP tlop;
          //const Dune::PDELab::OneStepThetaParameter<RF> eulerMethod(0.5);

          MBE       mbe;
          GO0       go0;
          GO1       go1;
          IGO       igo;

          P0FEM     p0fem;
          P0GFS     p0gfs;
          P0C       p0cg;
          LS        ls;
          PDESolver solver;
          OSM       osm;

          using ERTMatrixContainer = typename Traits::ERTMatrixContainer;
          std::shared_ptr<ERTMatrixContainer> ertMatrixContainer;
        public:

          /**
           * @brief Constructor
           */
          CrankNicolsonSolver(
              const Traits& traits,
              const EquationTraits<Traits,ModelType,DirectionType>& equationTraits_,
              const C& cg,
              ModelParameters<Traits,ModelType>& parameters_,
              RF Tend,
              std::shared_ptr<typename Traits::ERTMatrixContainer>& ertMatrixContainer_)
            : config(traits.config()), equationTraits(equationTraits_), parameters(parameters_),
            minStep(config.get<RF>("time.step_transport")),
            smallStep(0.),
            method(0.5),
            lop(traits,parameters), tlop(traits,parameters), mbe(9),
            go0(equationTraits.gfs(),cg,equationTraits.gfs(),cg,lop,mbe),
            go1(equationTraits.gfs(),cg,equationTraits.gfs(),cg,tlop,mbe), igo(go0,go1),
            p0fem(Dune::GeometryType(Dune::GeometryType::cube,Traits::GridTraits::dim)),
            p0gfs(equationTraits.gfs().gridView(),p0fem),
            //ls(igo,cg,p0gfs,p0cg,config.sub("amg")), solver(igo,ls,1e-6),
            ls(equationTraits.gfs()), solver(igo,ls,1e-6),
            osm(method,igo,solver),
            ertMatrixContainer(ertMatrixContainer_)
        {
          osm.setVerbosityLevel(config.get<int>("solver.verbosity",3));
        }

          /**
           * @brief Compute a time step
           */
          unsigned int step(RF timeFrom, RF timestep, GridVector& oldSolution, GridVector& newSolution)
          {
            osm.apply(timeFrom,timestep,oldSolution,newSolution);
            return 1;
          }

      };

    /**
     *  @brief Time step limitation due to CFL ignoring sign of timestep
     */
    template<class RF, class IGO> 
      class ScalingTimeController : public Dune::PDELab::TimeControllerInterface<RF>
    {
      public:

        typedef RF RealType;

      private:

        const RF   cfl;
        const IGO& igo;

      public:

        ScalingTimeController (RF cfl_, const IGO& igo_)
          : cfl(cfl_), igo(igo_)
        {}

        /**
         * @brief Reduce step size by given CFL factor
         */
        virtual RF suggestTimestep (RF time, RF givendt)
        {
          return cfl * igo.suggestTimestep(givendt);
        }
    };

    /**
     * @brief Solver for transient PDEs using explicit timestepping scheme
     */
    template<typename Traits, typename ModelType, typename DirectionType>
      class ExplicitLinearSolver
      {
        private:

          using RF = typename Traits::GridTraits::RangeField;
          using DF = typename Traits::GridTraits::DomainField;

          using DiscType   = typename EquationTraits<Traits,ModelType,DirectionType>::DiscretizationType;
          using GFS        = typename EquationTraits<Traits,ModelType,DirectionType>::GridFunctionSpace;
          using GridVector = typename EquationTraits<Traits,ModelType,DirectionType>::GridVector;
          using C          = typename GFS::template ConstraintsContainer<RF>::Type;

          using LOP  = SpatialOperator <Traits,ModelType,DiscType,DirectionType>;
          using TLOP = TemporalOperator<Traits,ModelType,DiscType,DirectionType>;

          using Method = typename EquationTraits<Traits,ModelType,DirectionType>::OneStepScheme;

          using MBE = Dune::PDELab::istl::BCRSMatrixBackend<>;
          using GO0 = Dune::PDELab::GridOperator<GFS,GFS, LOP,MBE,DF,RF,RF,C,C>;
          using GO1 = Dune::PDELab::GridOperator<GFS,GFS,TLOP,MBE,DF,RF,RF,C,C>;
          using IGO = Dune::PDELab::OneStepGridOperator<GO0,GO1,false>;

          using LS  = Dune::PDELab::ISTLBackend_OVLP_ExplicitDiagonal<GFS>;
          using TC  = ScalingTimeController<RF,IGO>;
          using OSM = Dune::PDELab::ExplicitOneStepMethod<RF,IGO,LS,GridVector,GridVector,TC>;

          const EquationTraits<Traits,ModelType,DirectionType>& equationTraits;
          ModelParameters<Traits,ModelType>&                    parameters;

          LOP lop;
          TLOP tlop;
          const Method method;

          MBE mbe;
          GO0 go0;
          GO1 go1;
          IGO igo;
          LS  ls;

          const RF cflFactor;
          TC       tc;
          OSM      osm;

          using ERTMatrixContainer = typename Traits::ERTMatrixContainer;
          std::shared_ptr<ERTMatrixContainer> ertMatrixContainer;

        public:

          /**
           * @brief Constructor
           */
          ExplicitLinearSolver(
              const Traits& traits,
              const EquationTraits<Traits,ModelType,DirectionType>& equationTraits_,
              C& cg,
              ModelParameters<Traits,ModelType>& parameters_,
              RF Tend,
              std::shared_ptr<typename Traits::ERTMatrixContainer>& ertMatrixContainer_
              )
            : equationTraits(equationTraits_), parameters(parameters_),
            lop(traits,parameters), tlop(traits,parameters), mbe(9),
            go0(equationTraits.gfs(),cg,equationTraits.gfs(),cg,lop ,mbe),
            go1(equationTraits.gfs(),cg,equationTraits.gfs(),cg,tlop,mbe),
            igo(go0,go1), ls(equationTraits.gfs()),
            cflFactor(traits.config().get<RF>("solver.cflFactor",0.999)),
            tc(cflFactor,igo), osm(method,igo,ls,tc),
            ertMatrixContainer(ertMatrixContainer_)
        {}

          /**
           * @brief Compute a time step
           */
          unsigned int step(RF timeFrom, RF timestep, GridVector& oldSolution, GridVector& newSolution)
          {
            unsigned int count = 0;

            RF time   = timeFrom;
            RF timeTo = timeFrom + timestep;
            RF timeTaken;

            GridVector intermediateSolution = oldSolution;
            while (std::abs(timeTo - time) > 1e-6)
            {
              timeTaken = osm.apply(time,timeTo - time,intermediateSolution,newSolution);

              intermediateSolution = newSolution;
              time += timeTaken;

              count++;
            }

            return count;
          }

      };
  }
}

#endif // DUNE_MODELLING_SOLVERS_HH
