// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FTG_TRANSPORT_HH
#define DUNE_FTG_TRANSPORT_HH

#include<dune/pdelab/common/referenceelements.hh>
#include<dune/pdelab/finiteelementmap/p0fem.hh>
#include<dune/pdelab/constraints/p0.hh>
#include<dune/pdelab/backend/istl.hh>

#include<dune/modelling/fluxreconstruction.hh>
#include<dune/modelling/solutionstorage.hh>
#include<dune/ftg/override/solvers_transport.hh>

namespace Dune {
  namespace Modelling {

    /**
     * @brief Parameter class for the solute transport equation
     */
    template<typename Traits>
      class ModelParameters<Traits, typename ModelTypes::Transport>
      : public ModelParametersBase<Traits>
      {
        using RF = typename Traits::GridTraits::RangeField;

        using ParameterList  = typename Traits::ParameterList;
        using ParameterField = typename ParameterList::SubRandomField;

        const Traits& traits;

        std::shared_ptr<SolutionStorage<Traits,ModelTypes::Transport,Direction::Forward> > forwardStorage;
        std::shared_ptr<ParameterField> porosityField;
        std::shared_ptr<const ModelParameters<Traits,ModelTypes::Groundwater> > groundwaterParams;

        RF c_in_private;
        RF t_in_private;

        public:
          unsigned int model_number = 0;

        ModelParameters(const Traits& traits_, const std::string& name)
          : ModelParametersBase<Traits>(name), traits(traits_)
        {
            c_in_private = traits.config().template get<RF>("tracer.c_in");
            t_in_private = traits.config().template get<RF>("tracer.t_in");
        }

        /**
         * @brief Model parameters should exist only once per (named) model
         */
        ModelParameters(const ModelParameters& other) = delete;

        /**
         * @brief Set internal storage object for forward solution
         */
        void setStorage(
            const std::shared_ptr<SolutionStorage<Traits,ModelTypes::Transport,Direction::Forward> > storage
            )
        {
          forwardStorage = storage;
        }

        /**
         * @brief Provide access to underlying parameter fields
         */
        void setParameterList(const std::shared_ptr<const typename Traits::ParameterList>& list)
        {
          porosityField = (*list).get("porosity");
        }

        /**
         * @brief Provide access to measurement storage object
         */
        void setMeasurementList(const std::shared_ptr<const typename Traits::MeasurementList>& list)
        {}

        RF timestep() const
        {
          return traits.config().template get<RF>("time.step_transport");
        }

        /**
         * @brief Minimum time stepsize accepted by model
         */
        RF minTimestep() const
        {
          return timestep();
        }

        /**
         * @brief Maximum time stepsize accepted by model
         */
        RF maxTimestep() const
        {
          return timestep();
        }

        /**
         * @brief Porosity value at given position
         */
        template<typename Domain>
          RF porosity(const Domain& x) const
          {
            // read constant porosity from file; if not there try looking for field
            typename Traits::GridTraits::Scalar output = -5.0;
            output = traits.config().template get<RF>("parameters.porosity",output);
            if (output == -5.0)
            {
              if (!porosityField)
                DUNE_THROW(Dune::Exception,"porosity field not set in transport model parameters");
              (*porosityField).evaluate(x,output);
            }
            return output[0];
          }
	
        /**
         * @brief Concentration at given position;
         */
        template<typename Element, typename Domain, typename Time>
          RF concentration(const Element& elem, const Domain& x, const Time& time) const
          {
            typename Traits::GridTraits::Scalar localConcentration;
            (*forwardStorage).value(time,elem,x,localConcentration);
            return localConcentration[0];
          }

        /**
         * @brief Normal convective flux across interface
         */
        template<typename Intersection, typename IDomain, typename Time>
          RF vNormal (const Intersection& is, const IDomain& x, const Time& t) const
          {
            if (!groundwaterParams)
              DUNE_THROW(Dune::Exception,"groundwater model parameters not set in transport model parameters");
            
            const auto& global = is.geometry().global(x);
            //if ((global[0] == 0.0) && (global[1] == 0.25) && (global[2]==-6.5))
            //{
            //  std::cout << "Test condition reached" << std::endl;
            //}
            RF v= (*groundwaterParams).flux(is,x,t)/porosity(global);
            return v; // division by porosity!
          }


        RF c_in () const
        {
          return c_in_private;
        }

        RF t_in () const
        {
          return t_in_private;
        }

        /**
         * @brief Constant for diffusive flux contribution
         */
        RF diffusion (const RF& v) const
          {
            RF dispersivity = traits.config().template get<RF>("parameters.dispersivity");
            RF molecular_D = traits.config().template get<RF>("parameters.molecular_D");
            return fabs(v)*dispersivity+molecular_D;
          }

        /**
         * @brief Maximum velocity in cell for CFL condition
         */
        template<typename Element, typename Time>
          RF maxVelocity(const Element& elem, const Time& time) const
          {
            if (!groundwaterParams)
              DUNE_THROW(Dune::Exception,"groundwater model parameters not set in transport model parameters");
            
            return (*groundwaterParams).maxFluxNorm(elem,time)/porosity(elem.geometry().center());
          }

        template<typename Element, typename Time>
          RF dt_CFL(const Element& elem, const Time& time) const
          {
            if (!groundwaterParams)
              DUNE_THROW(Dune::Exception,"groundwater model parameters not set in transport model parameters");
            
            return (*groundwaterParams).dt_CFL(elem,time)*porosity(elem.geometry().center());
          }

        /**
         * @brief Make ModelParameters of different model available
         */
        void registerModel(
            const std::string& name, 
            const std::shared_ptr<ModelParameters<Traits,ModelTypes::Groundwater> >& otherParams
            )
        {
          groundwaterParams = otherParams;
        }

        void registerModel(
            const std::string& name,
            const std::shared_ptr<ModelParameters<Traits,ModelTypes::Geoelectrics> >& otherParams
            )
        {}

        std::map<unsigned int, std::pair<RF, bool>> well_cells() const
        {
          return traits.read_well_cells();
        };
      };

    /**
     * @brief Equation traits class for forward solute transport equation
     */
    template<typename Traits, typename DirectionType>
      class EquationTraits<Traits, typename ModelTypes::Transport, DirectionType>
      {
        public:

          using GridTraits = typename Traits::GridTraits;

          using GridView    = typename GridTraits::GridView;
          using DomainField = typename GridTraits::DomainField;
          using RangeField  = typename GridTraits::RangeField;

          enum {dim = GridTraits::dim};

          using FEM = Dune::PDELab::P0LocalFiniteElementMap<DomainField,RangeField,dim>;
          using VBE = Dune::PDELab::istl::VectorBackend<Dune::PDELab::istl::Blocking::fixed,1>;
          using CON = Dune::PDELab::P0ParallelConstraints;

          using GridFunctionSpace = Dune::PDELab::GridFunctionSpace<GridView,FEM,CON,VBE>;
          using GridVector        = typename Dune::PDELab::Backend::Vector<GridFunctionSpace,RangeField>;

          using DiscretizationType = Discretization::CellCenteredFiniteVolume;
          template<typename... T>
            using StationarySolver = StationarySolverTransport<T...>; // dummy implementation
          template<typename... T>
            //using TransientSolver  = TransientSolverTransport_Explicit<T...>;
            using TransientSolver  = TransientSolverTransport_Implicit<T...>;
            //using TransientSolver  = TransientSolverTransport_CrankNicolson<T...>;
            //using OneStepScheme = Dune::PDELab::ExplicitEulerParameter<RangeField>;
            using OneStepScheme = Dune::PDELab::ImplicitEulerParameter<RangeField>;
            //using OneStepScheme = Dune::PDELab::OneStepThetaParameter<RangeField>;
          template<typename... T>
            using FluxReconstruction = RT0Reconstruction<T...>;
          template<typename... T>
            using StorageContainer = LastTwoContainer<T...>;
          template<typename... T>
            using TemporalInterpolation = PreviousTimestep<T...>;

        private:

          const FEM fem;
          const GridFunctionSpace space;

        public:

          EquationTraits(const Traits& traits)
            : fem(Dune::GeometryType(Dune::GeometryType::cube,dim)),
            space(traits.grid().levelGridView(0),fem)
          {
            space.ordering();
          }

          const GridFunctionSpace& gfs() const
          {
            return space;
          }

      };

    /**
     * @brief Class representing initial condition for solute concentration
     */
    template<typename Traits>
      class InitialValue<Traits, typename ModelTypes::Transport>
      : public InitialValueBase<Traits,InitialValue<Traits,ModelTypes::Transport> >
      {
        public:
          template<typename GV>
            InitialValue(
                const Traits& traits_,
                const GV& gv,
                const ModelParameters<Traits,ModelTypes::Transport>& parameters_
                )
            : InitialValueBase<Traits,InitialValue<Traits,ModelTypes::Transport> >(gv)
            {}

          template<typename Domain, typename Value>
            void evaluateGlobal(const Domain& x, Value& y) const
            {
              // initialize concentration
              y[0] = 0.;
            }
      };

    template<typename Traits, typename ModelType, typename DirectionType>
      class SourceTerm;

    /**
     * @brief Source term of solute transport equation
     */
    template<typename Traits>
      class SourceTerm<Traits, ModelTypes::Transport, Direction::Forward>
      {
        const Traits & traits;
        const ModelParameters<Traits,ModelTypes::Transport>& parameters;
        typename Traits::GridTraits::Grid::LeafGridView  lgv;
        using RF = typename Traits::GridTraits::RangeField;
        using IDomain = typename Traits::GridTraits::IDomain;
        //RF c_in;
        //RF t_in;
        std::map<unsigned int, std::pair<RF, bool>> well_cells;

        public:

          SourceTerm(const Traits& traits_, const ModelParameters<Traits,ModelTypes::Transport>& parameters_) : traits(traits_), parameters(parameters_), lgv(traits_.grid().leafGridView()) 
          {
            //c_in = traits.config().template get<RF>("tracer.c_in");
            //t_in = traits.config().template get<RF>("tracer.t_in");
            well_cells = parameters.well_cells();
          }

          template<typename Element, typename Domain>
            RF Q (const Element& elem, const Domain& x, bool & tracer_injection) const
            {
              // get the index of the current cell
              unsigned int current_index = lgv.indexSet().index(elem);
              auto temp = well_cells.find(current_index);

              RF rate = 0.0;

              if ( !(temp->first == well_cells.end()->first) ) //  current cell is a well cell 
              {  
                rate = (temp->second).first;
                if ((temp->second).second == true)
                  tracer_injection = true;
              }
              return rate;
            }
      };

    /**
     * @brief Define transport equation as a differential equation
     */
    template<typename Traits, typename DomainType, typename DirectionType>
      class Equation<Traits,ModelTypes::Transport,DomainType,DirectionType>
      : public DifferentialEquation<Traits,ModelTypes::Transport,DomainType,DirectionType>
      {
        using DifferentialEquation<Traits,ModelTypes::Transport,DomainType,DirectionType>::DifferentialEquation;
      };

    /**
     * @brief Spatial local operator of the convection-diffusion equation (CCFV version)
     */
    template<typename Traits, typename DirectionType>
      class SpatialOperator<Traits, typename ModelTypes::Transport,
            Discretization::CellCenteredFiniteVolume, DirectionType>
      : 
      public Dune::PDELab::NumericalJacobianSkeleton<SpatialOperator<Traits, ModelTypes::Transport, Discretization::CellCenteredFiniteVolume, DirectionType> >,
      public Dune::PDELab::NumericalJacobianVolume<SpatialOperator<Traits, ModelTypes::Transport, Discretization::CellCenteredFiniteVolume, DirectionType> >,
      public Dune::PDELab::NumericalJacobianBoundary<SpatialOperator<Traits, ModelTypes::Transport, Discretization::CellCenteredFiniteVolume, DirectionType> >,
      public Dune::PDELab::FullSkeletonPattern, 
      public Dune::PDELab::FullVolumePattern,
      public Dune::PDELab::LocalOperatorDefaultFlags,
      public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<typename Traits::GridTraits::RangeField>
      {
        using GridTraits = typename Traits::GridTraits;

        using RF      = typename GridTraits::RangeField;
        using Domain  = typename GridTraits::Domain;
        using IDomain = typename GridTraits::IDomain;

        public:

        // pattern assembly flags
        enum {doPatternVolume   = true};
        enum {doPatternSkeleton = true};

        // residual assembly flags
        enum {doAlphaVolume   = true};
        enum {doAlphaSkeleton = true};
        enum {doAlphaBoundary = true};
        enum {doLambdaVolume = true};
        //enum {doLambdaSkeleton = false};
        enum {doLambdaBoundary = true};

        private:

        const Traits& traits;

        const ModelParameters<Traits,ModelTypes::Transport>&              parameters;
        const Boundary       <Traits,ModelTypes::Transport,DirectionType> boundary;
        const SourceTerm     <Traits,ModelTypes::Transport,DirectionType> sourceTerm;

        RF time;

        //mutable bool firstStage = true;
        //mutable RF dtmin;
        //mutable RF dt_min_CFL;

        public:

        /**
         * @brief Constructor
         */
        SpatialOperator(
            const Traits& traits_,
            const ModelParameters<Traits,ModelTypes::Transport>& parameters_
            )
          : traits(traits_), parameters(parameters_), boundary(traits,parameters.name()), sourceTerm(traits,parameters)
        {}

        /**
         * @brief Volume integral depending on test and ansatz functions
         */
        template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
          void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
          {
            // contribution from sink term
            const Domain& cellCenterLocal = referenceElement(eg.geometry()).position(0,0);
            bool tracer_injection = false;
            RF Q = sourceTerm.Q(eg.entity(),cellCenterLocal,tracer_injection);
            
            const RF source = (Q < 0) ? -Q*x(lfsu,0)/parameters.porosity(cellCenterLocal) : 0;
            r.accumulate(lfsv,0,source);
          }

        /**
         * @brief Volume integral depending only on test functions (source term)
         */
        template<typename EG, typename LFSV, typename R>
          void lambda_volume (const EG& eg, const LFSV& lfsv, R& r) const
          {
            // contribution from source term
            
            RF source = 0.0;
            if (time <= parameters.t_in())
            {
              const Domain& cellCenterLocal = referenceElement(eg.geometry()).position(0,0);
              bool tracer_injection = false;
              RF Q = sourceTerm.Q(eg.entity(),cellCenterLocal,tracer_injection);
              if (tracer_injection==true)
                source = (Q > 0) ? -Q*parameters.c_in()/parameters.porosity(cellCenterLocal) : 0;
            } 
            r.accumulate(lfsv,0,source);
          }

        template<typename IG, typename LFSV, typename R>
          void lambda_boundary (const IG& ig, const LFSV& lfsv, R& r_s) const
          {
            /*
            // contribute to timestep calculation
            if (firstStage)
            {
              // compute minimum local mesh width h
              const RF faceVolume   = ig         .geometry().volume();
              const RF insideVolume = ig.inside().geometry().volume();
                  
              const RF distance = insideVolume/faceVolume;

              // compute maximum local velocity |v|
              const RF maxVelocity = parameters.maxVelocity(ig.inside(),time);
              RF dt_CFL = parameters.dt_CFL(ig.inside(),time);
              // update minimum of admissable timesteps
              dtmin = std::min(dtmin,distance/maxVelocity);
              dt_min_CFL = std::min(dt_min_CFL,dt_CFL);
            }
            */

            // evaluate boundary condition type
            const IDomain& faceCenterLocal = referenceElement(ig.geometry()).position(0,0);
            Dune::Modelling::BoundaryCondition::Type bc = boundary.bc(ig.intersection(),faceCenterLocal,time);

            // Neumann boundary condition
            if (Dune::Modelling::BoundaryCondition::isNeumann(bc))
            {
              const RF faceVolume = ig.geometry().volume();
              const RF j = boundary.j(ig.intersection(),faceCenterLocal,time);
              r_s.accumulate(lfsv,0, j * faceVolume);
            }
          }

        /**
         * @brief Skeleton integral depending on test and ansatz functions
         */
        // each face is only visited ONCE!
        template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
          void alpha_skeleton (const IG& ig, 
              const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
              const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n, 
              R& r_s, R& r_n) const
          {
            /*
            // contribute to timestep calculation
            if (firstStage)
            {
              // compute minimum local mesh width h
              const RF faceVolume    = ig          .geometry().volume();
              const RF insideVolume  = ig.inside() .geometry().volume();
              const RF outsideVolume = ig.outside().geometry().volume();
                  
              const RF distance = havg(insideVolume,outsideVolume)/faceVolume;

              // compute maximum local velocity |v|
              const RF insideMaxVelocity  = parameters.maxVelocity(ig.inside() ,time);
              const RF outsideMaxVelocity = parameters.maxVelocity(ig.outside(),time);
              
              const RF maxVelocity = std::max(insideMaxVelocity,outsideMaxVelocity);

              RF inside_dt_CFL = parameters.dt_CFL(ig.inside(),time);
              RF outside_dt_CFL = parameters.dt_CFL(ig.outside(),time);

              RF dt_CFL = std::min(inside_dt_CFL,outside_dt_CFL);
              
              // update minimum of admissable timesteps
              dt_min_CFL = std::min(dt_min_CFL,dt_CFL);
              dtmin = std::min(dtmin,distance/maxVelocity);
            }
            */

            const RF normalFlux = skeletonNormalFlux(ig.intersection(),x_s(lfsu_s,0),x_n(lfsu_n,0),time);
            const RF faceVolume = ig.geometry().volume();

            r_s.accumulate(lfsv_s,0,   normalFlux * faceVolume);
            r_n.accumulate(lfsv_n,0, - normalFlux * faceVolume);
          }

        /**
         * @brief Boundary integral depending on test and ansatz functions
         */
        // We put the Dirichlet evaluation also in the alpha term to save some geometry evaluations
        template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
          void alpha_boundary (const IG& ig, 
              const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
              R& r_s) const
          {
            /*
            // contribute to timestep calculation
            if (firstStage)
            {
              // compute minimum local mesh width h
              const RF faceVolume   = ig         .geometry().volume();
              const RF insideVolume = ig.inside().geometry().volume();
                  
              const RF distance = insideVolume/faceVolume;

              // compute maximum local velocity |v|
              const RF maxVelocity = parameters.maxVelocity(ig.inside(),time);
              RF dt_CFL = parameters.dt_CFL(ig.inside(),time);
              // update minimum of admissable timesteps
              dtmin = std::min(dtmin,distance/maxVelocity);
              dt_min_CFL = std::min(dt_min_CFL,dt_CFL);
            }
            */

            // evaluate boundary condition type           
            const IDomain& faceCenterLocal = referenceElement(ig.geometry()).position(0,0);
            Dune::Modelling::BoundaryCondition::Type bc = boundary.bc(ig.intersection(),faceCenterLocal,time);

            // Dirichlet boundary condition
            if (Dune::Modelling::BoundaryCondition::isDirichlet(bc))
            {
              const RF normalFlux = dirichletNormalFlux(ig.intersection(),x_s(lfsu_s,0),time);
              const RF faceVolume = ig.geometry().volume();

              r_s.accumulate(lfsv_s,0, normalFlux * faceVolume);
            }
          }

        /**
         * @brief Set up CFL timestep calculation in first stage
         */
        void preStage(RF time, int r)
        {
          /*
          if (r == 1)
          {
            //firstStage = true; <- not needed here, as groundwater flow is stationary; one evaluation sufficient
            dtmin = std::numeric_limits<RF>::max();
            dt_min_CFL = std::numeric_limits<RF>::max();
          }
          else
            firstStage = false;
          */
        }

        /**
         * @brief Set time for subsequent evaluation
         */
        void setTime (const RF t)
        {
          time = t;
        }

        /**
         * @brief Flux across intersection based on DGF values
         */
        template<typename Intersection, typename IDomain, typename ScalarDGF, typename VectorDGF>
          RF normalFlux(
              const Intersection& is,
              const IDomain& x,
              const ScalarDGF& scalar,
              const VectorDGF& localGrad,
              RF time
              ) const
          {
            if (is.neighbor())
            {
              const Domain& cellCenterInside  = referenceElement(is.inside() .geometry()).position(0,0);
              const Domain& cellCenterOutside = referenceElement(is.outside().geometry()).position(0,0);

              typename Traits::GridTraits::Scalar innerValue, outerValue;
              (*scalar).evaluate(is.inside(), cellCenterInside,innerValue);
              (*scalar).evaluate(is.outside(),cellCenterOutside,outerValue);

              return skeletonNormalFlux(is,innerValue[0],outerValue[0],time);
            }
            else if (is.boundary())
            {
              const IDomain& faceCenterLocal = referenceElement(is.geometry()).position(0,0);

              // evaluate boundary condition type
              Dune::Modelling::BoundaryCondition::Type bc = boundary.bc(is,faceCenterLocal,time);

              if (Dune::Modelling::BoundaryCondition::isDirichlet(bc))
              {
                const Domain& cellCenterInside = referenceElement(is.inside().geometry()).position(0,0);
                typename Traits::GridTraits::Scalar innerValue;
                (*scalar).evaluate(is.inside(), cellCenterInside,innerValue);
                return dirichletNormalFlux(is,innerValue[0],time);
              }
              else if (Dune::Modelling::BoundaryCondition::isNeumann(bc))
                return boundary.j(is,faceCenterLocal,time);
              else
                DUNE_THROW(Dune::Exception,"unknown boundary condition type");
            }
            else
              return 0.;
          }

        /**
         * @brief Derivative of flux across intersection based on DGF values
         */
        template<typename Intersection, typename IDomain, typename ScalarDGF, typename VectorDGF>
          RF normalFluxDeriv(
              const Intersection& is,
              const IDomain& x,
              const ScalarDGF& scalar,
              const VectorDGF& localGrad,
              RF time
              ) const
          {
            DUNE_THROW(Dune::NotImplemented,
                "derivative of normal flux not implemented for groundwater flow equation");
          }

        /**
         * @brief Derivative across intersection based on DGF values
         */
        template<typename Intersection, typename IDomain, typename ScalarDGF, typename VectorDGF>
          RF normalDerivative(
              const Intersection& is,
              const IDomain& x,
              const ScalarDGF& scalar,
              const VectorDGF& localGrad,
              RF time
              ) const
          {
            if (is.neighbor())
            {
              const Domain& cellCenterInside  = referenceElement(is.inside() .geometry()).position(0,0);
              const Domain& cellCenterOutside = referenceElement(is.outside().geometry()).position(0,0);

              typename Traits::GridTraits::Scalar innerValue, outerValue;
              (*scalar).evaluate(is.inside(), cellCenterInside,innerValue);
              (*scalar).evaluate(is.outside(),cellCenterOutside,outerValue);

              return skeletonNormalDerivative(is,innerValue[0],outerValue[0]);
            }
            else if (is.boundary())
            {
              const IDomain& faceCenterLocal = referenceElement(is.geometry()).position(0,0);

              // evaluate boundary condition type
              Dune::Modelling::BoundaryCondition::Type bc = boundary.bc(is,faceCenterLocal,time);

              if (Dune::Modelling::BoundaryCondition::isDirichlet(bc))
              {
                const Domain& cellCenterInside = referenceElement(is.inside().geometry()).position(0,0);
                typename Traits::GridTraits::Scalar innerValue;
                (*scalar).evaluate(is.inside(), cellCenterInside,innerValue);
                return dirichletNormalDerivative(is,innerValue[0],time);
              }
              else if (Dune::Modelling::BoundaryCondition::isNeumann(bc))
              {
                const Domain& cellCenterInside = referenceElement(is.inside().geometry()).position(0,0);
                const RF boundaryFlux = boundary.j(is,faceCenterLocal,time);
                typename Traits::GridTraits::Scalar innerValue;
                (*scalar).evaluate(is.inside(), cellCenterInside,innerValue);
                const RF v = parameters.vNormal(is,faceCenterLocal,time);
                const RF D_inside = parameters.diffusion(v);
                // revert flux computation to obtain gradient on boundary
                return - (boundaryFlux - v * innerValue[0]) / D_inside;
              }
              else
                DUNE_THROW(Dune::Exception,"unknown boundary condition type");
            }
            else
              return 0.;
          }

        /**
         * @brief Suggest timestep based on CFL condition
         */
        /*RF suggestTimestep(RF dt) const
        {
            //std::cout << "Suggested time step: " << traits.comm().min(dt_min_CFL) << "s" << std::endl; // <- the other CFL implementation works better
            if (firstStage)
            {
              std::cout << "CFL-based time step suggestion: " << traits.comm().min(dtmin) << "s" << std::endl;
              firstStage = false;
            }
            return traits.comm().min(dt_min_CFL);
        }*/

        private:

        /**
         * @brief Flux in normal direction across interface
         */
        template<typename Intersection>
          RF skeletonNormalFlux(const Intersection& is, RF innerValue, RF outerValue, RF time) const
          {
            // geometry information
            const IDomain& faceCenterLocal  = referenceElement(is.geometry()).position(0,0);

            // advection velocity
            const RF v = parameters.vNormal(is,faceCenterLocal,time);

            // upwinding
            const RF upwindValue = (v >= 0) ? innerValue : outerValue;

            // conductivity
            const RF D_inside  = parameters.diffusion(v);
            const RF D_outside = parameters.diffusion(v);
            const RF D = havg(D_inside,D_outside);

            return v * upwindValue - D * skeletonNormalDerivative(is,innerValue,outerValue);
          }

        /**
         * @brief Derivative of solution in normal direction across interface
         */
        template<typename Intersection>
          RF skeletonNormalDerivative(const Intersection& is, RF innerValue, RF outerValue) const
          {
            // geometry information
            const RF faceVolume    = is          .geometry().volume();
            const RF insideVolume  = is.inside() .geometry().volume();
            const RF outsideVolume = is.outside().geometry().volume();

            const RF distance = havg(insideVolume,outsideVolume)/faceVolume;

            // approximation of gradient
            const RF w = (outerValue - innerValue)/distance;

            return w;
          }

        /**
         * @brief Flux in normal direction on Dirichlet boundary
         */
        template<typename Intersection>
          RF dirichletNormalFlux(const Intersection& is, RF innerValue, RF time) const
          {
            // geometry information
            const IDomain& faceCenterLocal  = referenceElement(is.geometry()).position(0,0);
            //const Domain&  cellCenterInside = referenceElement(is.inside().geometry()).position(0,0);
            
            // advection velocity
            const RF v = parameters.vNormal(is,faceCenterLocal,time);

            // upwinding
            const RF upwindValue = (v >= 0) ? innerValue : boundary.g(is,faceCenterLocal,time);

            // diagonal diffusion tensor
            const RF D_inside = parameters.diffusion(v);

            return v * upwindValue - D_inside * dirichletNormalDerivative(is,innerValue,time);
          }

        /**
         * @brief Derivative of solution in normal direction on Dirichlet boundary
         */
        template<typename Intersection>
          RF dirichletNormalDerivative(const Intersection& is, RF innerValue, RF time) const
          {
            // geometry information
            const IDomain& faceCenterLocal  = referenceElement(is.geometry()).position(0,0);

            const RF faceVolume    = is          .geometry().volume();
            const RF insideVolume  = is.inside() .geometry().volume();

            const RF distance = insideVolume/faceVolume;

            // approximation of gradient
            const RF g = boundary.g(is,faceCenterLocal,time);
            const RF w = (g - innerValue)/(distance/2.);

            return w;
          }

        /**
         * @brief Harmonic average
         */
        template<typename T>
          T havg (T a, T b) const
          {
            T eps = 1e-30;
            return 2./(1./(a + eps) + 1./(b + eps));
          }

      };

    /**
     * @brief Temporal local operator of the convection-diffusion equation (CCFV version)
     */
    template<typename Traits, typename DirectionType>
      class TemporalOperator<Traits, typename ModelTypes::Transport,
            Discretization::CellCenteredFiniteVolume, DirectionType>
      : public Dune::PDELab::NumericalJacobianVolume
      <TemporalOperator<Traits, ModelTypes::Transport, Discretization::CellCenteredFiniteVolume, DirectionType> >,
      public Dune::PDELab::FullVolumePattern,
      public Dune::PDELab::LocalOperatorDefaultFlags,
      public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<typename Traits::GridTraits::RangeField>

      {
        using RF = typename Traits::GridTraits::RangeField;

        public:

        // pattern assembly flags
        enum {doPatternVolume = true};

        // residual assembly flags
        enum {doAlphaVolume = true};

        private:

        public:

        /**
         * @brief Constructor
         */
        TemporalOperator(
            const Traits& traits,
            const ModelParameters<Traits,ModelTypes::Transport>& parameters
            ) 
        {}

        /**
         * @brief Volume integral depending on test and ansatz functions
         */
        template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
          void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
          {
            r.accumulate(lfsv,0, x(lfsu,0) * eg.geometry().volume());
          }

        RF suggestTimestep (RF dt) const
        {
          return std::numeric_limits<RF>::max();
        }

      };
  }
}

#endif // DUNE_FTG_TRANSPORT_HH
