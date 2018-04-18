// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MODELLING_TRANSPORT_HH
#define DUNE_MODELLING_TRANSPORT_HH

#include<dune/pdelab/common/referenceelements.hh>
#include<dune/pdelab/finiteelementmap/p0fem.hh>
#include<dune/pdelab/constraints/p0.hh>
#include<dune/pdelab/backend/istl.hh>

#include<dune/modelling/fluxreconstruction.hh>
#include<dune/modelling/solutionstorage.hh>
#include<dune/modelling/solvers.hh>

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
        std::shared_ptr<SolutionStorage<Traits,ModelTypes::Transport,Direction::Adjoint> > adjointStorage;

        std::shared_ptr<ParameterField> porosityField;

        std::shared_ptr<const ModelParameters<Traits,ModelTypes::Groundwater> > groundwaterParams;

        public:

        ModelParameters(const Traits& traits_, const std::string& name)
          : ModelParametersBase<Traits>(name), traits(traits_)
        {}

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
         * @brief Set internal storage object for adjoint solution
         */
        void setStorage(
            const std::shared_ptr<SolutionStorage<Traits,ModelTypes::Transport,Direction::Adjoint> > storage
            )
        {
          adjointStorage = storage;
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
        {
          //std::cout << "this is the setMeasurementList method of the transport model" << std::endl;
        }

        /**
         * @brief Minimum time stepsize accepted by model
         */
        RF minTimestep() const
        {
          return traits.config().template get<RF>("time.minStep");
        }

        /**
         * @brief Maximum time stepsize accepted by model
         */
        RF maxTimestep() const
        {
          return traits.config().template get<RF>("time.maxStep");
        }

        /**
         * @brief Porosity value at given position
         */
        template<typename Domain>
          RF porosity(const Domain& x) const
          {
            if (!porosityField)
              DUNE_THROW(Dune::Exception,"porosity field not set in transport model parameters");

            typename Traits::GridTraits::Scalar output;
            (*porosityField).evaluate(x,output);
            return output[0];
          }
	
	      /**
         * @brief Concentration at given position;
         */
        template<typename Element, typename Domain, typename Time>
          RF concentration(const Element& elem, const Domain& x, const Time& time) const
          {
            typename Traits::GridTraits::Scalar localConcentration;
            //const auto& global  = elem.geometry().global(x);
            //const auto& xInside = elem.geometry().local(global);
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
            return (*groundwaterParams).flux(is,x,t)/porosity(global); // !!! division by porosity !!!
          }

        /**
         * @brief Constant for diffusive flux contribution
         */
        template<typename Element, typename Domain, typename Time>
          RF diffusion (const Element& elem, const Domain& x, const Time& t) const
          {
            return traits.config().template get<RF>("parameters.diffusion");
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

        /**
         * @brief Adjoint source term based on measurements
         */
        template<typename Element, typename Domain, typename Time>
          RF adjointSource(const Element& elem, const Domain& x, const Time& t) const
          {
            // only needed for adjoint
            return 0.;
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
        {
          // only needed for adjoint?
        }
        
        std::vector<int> tracer_cell_indices() const
        {
          return traits.read_tracer_cell_indices();
        };

      };

    /**
     * @brief Equation traits class for forward / adjoint solute transport equation
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

          // use linear solver in stationary case,
          // explicit linear solver for transient case
          template<typename... T>
            using StationarySolver = StationaryLinearSolver<T...>;
          template<typename... T>
            using TransientSolver  = ExplicitLinearSolver<T...>;

          // use explicit Euler for timestepping
          // alternative: Heun
          using OneStepScheme = Dune::PDELab::ExplicitEulerParameter<RangeField>;

          // use RT0 flux reconstruction
          // alternatives: BDM1 and RT1
          template<typename... T>
            using FluxReconstruction = RT0Reconstruction<T...>;

          // store complete space-time solution
          // alternative: only store last two steps
          template<typename... T>
            using StorageContainer = FullContainer<T...>;

          // use previous timestep when interpolating stored solution
          // alternatives: NextTimestep and LinearInterpolation
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
        /*
        const ModelParameters<Traits,ModelTypes::Transport>& parameters;
        const Traits& traits;
        typename Traits::GridTraits::Grid::LeafGridView lgv;
        using RF = typename Traits::GridTraits::RangeField;
        
        std::vector<double> x_min;
        std::vector<double> x_max;
        std::vector<double> y_min;
        std::vector<double> y_max;
        std::vector<double> z_min;
        std::vector<double> z_max;
        std::vector<double> volume;
        double totalvolume = 0.;
        RF tracermass;
        */
        public:
          template<typename GV>
            InitialValue(
                const Traits& traits_,
                const GV& gv,
                const ModelParameters<Traits,ModelTypes::Transport>& parameters_
                )
            : InitialValueBase<Traits,InitialValue<Traits,ModelTypes::Transport> >(gv)//, parameters(parameters_), traits(traits_), lgv(traits.grid().leafGridView())
            {
              /*
              // as the evaluateGlobal function only receives global coordinates, we need to access all injection cells and get their global
              // coordinate limits...
              tracermass = traits.config().template get<RF>("tracer.mass");
              auto &index_set = lgv.indexSet();
              for (const auto & elem : elements (lgv))
              {
                int the_index = index_set.index(elem);
                for (unsigned int i = 0; i!=parameters.well_in_cells().size();i++)
                {
                  if (the_index == parameters.well_in_cells()[i]) 
                  {
                    //std::cout << "the index is " << the_index << std::endl;
                    x_min.push_back(elem.geometry().corner(0)[0]);
                    x_max.push_back(elem.geometry().corner(1)[0]);
                    y_min.push_back(elem.geometry().corner(0)[1]);
                    y_max.push_back(elem.geometry().corner(2)[1]);
                    z_min.push_back(elem.geometry().corner(0)[2]);
                    z_max.push_back(elem.geometry().corner(4)[2]);
                    volume.push_back(elem.geometry().volume());
                    totalvolume += elem.geometry().volume();
                  }
                }
              }
              */
            }

          template<typename Domain, typename Value>
            void evaluateGlobal(const Domain& x, Value& y) const
            {
              // initialize concentration
              y[0] = 0.;
              /*
              // if we have a tracer injection well, we have to set y<>0 according to the user specified tracer mass, porosity, etc.
              for (unsigned int i = 0; i!=x_min.size();i++)
              {
                if (x[0] >= x_min[i] && x[0] < x_max[i]) 
                {
                  if (x[1] >= y_min[i] && x[1] < y_max[i]) 
                  {
                    if (x[2] >= z_min[i] && x[2] < z_max[i]) 
                    {
                      // this gives problems if run in parallel and set to nonzero... why?
                      y[0] = 10.;//tracermass/totalvolume/parameters.porosity(x);
                      break;
                    }
                  }
                }
              } */
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
        RF tracermass;
        
        public:

          SourceTerm(const Traits& traits_, const ModelParameters<Traits,ModelTypes::Transport>& parameters_) : traits(traits_), parameters(parameters_), lgv(traits_.grid().leafGridView()) 
          {
            tracermass = traits.config().template get<RF>("tracer.mass");
          }

          template<typename Element, typename Domain, typename Value, typename Time>
            auto q (const Element& elem, const Domain& x, const Value& value, const Time& t) const
            {
              
              // if we are in one of the injection cells... and time is beginning
              // -> add a source term in such a way that the released mass equals tracer.mass
              if (t == 0)
              {
                // get the index of the current cell
                int the_index = lgv.indexSet().index(elem);

                for (unsigned int i = 0; i!=parameters.tracer_cell_indices().size();i++)
                {
                  // check if the index matches with one of the injection well indices
                  if (the_index == parameters.tracer_cell_indices()[i]) 
                  {
                    // if so, evaluate the water flux leaving the cell
                    RF v = 0.;
                    for (const auto& intersection : intersections(lgv,elem))
                    {
                      const IDomain& faceCenterLocal = referenceElement(intersection.geometry()).position(0,0);
                      v += parameters.porosity(intersection.geometry().global(faceCenterLocal)) * intersection.geometry().volume() * parameters.vNormal(intersection,faceCenterLocal,t);
                    }
                    return tracermass/v;  //TODO this is not yet correct
                  }
                }
                // if we end up here, t=0 but we are not in an injection cell -> return 0.
                return 0.;
              } else 
              {
                // no source term
                return 0.;
              }
            }
      };

    /**
     * @brief Source term of adjoint solute transport equation
     */
    /*
    template<typename Traits>
      class SourceTerm<Traits, ModelTypes::Transport, Direction::Adjoint>
      {
        const ModelParameters<Traits,ModelTypes::Transport>& parameters;

        public:

        SourceTerm(const ModelParameters<Traits,ModelTypes::Transport>& parameters_)
          : parameters(parameters_)
        {}

        template<typename Element, typename Domain, typename Value, typename Time>
          auto q (const Element& elem, const Domain& x, const Value& value, const Time& t) const
          {
            // adjoint source depends on parameters / measurements
            return parameters.adjointSource(elem,x,t);
          }
      };
    */
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
      : public Dune::PDELab::NumericalJacobianSkeleton<SpatialOperator
      <Traits, ModelTypes::Transport, Discretization::CellCenteredFiniteVolume, DirectionType> >,
      public Dune::PDELab::NumericalJacobianBoundary
        <SpatialOperator<Traits, ModelTypes::Transport, Discretization::CellCenteredFiniteVolume, DirectionType> >,
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

        private:

        const Traits& traits;

        const ModelParameters<Traits,ModelTypes::Transport>&              parameters;
        const Boundary       <Traits,ModelTypes::Transport,DirectionType> boundary;
        const SourceTerm     <Traits,ModelTypes::Transport,DirectionType> sourceTerm;

        RF adjointSign;
        RF time;

        mutable bool firstStage;
        mutable RF   dtmin;

        public:

        /**
         * @brief Constructor
         */
        SpatialOperator(
            const Traits& traits_,
            const ModelParameters<Traits,ModelTypes::Transport>& parameters_
            )
          : traits(traits_), parameters(parameters_), boundary(traits,parameters.name()), sourceTerm(traits,parameters)
        {
          if (DirectionType::isAdjoint())
            adjointSign = -1.;
          else
            adjointSign = 1.;
        }

        /**
         * @brief Volume integral depending on test and ansatz functions
         */
        template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
          void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
          {
            // contribution from source term
            const Domain& cellCenterLocal = referenceElement(eg.geometry()).position(0,0);
            r.accumulate(lfsv,0,-sourceTerm.q(eg.entity(),cellCenterLocal,x(lfsu,0),time));
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

              // update minimum of admissable timesteps
              dtmin = std::min(dtmin,distance/maxVelocity);
            }

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
            // contribute to timestep calculation
            if (firstStage)
            {
              // compute minimum local mesh width h
              const RF faceVolume   = ig         .geometry().volume();
              const RF insideVolume = ig.inside().geometry().volume();
                  
              const RF distance = insideVolume/faceVolume;

              // compute maximum local velocity |v|
              const RF maxVelocity = parameters.maxVelocity(ig.inside(),time);
              
              // update minimum of admissable timesteps
              dtmin = std::min(dtmin,distance/maxVelocity);
            }
            
            const IDomain& faceCenterLocal = referenceElement(ig.geometry()).position(0,0);
            const RF faceVolume = ig.geometry().volume();

            // evaluate boundary condition type
            Dune::Modelling::BoundaryCondition::Type bc = boundary.bc(ig.intersection(),faceCenterLocal,time);

            // Neumann boundary condition
            if (Dune::Modelling::BoundaryCondition::isNeumann(bc))
            {
              const RF j = boundary.j(ig.intersection(),faceCenterLocal,time);
              r_s.accumulate(lfsu_s,0, j * faceVolume);
              
              return;
            }

            // Dirichlet boundary condition
            if (Dune::Modelling::BoundaryCondition::isDirichlet(bc))
            {
              const RF normalFlux = dirichletNormalFlux(ig.intersection(),x_s(lfsu_s,0),time);
              const RF faceVolume = ig.geometry().volume();

              r_s.accumulate(lfsv_s,0, normalFlux * faceVolume);

              return;
            }

            // Outflow boundary condition

            // advection velocity
            const RF v = parameters.vNormal(ig.intersection(),faceCenterLocal,time);
            r_s.accumulate(lfsu_s,0, v * x_s(lfsu_s,0) * faceVolume);
          }

        /**
         * @brief Set up CFL timestep calculation in first stage
         */
        void preStage(RF time, int r)
        {
          if (r == 1)
          {
            firstStage = true;
            dtmin = std::numeric_limits<RF>::max();
          }
          else
            firstStage = false;
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
                const RF D_inside = parameters.diffusion(is.inside(),cellCenterInside,time);
                typename Traits::GridTraits::Scalar innerValue;
                (*scalar).evaluate(is.inside(), cellCenterInside,innerValue);
                const RF v = adjointSign * parameters.vNormal(is,faceCenterLocal,time);
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
        RF suggestTimestep(RF dt) const
        {
          if (dt*dtmin > 0.)
            return   traits.comm().max(dtmin);
          else
            return - traits.comm().max(dtmin);
        }

        private:

        /**
         * @brief Flux in normal direction across interface
         */
        template<typename Intersection>
          RF skeletonNormalFlux(const Intersection& is, RF innerValue, RF outerValue, RF time) const
          {
            // geometry information
            const IDomain& faceCenterLocal  = referenceElement(is.geometry()).position(0,0);
            const Domain& cellCenterInside  = referenceElement(is.inside() .geometry()).position(0,0);
            const Domain& cellCenterOutside = referenceElement(is.outside().geometry()).position(0,0);

            // advection velocity
            const RF v = adjointSign * parameters.vNormal(is,faceCenterLocal,time);

            // upwinding
            const RF upwindValue = (v >= 0) ? innerValue : outerValue;

            // conductivity
            const RF D_inside  = parameters.diffusion(is.inside(), cellCenterInside, time);
            const RF D_outside = parameters.diffusion(is.outside(),cellCenterOutside,time);
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
            const Domain&  cellCenterInside = referenceElement(is.inside().geometry()).position(0,0);
            
            // advection velocity
            const RF v = parameters.vNormal(is,faceCenterLocal,time);

            // upwinding
            const RF upwindValue = (v >= 0) ? innerValue : boundary.g(is,faceCenterLocal,time);

            // diagonal diffusion tensor
            const RF D_inside = parameters.diffusion(is.inside(),cellCenterInside,time);

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

        RF adjointSign;

        public:

        /**
         * @brief Constructor
         */
        TemporalOperator(
            const Traits& traits,
            const ModelParameters<Traits,ModelTypes::Transport>& parameters
            ) 
        {
          if (DirectionType::isAdjoint())
            adjointSign = -1.;
          else
            adjointSign = 1.;
        }

        /**
         * @brief Volume integral depending on test and ansatz functions
         */
        template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
          void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
          {
            r.accumulate(lfsv,0, adjointSign * x(lfsu,0) * eg.geometry().volume());
          }

        RF suggestTimestep (RF dt) const
        {
          return std::numeric_limits<RF>::max();
        }

      };

    /**
     * @brief Class providing sensitivity computation for solute transport
     */
    template<typename Traits>
      class SensitivityComp<Traits,ModelTypes::Transport>
      {
        public:

          SensitivityComp(
              const Traits& traits,
              const ModelParameters<Traits,ModelTypes::Transport>& parameters
              )
          {}

            void initialize(const std::shared_ptr<typename Traits::SensitivityList>& sensitivityList)
            {
            }

            template<typename Time>
            void extract(const Time& first, const Time& last)
            {
              std::cout << "(would extract sensitivity from " << first
                << " to " << last << ")" << std::endl;
            }

      };

  }
}

#endif // DUNE_MODELLING_TRANSPORT_HH
