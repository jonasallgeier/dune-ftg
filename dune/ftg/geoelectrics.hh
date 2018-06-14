// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FTG_GEOELECTRICS_HH
#define DUNE_FTG_GEOELECTRICS_HH

#include<dune/pdelab/common/referenceelements.hh>
#include<dune/pdelab/finiteelementmap/p0fem.hh>
#include<dune/pdelab/constraints/p0.hh>
#include<dune/pdelab/backend/istl.hh>

#include<dune/modelling/fluxreconstruction.hh>
#include<dune/modelling/solutionstorage.hh>
#include<dune/ftg/override/solvers.hh>
#include<dune/grid/utility/hierarchicsearch.hh>

namespace Dune {
  namespace Modelling {
    /**
     * @brief Parameter class for the geoelectrics equation
     */
    template<typename Traits>
      class ModelParameters<Traits, ModelTypes::Geoelectrics>
      : public ModelParametersBase<Traits>
      {
        using RF = typename Traits::GridTraits::RangeField;
        using ParameterList   = typename Traits::ParameterList;
        using ParameterField  = typename ParameterList::SubRandomField;
        using MeasurementList = typename Traits::MeasurementList;

        const Traits& traits;

        std::shared_ptr<SolutionStorage<Traits,ModelTypes::Geoelectrics,Direction::Forward> > forwardStorage;
        //std::shared_ptr<SolutionStorage<Traits,ModelTypes::Geoelectrics,Direction::Adjoint> > adjointStorage;

        std::shared_ptr<ParameterField> sigma_bgField;
        std::shared_ptr<ParameterField> kappaField;

        std::shared_ptr<const ModelParameters<Traits,ModelTypes::Transport> > transportParams;
        
        unsigned int electrode_cell_indx;
        public:
          unsigned int model_number;
        
        ModelParameters(const Traits& traits_, const std::string& name)
          : ModelParametersBase<Traits>(name), traits(traits_)
        {
          // access the name of the model to know which electrode is active
          std::string model_name = name;
          std::string common_base = "ERT_";
          auto start_position_to_erase = model_name.find(common_base);
          // get rid of "ERT_" to get the number n of the model "ERT_n"
          model_number = std::stoi(model_name.erase(start_position_to_erase, common_base.size())); 
          std::string model_name_temp = name;
                 
          electrode_cell_indx = traits.read_electrode_cell_indices()[model_number];
        }

        /**
         * @brief Model parameters should exist only once per (named) model
         */
        ModelParameters(const ModelParameters& other) = delete;

        /**
         * @brief Set internal storage object for forward solution
         */
        void setStorage(
            const std::shared_ptr<SolutionStorage<Traits,ModelTypes::Geoelectrics,Direction::Forward> > storage
            )
        {
          forwardStorage = storage;
        }

        /**
         * @brief Set internal storage object for adjoint solution
         */
        /*void setStorage(
            const std::shared_ptr<SolutionStorage<Traits,ModelTypes::Geoelectrics,Direction::Adjoint> > storage
            )
        {
          adjointStorage = storage;
        }*/

        /**
         * @brief Provide access to underlying parameter fields
         */
        void setParameterList(const std::shared_ptr<const ParameterList>& list)
        {
          sigma_bgField = (*list).get("sigma_bg");
          kappaField = (*list).get("kappa");
        }

        /**
         * @brief Provide access to measurement storage object
         */
        void setMeasurementList(const std::shared_ptr<const MeasurementList>& list)
        {
        }

        RF timestep() const
        {
          if (traits.basePotentialEvaluation)
            return traits.config().template get<RF>("time.end"); // max time, because stationary
          return traits.config().template get<RF>("time.step_ERT");
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

        // the cell with this index contains the electrode for this model
        unsigned int electrode_cell_index() const
        {
          return electrode_cell_indx;
        };
        
        // provide an index set to compare cells with the electrode cell
        auto& index_set() const
        {
          return traits.grid().leafGridView().indexSet();
        };

        template<typename Element, typename Domain, typename Time>
          RF potential(const Element& elem, const Domain& x, const Time& time) const
          {
            typename Traits::GridTraits::Scalar localConcentration;
            (*forwardStorage).value(time,elem,x,localConcentration);
            return localConcentration[0];
          }

        // electrical conductivity based on concentration transformation in each element
        template<typename Element, typename Domain, typename Time>
          RF cond (const Element& elem, const Domain& x, const Time& time) const
          {
            // read constant kappa from file; if not there try looking for field
            typename Traits::GridTraits::Scalar kappa = -5.0; // TODO there must be a better way to do this!
            kappa = traits.config().template get<RF>("parameters.kappa",kappa);
            if (kappa == -5.0)
            {
              if (!kappaField)
              {
                std::cout << "ModelParameters::cond " << this->name() << " " << this << std::endl;
                DUNE_THROW(Dune::Exception,"kappa field not set in geoelectrics model parameters");
              }
              (*kappaField).evaluate(elem,x,kappa);
            }

            if (!sigma_bgField)
            {
              std::cout << "ModelParameters::cond " << this->name() << " " << this << std::endl;
              DUNE_THROW(Dune::Exception,"sigma_bg not set in geoelectrics model parameters");
            }
            typename Traits::GridTraits::Scalar sigma_bg;
            (*sigma_bgField).evaluate(elem,x,sigma_bg);           

            // transform concentration to an electric conductivity, unless base potential evaluation
            if (traits.basePotentialEvaluation)
              return sigma_bg[0];
            return (*transportParams).concentration(elem,x,time)*kappa[0]+sigma_bg[0];
          }

        /**
         * @brief Water flux across given interface
         */
        template<typename Intersection, typename IDomain, typename Time>
          RF flux(const Intersection& is, const IDomain& x, const Time& time) const
          {
            typename Traits::GridTraits::Vector localFlux;
            const auto& global  = is.geometry().global(x);
            const auto& xInside = is.inside().geometry().local(global);
            (*forwardStorage).flux(time,is.inside(),xInside,localFlux);

            return localFlux * is.unitOuterNormal(x);
          }

        /**
         * @brief Maximum water flux on element (for CFL condition)
         */
        template<typename Element, typename Time>
          RF maxFluxNorm(const Element& elem, const Time& time) const
          {
            RF output = 0.;
            typename Traits::GridTraits::Vector localFlux;

            Dune::GeometryType type = elem.geometry().type();
            const auto& rule = Dune::QuadratureRules
              <typename Traits::GridTraits::DomainField, Traits::GridTraits::dim>::rule(type, 1);

            // loop over quadrature points 
            for(const auto& point : rule) 
            {
              const auto& x = point.position();
              (*forwardStorage).flux(time,elem,x,localFlux);
              RF norm = 0.;
              for (const auto& entry : localFlux)
                norm += entry*entry;
              norm = std::sqrt(norm);

              output = std::max(output,norm);
            }

            return output;
          }

        /**
         * @brief Adjoint source term based on measurements
         */
        /*template<typename Element, typename Domain, typename Time>
          RF adjointSource(const Element& elem, const Domain& x, const Time& t) const
          {
            return 0.;
          }*/

        /**
         * @brief Adjoint flux term based on measurements
         */
        /*template<typename Intersection, typename IDomain, typename Time>
          RF adjointFlux(const Intersection& is, const IDomain& x, const Time& t) const
          {
            return 0.;
          }*/

        /**
         * @brief Make ModelParameters of different model available
         */
        void registerModel(
            const std::string& name,
            const std::shared_ptr<ModelParameters<Traits,ModelTypes::Transport> >& otherParams
            )
        {
       	  transportParams = otherParams;  
        }
        void registerModel(
            const std::string& name,
            const std::shared_ptr<ModelParameters<Traits,ModelTypes::Moments_ERT> >& otherParams
            )
        {
        }

      };

    /**
     * @brief Equation traits class for forward geoelectrics equation
     */
    template<typename Traits, typename DirectionType>
      class EquationTraits<Traits, ModelTypes::Geoelectrics, DirectionType>
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
          // implicit linear solver for transient case
          template<typename... T>
            using StationarySolver = StationaryLinearSolver_CG_AMG_SSOR_reuse_matrix<T...>; // solver is modified for ERT matrix storage
          template<typename... T>
            using TransientSolver  = ImplicitLinearSolver<T...>;

          // use implicit euler for timestepping
          // alternative: Alexander2
          using OneStepScheme = Dune::PDELab::ImplicitEulerParameter<RangeField>;

          // use RT0 flux reconstruction
          // alternatives: BDM1 and RT1
          template<typename... T>
            using FluxReconstruction = RT0Reconstruction<T...>;

          // store complete space-time solution
          // alternative: only store last two steps
          template<typename... T>
            using StorageContainer = LastTwoContainer<T...>;

          // use next timestep when interpolating stored solution
          // alternatives: PreviousTimestep and LinearInterpolation
          template<typename... T>
            using TemporalInterpolation = NextTimestep<T...>;

        private:

          const FEM fem;
          GridFunctionSpace space;

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
     * @brief Class representing guess condition for geoelectrics
     */
    template<typename Traits>
      class InitialValue<Traits, ModelTypes::Geoelectrics>
      : public InitialValueBase<Traits,InitialValue<Traits,ModelTypes::Geoelectrics> >
      {
        public:

          template<typename GV>
            InitialValue(
                const Traits& traits,
                const GV& gv,
                const ModelParameters<Traits,ModelTypes::Geoelectrics>& parameters
                )
            : InitialValueBase<Traits,InitialValue<Traits,ModelTypes::Geoelectrics> >(gv)
            {}

          template<typename Domain, typename Value>
            void evaluateGlobal(const Domain& x, Value& y) const
            {
              // homogeneous guess conditions
              y = 0.;
            }
      };

    template<typename Traits, typename ModelType, typename DirectionType>
      class SourceTerm;

    /**
     * @brief Source term of geoelectrics equation
     */
    template<typename Traits>
      class SourceTerm<Traits, ModelTypes::Geoelectrics, Direction::Forward>
      {
        const ModelParameters<Traits,ModelTypes::Geoelectrics>& parameters;
        public:
          
          SourceTerm(const ModelParameters<Traits,ModelTypes::Geoelectrics>& parameters_)
          : parameters(parameters_)
          {}
          
          template<typename Element, typename Domain, typename Time>
            auto q (const Element& elem, const Domain& x, const Time& t) const
            {
              using RF = typename Traits::GridTraits::RangeField;
              
              // initialize electric current
              RF I = 0.0;
              
              unsigned int source_index = parameters.electrode_cell_index();   // this is the index of the cell containing the electrode
              unsigned int current_index = parameters.index_set().index(elem); // this is the index of the cell we are looking at right now
              
              // check if the current cell contains the electrode, if so -> give it a source term
              if (source_index == current_index)
              {
                I = 1.0;
              }
              return I;
            }
      };

    /**
     * @brief Source term of adjoint geoelectrics equation
     */
    /*template<typename Traits>
      class SourceTerm<Traits, ModelTypes::Geoelectrics, Direction::Adjoint>
      {
        const ModelParameters<Traits,ModelTypes::Geoelectrics>& parameters;

        public:

        SourceTerm(const ModelParameters<Traits,ModelTypes::Geoelectrics>& parameters_)
          : parameters(parameters_)
        {}

        template<typename Element, typename Domain, typename Time>
          auto q (const Element& elem, const Domain& x, const Time& t) const
          {
            // adjoint source depends on parameters / measurements
            return parameters.adjointSource(elem,x,t);
          }
      };*/

    /**
     * @brief Define geoelectrics equation as a differential equation
     */
    template<typename Traits, typename DomainType, typename DirectionType>
      class Equation<Traits,ModelTypes::Geoelectrics,DomainType,DirectionType>
      : public DifferentialEquation<Traits,ModelTypes::Geoelectrics,DomainType,DirectionType>
      {
        using DifferentialEquation<Traits,ModelTypes::Geoelectrics,DomainType,DirectionType>::DifferentialEquation;
      };

    /**
     * @brief Spatial local operator of the geoelectrics equation (CCFV version)
     */
    template<typename Traits, typename DirectionType>
      class SpatialOperator<Traits, ModelTypes::Geoelectrics, Discretization::CellCenteredFiniteVolume, DirectionType>
      : public Dune::PDELab::NumericalJacobianSkeleton
      <SpatialOperator<Traits, ModelTypes::Geoelectrics, Discretization::CellCenteredFiniteVolume, DirectionType> >,
      public Dune::PDELab::NumericalJacobianBoundary
        <SpatialOperator<Traits, ModelTypes::Geoelectrics, Discretization::CellCenteredFiniteVolume, DirectionType> >,
      public Dune::PDELab::FullSkeletonPattern, 
      public Dune::PDELab::FullVolumePattern,
      public Dune::PDELab::LocalOperatorDefaultFlags,
      public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<typename Traits::GridTraits::RangeField>
      {
        using GridTraits = typename Traits::GridTraits;

        enum {dim = GridTraits::dim};

        using RF      = typename GridTraits::RangeField;
        using DF      = typename GridTraits::DomainField;
        using Domain  = typename GridTraits::Domain;
        using IDomain = typename GridTraits::IDomain;

        public:

        // pattern assembly flags
        enum {doPatternVolume   = true};
        enum {doPatternSkeleton = true};

        // residual assembly flags
        enum {doAlphaSkeleton  = true};
        enum {doAlphaBoundary  = true};
        enum {doLambdaVolume   = true};
        //enum {doLambdaSkeleton = DirectionType::isAdjoint()};
        enum {doLambdaSkeleton = false};
        enum {doLambdaBoundary = true};

        private:

        const ModelParameters<Traits,ModelTypes::Geoelectrics>&              parameters;
        const Boundary       <Traits,ModelTypes::Geoelectrics,DirectionType> boundary;
        const SourceTerm     <Traits,ModelTypes::Geoelectrics,DirectionType> sourceTerm;

        RF time;

        public:

        /**
         * @brief Constructor
         */
        SpatialOperator(
            const Traits& traits,
            const ModelParameters<Traits,ModelTypes::Geoelectrics>& parameters_
            )
          : parameters(parameters_), boundary(traits,parameters.name()), sourceTerm(parameters)
        {}

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
            const RF normalFlux = skeletonNormalFlux(ig.intersection(),x_s(lfsu_s,0),x_n(lfsu_n,0),time);
            const RF faceVolume = ig.geometry().volume();

            r_s.accumulate(lfsv_s,0,   normalFlux * faceVolume);
            r_n.accumulate(lfsv_n,0, - normalFlux * faceVolume);
          }

        /**
         * @brief Skeleton integral depending on test and ansatz functions
         */
        template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
          void alpha_boundary (const IG& ig, 
              const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
              R& r_s) const
          {
            const IDomain& faceCenterLocal = referenceElement(ig.geometry()).position(0,0);

            // evaluate boundary condition type
            Dune::Modelling::BoundaryCondition::Type bc = boundary.bc(ig.intersection(),faceCenterLocal,time);
            if (!Dune::Modelling::BoundaryCondition::isDirichlet(bc))
              return;

            // Dirichlet boundary conditions
            const RF normalFlux = dirichletNormalFlux(ig.intersection(),x_s(lfsu_s,0),time);
            const RF faceVolume = ig.geometry().volume();

            r_s.accumulate(lfsv_s,0, normalFlux * faceVolume);
          }

        /**
         * @brief Volume integral depending only on test functions (source term)
         */
        template<typename EG, typename LFSV, typename R>
          void lambda_volume (const EG& eg, const LFSV& lfsv, R& r) const
          {
            const Domain& cellCenterLocal = referenceElement(eg.geometry()).position(0,0);
            r.accumulate(lfsv,0, - sourceTerm.q(eg.entity(),cellCenterLocal,time));

            /*if (DirectionType::isAdjoint())
            {
              RF q = 0.;

              Dune::GeometryType geometrytype = eg.geometry().type();
              /// @todo order
              const Dune::QuadratureRule<DF, dim>& rule = Dune::QuadratureRules<DF, dim>::rule(geometrytype, 2);

              // loop over quadrature points 
              for(const auto& point : rule) 
              {
                const Domain& x = point.position();
                const RF factor = point.weight() * eg.geometry().integrationElement(x);

                q += sourceTerm.q(eg.entity(),x,time) * factor;
              }

              r.accumulate(lfsv,0, -q);
            }*/
          }

        /**
         * @brief Skeleton integral independent of ansatz functions
         */
        // each face is only visited ONCE!
        /*template<typename IG, typename LFSV, typename R>
          void lambda_skeleton (const IG& ig, const LFSV& lfsv_s, const LFSV& lfsv_n, R& r_s, R& r_n) const
          {
            RF q = 0.;

            Dune::GeometryType geometrytype = ig.geometry().type();
            /// @todo order
            const Dune::QuadratureRule<DF, dim-1>& rule = Dune::QuadratureRules<DF, dim-1>::rule(geometrytype, 1);

            // loop over quadrature points 
            for(const auto& point : rule) 
            {
              const IDomain& x = point.position();
              const RF factor  = point.weight() * ig.geometry().integrationElement(x);

              q += parameters.adjointFlux(ig.intersection(),x,time) * factor;
            }

            r_s.accumulate(lfsv_s,0, -q);
            r_n.accumulate(lfsv_n,0,  q);
          }*/

        /**
         * @brief Boundary integral independent of ansatz functions
         */
        template<typename IG, typename LFSV, typename R>
          void lambda_boundary (const IG& ig, const LFSV& lfsv, R& r_s) const
          {
            const IDomain& faceCenterLocal = referenceElement(ig.geometry()).position(0,0);

            // evaluate boundary condition type
            Dune::Modelling::BoundaryCondition::Type bc = boundary.bc(ig.intersection(),faceCenterLocal,time);
            if (!Dune::Modelling::BoundaryCondition::isNeumann(bc))
              return;

            // Neumann boundary conditions
            const RF faceVolume = ig.geometry().volume();
            const RF j = boundary.j(ig.intersection(),faceCenterLocal,time);

            r_s.accumulate(lfsv,0, j * faceVolume);
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
                "derivative of normal flux not implemented for geoelectrics flow equation");
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
                const RF K_inside = parameters.cond(is.inside(),cellCenterInside,time);
                // revert flux computation to obtain gradient on boundary
                return - boundaryFlux / K_inside;
              }
              else
                DUNE_THROW(Dune::Exception,"unknown boundary condition type");
            }
            else
              return 0.;
          }

        private:

        /**
         * @brief Flux in normal direction across interface
         */
        template<typename Intersection>
          RF skeletonNormalFlux(const Intersection& is, RF innerValue, RF outerValue, RF time) const
          {
            // geometry information
            const Domain& cellCenterInside  = referenceElement(is.inside() .geometry()).position(0,0);
            const Domain& cellCenterOutside = referenceElement(is.outside().geometry()).position(0,0);

            // conductivity
            const RF K_inside  = parameters.cond(is.inside(), cellCenterInside, time);
            const RF K_outside = parameters.cond(is.outside(),cellCenterOutside,time);
            const RF K = havg(K_inside,K_outside);

            return - K * skeletonNormalDerivative(is,innerValue,outerValue);
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
            const Domain&  cellCenterInside = referenceElement(is.inside().geometry()).position(0,0);

            // conductivity
            const RF K_inside = parameters.cond(is.inside(),cellCenterInside,time);

            return - K_inside * dirichletNormalDerivative(is,innerValue,time);
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
            const T eps = 1e-30;
            return 2./(1./(a + eps) + 1./(b + eps));
          }

      };
  }
}

#endif // DUNE_FTG_GEOELECTRICS_HH
