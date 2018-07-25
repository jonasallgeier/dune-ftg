// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FTG_MOMENTS_ERT_HH
#define DUNE_FTG_MOMENTS_ERT_HH

#include<dune/pdelab/common/referenceelements.hh>
#include<dune/pdelab/finiteelementmap/p0fem.hh>
#include<dune/pdelab/constraints/p0.hh>
#include<dune/pdelab/backend/istl.hh>

#include<dune/modelling/fluxreconstruction.hh>
#include<dune/modelling/solutionstorage.hh>
#include<dune/ftg/override/solvers_moments_ERT.hh>

namespace Dune {
  namespace Modelling {

    /**
     * @brief Parameter class for the solute transport equation
     */
    template<typename Traits>
      class ModelParameters<Traits, typename ModelTypes::Moments_ERT>
      : public ModelParametersBase<Traits>
      {
        using RF = typename Traits::GridTraits::RangeField;

        using ParameterList  = typename Traits::ParameterList;
        using ParameterField = typename ParameterList::SubRandomField;

        const Traits& traits;

        std::shared_ptr<SolutionStorage<Traits,ModelTypes::Moments_ERT,Direction::Forward> > forwardStorage;

        std::shared_ptr<ParameterField> sigma_bgField;
        std::shared_ptr<ParameterField> kappaField;

        std::shared_ptr<const ModelParameters<Traits,ModelTypes::Moments_c> > c_momentParams;
        std::shared_ptr<const ModelParameters<Traits,ModelTypes::Geoelectrics> > geoelectricsParams;
        std::shared_ptr<const ModelParameters<Traits,ModelTypes::Groundwater> > groundwaterParams;

        public:
          unsigned int model_number;
          unsigned int k;

        ModelParameters(const Traits& traits_, const std::string& name)
          : ModelParametersBase<Traits>(name), traits(traits_)
        {
          std::string str = name;
          std::string common_base = "momentsERT_";
          auto start_position_to_erase = str.find(common_base);
          str.erase(start_position_to_erase, common_base.size());      
          std::replace(str.begin(), str.end(), '_', ' ');


          std::stringstream ss; 
          ss << str;
          std::string temp;

          ss >> temp; // get the ERT model number
          std::stringstream(temp) >> model_number;
          ss >> temp; // get the moment number k
          std::stringstream(temp) >> k;
        }

        /**
         * @brief Model parameters should exist only once per (named) model
         */
        ModelParameters(const ModelParameters& other) = delete;

        /**
         * @brief Set internal storage object for forward solution
         */
        void setStorage(
            const std::shared_ptr<SolutionStorage<Traits,ModelTypes::Moments_ERT,Direction::Forward> > storage
            )
        {
          forwardStorage = storage;
        }

        /**
         * @brief Provide access to underlying parameter fields
         */
        void setParameterList(const std::shared_ptr<const typename Traits::ParameterList>& list)
        {
          sigma_bgField = (*list).get("sigma_bg");
          kappaField = (*list).get("kappa");
        }

        /**
         * @brief Provide access to measurement storage object
         */
        void setMeasurementList(const std::shared_ptr<const typename Traits::MeasurementList>& list)
        {
        }

        RF timestep() const
        {
          return traits.config().template get<RF>("time.end"); // use max time, because stationary
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

        template<typename Element, typename Domain, typename Time>
          RF sigma (const Element& elem, const Domain& x, const Time& time) const
          {
            if (!sigma_bgField)
            {
              std::cout << "ModelParameters::sigma " << this->name() << " " << this << std::endl;
              DUNE_THROW(Dune::Exception,"sigma_bg field not set in geoelectrics model parameters");
            }
            typename Traits::GridTraits::Scalar sigma_bg;
            (*sigma_bgField).evaluate(elem,x,sigma_bg);
            return sigma_bg[0];
          }

        template<typename Element, typename Domain, typename Time>
          RF kappa (const Element& elem, const Domain& x, const Time& time) const
          {
            // read constant kappa from file; if not there try looking for field
            typename Traits::GridTraits::Scalar kappa = -5.0;
            kappa = traits.config().template get<RF>("parameters.kappa",kappa);
            if (kappa == -5.0)
            {
              if (!kappaField)
              {
                std::cout << "ModelParameters::kappa " << this->name() << " " << this << std::endl;
                DUNE_THROW(Dune::Exception,"kappa field not set in geoelectrics model parameters");
              }
              (*kappaField).evaluate(elem,x,kappa);
            }
            return kappa[0];
          }

        template<typename Element, typename Domain, typename Time>
          RF c_moment(const Element& elem, const Domain& x, const Time& time) const
          {
            if (!c_momentParams)
              DUNE_THROW(Dune::Exception,"moment k-1 model parameters not set in moment k model parameters");
            
            return (*c_momentParams).moment(elem,x,time);
          }

        template<typename Intersection, typename IDomain, typename Time>
          RF vNormal (const Intersection& is, const IDomain& x, const Time& t) const
          {
            if (!groundwaterParams)
              DUNE_THROW(Dune::Exception,"groundwater model parameters not set in transport model parameters");
            
            return (*groundwaterParams).flux(is,x,t);
          }

        template<typename Element, typename Domain, typename Time>
          RF el_potential(const Element& elem, const Domain& x, const Time& time) const
          {
            if (!geoelectricsParams)
              DUNE_THROW(Dune::Exception,"geoelectrics parameters not set in ERT moment model parameters");
            
            return (*geoelectricsParams).potential(elem,x,time);
          }

        template<typename Element, typename Domain, typename Time>
          auto potential_gradient(const Element& elem, const Domain& x, const Time& time) const
          {
            if (!geoelectricsParams)
              DUNE_THROW(Dune::Exception,"geoelectrics parameters not set in ERT moment model parameters");
            
            return (*geoelectricsParams).potential_gradient(elem,x,time);
          }


        /**
         * @brief Make ModelParameters of different model available
         */
        void registerModel(
            const std::string& name, 
            const std::shared_ptr<ModelParameters<Traits,ModelTypes::Moments_c> >& otherParams
            )
        {
          c_momentParams = otherParams;
        }

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
          geoelectricsParams = otherParams;
        }

        auto& index_set() const
        {
          return traits.grid().leafGridView().indexSet();
        };

        std::map<unsigned int, std::pair<RF, bool>> well_cells() const
        {
          return traits.read_well_cells();
        };
      };

    /**
     * @brief Equation traits class for forward moments ERT
     */
    template<typename Traits, typename DirectionType>
      class EquationTraits<Traits, typename ModelTypes::Moments_ERT, DirectionType>
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
            using StationarySolver = StationarySolverMomentsERT_CG_AMG_SSOR_reuse_matrix<T...>;
            //using StationarySolver = StationaryLinearSolver_BCGS_AMG_ILU0<T...>;
          template<typename... T>
            using TransientSolver  = TransientSolverMomentsERT<T...>; // dummy implementation
            using OneStepScheme = Dune::PDELab::ImplicitEulerParameter<RangeField>;
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
      class InitialValue<Traits, typename ModelTypes::Moments_ERT>
      : public InitialValueBase<Traits,InitialValue<Traits,ModelTypes::Moments_ERT> >
      {
        public:
          template<typename GV>
            InitialValue(
                const Traits& traits_,
                const GV& gv,
                const ModelParameters<Traits,ModelTypes::Moments_ERT>& parameters_
                )
            : InitialValueBase<Traits,InitialValue<Traits,ModelTypes::Moments_ERT> >(gv)
            {}

          template<typename Domain, typename Value>
            void evaluateGlobal(const Domain& x, Value& y) const
            {
              // homogeneous initial/guess conditions
              y= 0.0;
            }
      };

    /**
     * @brief Define transport equation as a differential equation
     */
    template<typename Traits, typename DomainType, typename DirectionType>
      class Equation<Traits,ModelTypes::Moments_ERT,DomainType,DirectionType>
      : public DifferentialEquation<Traits,ModelTypes::Moments_ERT,DomainType,DirectionType>
      {
        using DifferentialEquation<Traits,ModelTypes::Moments_ERT,DomainType,DirectionType>::DifferentialEquation;
      };

    /**
     * @brief Spatial local operator of the convection-diffusion equation (CCFV version)
     */
    template<typename Traits, typename DirectionType>
      class SpatialOperator<Traits, typename ModelTypes::Moments_ERT, Discretization::CellCenteredFiniteVolume, DirectionType>
      : public Dune::PDELab::NumericalJacobianSkeleton<SpatialOperator<Traits, ModelTypes::Moments_ERT, Discretization::CellCenteredFiniteVolume, DirectionType> >,
      public Dune::PDELab::NumericalJacobianBoundary<SpatialOperator<Traits, ModelTypes::Moments_ERT, Discretization::CellCenteredFiniteVolume, DirectionType> >,
      public Dune::PDELab::NumericalJacobianVolume<SpatialOperator<Traits, ModelTypes::Moments_ERT, Discretization::CellCenteredFiniteVolume, DirectionType> >,
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
        enum {doAlphaVolume    = false};
        enum {doAlphaSkeleton  = true};
        enum {doAlphaBoundary  = true};
        enum {doLambdaVolume   = false};
        enum {doLambdaSkeleton = true};
        enum {doLambdaBoundary = true};

        private:

        const Traits& traits;

        const ModelParameters<Traits,ModelTypes::Moments_ERT>&              parameters;
        const Boundary       <Traits,ModelTypes::Moments_ERT,DirectionType> boundary;

        typename Traits::GridTraits::Grid::LeafGridView  lgv;
        std::map<unsigned int, std::pair<RF, bool>> well_cells;

        RF time;

        public:

        /**
         * @brief Constructor
         */
        SpatialOperator(
            const Traits& traits_,
            const ModelParameters<Traits,ModelTypes::Moments_ERT>& parameters_
            )
          : traits(traits_), parameters(parameters_), boundary(traits,parameters.name()), lgv(traits_.grid().leafGridView()) 
        {
          well_cells = parameters.well_cells();
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
            const RF faceVolume = ig.geometry().volume();
            const Domain& cellCenterInside  = referenceElement(ig.intersection().inside() .geometry()).position(0,0);
            const Domain& cellCenterOutside = referenceElement(ig.intersection().outside().geometry()).position(0,0);

            const RF sigma_inside  = parameters.sigma(ig.intersection().inside(),cellCenterInside,time);
            const RF sigma_outside = parameters.sigma(ig.intersection().outside(),cellCenterOutside,time);
            const RF sigma_0 = havg(sigma_inside,sigma_outside);

            RF left_part = - sigma_0 * skeletonNormalDerivative(ig.intersection(),x_s(lfsu_s,0),x_n(lfsu_n,0)) * faceVolume;

            r_s.accumulate(lfsv_s,0,  left_part);
            r_n.accumulate(lfsv_n,0, -left_part);
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
            // evaluate boundary condition type
            const IDomain& faceCenterLocal = referenceElement(ig.geometry()).position(0,0);
            Dune::Modelling::BoundaryCondition::Type bc = boundary.bc(ig.intersection(),faceCenterLocal,time);
            if (!Dune::Modelling::BoundaryCondition::isDirichlet(bc))
              return;

            // Dirichlet boundary conditions
            const RF faceVolume = ig.geometry().volume();
            const Domain&  cellCenterInside = referenceElement(ig.intersection().inside().geometry()).position(0,0);

            const RF sigma_0 = parameters.sigma(ig.intersection().inside(),cellCenterInside,time);
            
            RF left_part = - sigma_0 * dirichletNormalDerivative(ig.intersection(),x_s(lfsu_s,0),time)*faceVolume;

            r_s.accumulate(lfsv_s,0, left_part);
          }

        template<typename IG, typename LFSV, typename R>
          void lambda_boundary (const IG& ig, const LFSV& lfsv, R& r_s) const
          {
            // evaluate boundary condition type
            const IDomain& faceCenterLocal = referenceElement(ig.geometry()).position(0,0);
            Dune::Modelling::BoundaryCondition::Type bc = boundary.bc(ig.intersection(),faceCenterLocal,time);

            // Neumann boundary condition
            if (Dune::Modelling::BoundaryCondition::isNeumann(bc))
            {
              const RF j = boundary.j(ig.intersection(),faceCenterLocal,time);
              const RF faceVolume = ig.geometry().volume();
              r_s.accumulate(lfsv,0, j * faceVolume);
              return;
            }
            if (Dune::Modelling::BoundaryCondition::isDirichlet(bc))
            {
              const RF faceVolume = ig.geometry().volume();
              const Domain& cellCenterInside  = referenceElement(ig.intersection().inside() .geometry()).position(0,0);
              
              const RF c_moment_inside  = parameters.c_moment(ig.intersection().inside(),cellCenterInside,time);
              //const RF c_moment_outside = parameters.c_moment(ig.intersection().outside(),cellCenterOutside,time);
              const RF v = parameters.vNormal(ig.intersection(),faceCenterLocal,time);
              const RF c_moment = (v >= 0) ? c_moment_inside : 0; // upwinding
              //const RF c_moment = c_moment_inside;
              const RF kappa  = parameters.kappa(ig.intersection().inside(),cellCenterInside,time);

              const Domain& faceCenterGlobal = ig.intersection().geometry().global(faceCenterLocal);
              const auto& faceCenterLocal = ig.intersection().inside().geometry().local(faceCenterGlobal);

              auto gradient_potential = parameters.potential_gradient(ig.intersection().inside(),faceCenterLocal,time);

              auto gradient_normal = ig.centerUnitOuterNormal()*gradient_potential;
              RF right_part = - kappa * c_moment * gradient_normal * faceVolume;
              r_s.accumulate(lfsv,0, right_part);
              return;
            }
          }

        template<typename IG, typename LFSV, typename R>
          void lambda_skeleton (const IG& ig, const LFSV& lfsv_s, const LFSV& lfsv_n, R& r_s, R& r_n) const
          {
            const RF faceVolume = ig.geometry().volume();
            const IDomain& faceCenterLocal = referenceElement(ig.geometry()).position(0,0);
            const Domain& cellCenterInside  = referenceElement(ig.intersection().inside() .geometry()).position(0,0);
            const Domain& cellCenterOutside = referenceElement(ig.intersection().outside().geometry()).position(0,0);
            
            const RF c_moment_inside  = parameters.c_moment(ig.intersection().inside(),cellCenterInside,time);
            const RF c_moment_outside = parameters.c_moment(ig.intersection().outside(),cellCenterOutside,time);
            const RF v = parameters.vNormal(ig.intersection(),faceCenterLocal,time);
            const RF c_moment = (v >= 0) ? c_moment_inside : c_moment_outside; // upwinding
            //const RF c_moment = havg(c_moment_inside,c_moment_outside);

            const RF kappa_inside  = parameters.kappa(ig.intersection().inside(),cellCenterInside,time);
            const RF kappa_outside = parameters.kappa(ig.intersection().outside(),cellCenterOutside,time);
            const RF kappa = havg(kappa_inside,kappa_outside);     

            const RF potential_inside  = parameters.el_potential(ig.intersection().inside(),cellCenterInside,time);
            const RF potential_outside = parameters.el_potential(ig.intersection().outside(),cellCenterOutside,time);

            RF right_part = - kappa * c_moment * skeletonNormalDerivative(ig.intersection(),potential_inside,potential_outside) * faceVolume;

            r_s.accumulate(lfsv_s,0,  right_part);
            r_n.accumulate(lfsv_n,0, -right_part);
          }


        /**
         * @brief Set time for subsequent evaluation
         */
        void setTime (const RF t)
        {
          time = t;
        }

        private:

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
  }
}

#endif // DUNE_FTG_MOMENTS_ERT_HH
