// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FTG_GROUNDWATER_HH
#define DUNE_FTG_GROUNDWATER_HH

#include<dune/pdelab/common/referenceelements.hh>
#include<dune/pdelab/finiteelementmap/p0fem.hh>
#include<dune/pdelab/constraints/p0.hh>
#include<dune/pdelab/backend/istl.hh>

#include<dune/modelling/fluxreconstruction.hh>
#include<dune/modelling/solutionstorage.hh>
#include<dune/ftg/override/solvers.hh>

namespace Dune {
  namespace Modelling {
    /**
     * @brief Parameter class for the groundwater flow equation
     */
    template<typename Traits>
      class ModelParameters<Traits, ModelTypes::Groundwater>
      : public ModelParametersBase<Traits>
      {
        using RF = typename Traits::GridTraits::RangeField;
        using DF = typename Traits::GridTraits::DomainField;
        using ParameterList   = typename Traits::ParameterList;
        using ParameterField  = typename ParameterList::SubRandomField;
        using MeasurementList = typename Traits::MeasurementList;

        const Traits& traits;

        std::shared_ptr<SolutionStorage<Traits,ModelTypes::Groundwater,Direction::Forward> > forwardStorage;

        std::shared_ptr<ParameterField> conductivityField;

        public:
          unsigned int model_number = 0; // we need this, as the number is requested for ERT models

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
            const std::shared_ptr<SolutionStorage<Traits,ModelTypes::Groundwater,Direction::Forward> > storage
            )
        {
          forwardStorage = storage;
        }

        /**
         * @brief Provide access to underlying parameter fields
         */
        void setParameterList(const std::shared_ptr<const ParameterList>& list)
        {
          conductivityField = (*list).get("conductivity");
        }

        /**
         * @brief Provide access to measurement storage object
         */
        void setMeasurementList(const std::shared_ptr<const MeasurementList>& list)
        {}

        RF timestep() const
        {
          return traits.config().template get<RF>("time.end");
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
         * @brief Conductivity value at given position
         */
        template<typename Element, typename Domain, typename Time>
          RF cond(const Element& elem, const Domain& x, const Time& t) const
          {
            
            if (!conductivityField)
            {
              std::cout << "ModelParameters::cond " << this->name() << " " << this << std::endl;
              DUNE_THROW(Dune::Exception,"conductivity field not set in groundwater model parameters");
            }

            unsigned int current_index = index_set().index(elem);
            std::map<unsigned int,std::pair<RF,bool>> wellcells = well_cells(); // this step is necessary! direct usage of well_cells() leads to mismatch!
            typename std::map<unsigned int,std::pair<RF,bool>>::iterator temp = wellcells.find(current_index);

            if ( temp != wellcells.end() ) 
            {
              // a well cell -> return very high conductivity to achieve "shortcut" or very low conductivity to achieve "packer"
              typename Traits::GridTraits::Scalar value;
              value = 10.0;
              if ( (temp->second).first == 0.0)
                value = 0.0;
              return value[0];
            }
            // no well cell
            typename Traits::GridTraits::Scalar value;
            (*conductivityField).evaluate(elem,x,value);
            return value[0];
          }

        /**
         * @brief Conductivity value at given position
         */
        template<typename IS, typename Time>
          RF cond_intersection(const IS& is, const Time& t) const
          {
            
            if (!conductivityField)
            {
              std::cout << "ModelParameters::cond " << this->name() << " " << this << std::endl;
              DUNE_THROW(Dune::Exception,"conductivity field not set in groundwater model parameters");
            }

            std::vector<int> cell_status; // 0: normal cell; 1: well cell; -1 packer cell
            for (unsigned int i =0; i<2;i++)
            {
              auto elem = is.inside();
              if (i == 1)
                elem = is.outside();

              unsigned int current_index = index_set().index(elem);
              std::map<unsigned int,std::pair<RF,bool>> wellcells = well_cells(); // this step is necessary!
              typename std::map<unsigned int,std::pair<RF,bool>>::iterator temp = wellcells.find(current_index);

              if ( temp != wellcells.end() ) 
              {
                if ( (temp->second).first == 0.0)
                  cell_status.push_back(-1);
                else
                  cell_status.push_back(1);
              } else {cell_status.push_back(0);}
            }
            
            if ((cell_status[0]==-1) || (cell_status[1]==-1) )
            {
              typename Traits::GridTraits::Scalar value = 0.0;
              return value[0];  // if one cell is packer: K=0
            }
            if ((cell_status[0]==1) && (cell_status[1]==1))
            {
              typename Traits::GridTraits::Scalar value = 1.0;
              return value[0];  // of both cells are well cells: K=1
            }

            // all other cases: K = havg(K1,K2)
            std::vector<typename Traits::GridTraits::Scalar> values;                          
            for (unsigned int i =0; i<2;i++)
            {
              auto elem = is.inside();
              if (i == 1)
                elem = is.outside();

              auto x = referenceElement(elem.geometry()).position(0,0);
              typename Traits::GridTraits::Scalar tmp;
              (*conductivityField).evaluate(elem,x,tmp);
              values.push_back(tmp);
            }
            const typename Traits::GridTraits::Scalar eps = 1e-30;
            return 2./(1./(values[0] + eps) + 1./(values[1] + eps));
          }

        template<typename Element, typename Domain>
          RF head(const Element& elem, const Domain& x) const
          {
            typename Traits::GridTraits::Scalar localHead;
            const auto& global  = elem.geometry().global(x);
            const auto& xInside = elem.geometry().local(global);
            (*forwardStorage).value(traits.config().template get<RF>("time.end"),elem,xInside,localHead);
            return localHead;
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

          // this does not work correctly! use the other CFL implementation...
          template<typename Element, typename Time>
          RF dt_CFL(const Element& elem, const Time& time) const
          {
            RF output = 0.;
            typename Traits::GridTraits::Vector localFlux;

            Dune::GeometryType type = elem.geometry().type();
            const auto& rule = Dune::QuadratureRules
              <typename Traits::GridTraits::DomainField, Traits::GridTraits::dim>::rule(type, 1);

            std::vector<DF> dx;
            dx.push_back(elem.geometry().corner(1)[0]-elem.geometry().corner(0)[0]); // delta_x
            dx.push_back(elem.geometry().corner(2)[1]-elem.geometry().corner(0)[1]); // delta_y
            dx.push_back(elem.geometry().corner(4)[2]-elem.geometry().corner(0)[2]); // delta_z

            // loop over quadrature points 
            for(const auto& point : rule) 
            {
              const auto& x = point.position();
              (*forwardStorage).flux(time,elem,x,localFlux);
              RF norm = 0.;
              for (unsigned int i=0; i<localFlux.size(); i++)
              {
                norm += std::abs(localFlux[i])/dx[i];
              }
              norm = 1/norm;

              output = std::max(output,norm);
            }

            return output;
          }

        /**
         * @brief Make ModelParameters of different model available && make this model available for the others!
         */
        void registerModel(
            const std::string& name,
            const std::shared_ptr<ModelParameters<Traits,ModelTypes::Transport> >& otherParams
            )
        {
        }

	    void registerModel(
            const std::string& name,
            const std::shared_ptr<ModelParameters<Traits,ModelTypes::Moments_c> >& otherParams
            )
        {
        }

	    void registerModel(
            const std::string& name,
            const std::shared_ptr<ModelParameters<Traits,ModelTypes::Moments_ERT> >& otherParams
            )
        {
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
     * @brief Equation traits class for forward groundwater flow equation
     */
    template<typename Traits, typename DirectionType>
      class EquationTraits<Traits, ModelTypes::Groundwater, DirectionType>
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
            using StationarySolver = StationaryLinearSolver_CG_AMG_SSOR<T...>;
          template<typename... T>
            using TransientSolver  = ImplicitLinearSolver<T...>;
          using OneStepScheme = Dune::PDELab::ImplicitEulerParameter<RangeField>;
          template<typename... T>
            using FluxReconstruction = RT0Reconstruction<T...>;
          template<typename... T>
            using StorageContainer = LastTwoContainer<T...>;
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
     * @brief Class representing guess condition for hydraulic head
     */
    template<typename Traits>
      class InitialValue<Traits, ModelTypes::Groundwater>
      : public InitialValueBase<Traits,InitialValue<Traits,ModelTypes::Groundwater> >
      {
        public:

          template<typename GV>
            InitialValue(
                const Traits& traits,
                const GV& gv,
                const ModelParameters<Traits,ModelTypes::Groundwater>& parameters
                )
            : InitialValueBase<Traits,InitialValue<Traits,ModelTypes::Groundwater> >(gv)
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
     * @brief Source term of groundwater flow equation
     */
    template<typename Traits>
      class SourceTerm<Traits, ModelTypes::Groundwater, Direction::Forward>
      {
          const ModelParameters<Traits,ModelTypes::Groundwater>& parameters;
          using RF  = typename Traits::GridTraits::RangeField;
          std::map<unsigned int, std::pair<RF, bool>> well_cells;
        public:

          SourceTerm(const ModelParameters<Traits,ModelTypes::Groundwater>& parameters_)
          : parameters(parameters_)
          {
            well_cells = parameters.well_cells();
          }

          template<typename Element, typename Domain, typename Time>
            auto q (const Element& elem, const Domain& x, const Time& t) const
            {
              unsigned int current_index = parameters.index_set().index(elem);
              auto temp = well_cells.find(current_index);

              if ( !(temp->first == well_cells.end()->first) ) 
              {
                return (temp->second).first; // a well cell, return extraction rate
              }
              return 0.0; // no well cell
            }
      };

    /**
     * @brief Define groundwater equation as a differential equation
     */
    template<typename Traits, typename DomainType, typename DirectionType>
      class Equation<Traits,ModelTypes::Groundwater,DomainType,DirectionType>
      : public DifferentialEquation<Traits,ModelTypes::Groundwater,DomainType,DirectionType>
      {
        using DifferentialEquation<Traits,ModelTypes::Groundwater,DomainType,DirectionType>::DifferentialEquation;
      };

    /**
     * @brief Spatial local operator of the groundwater flow equation (CCFV version)
     */
    template<typename Traits, typename DirectionType>
      class SpatialOperator<Traits, ModelTypes::Groundwater, Discretization::CellCenteredFiniteVolume, DirectionType>
      : public Dune::PDELab::NumericalJacobianSkeleton
        <SpatialOperator<Traits, ModelTypes::Groundwater, Discretization::CellCenteredFiniteVolume, DirectionType> >,
      public Dune::PDELab::NumericalJacobianBoundary
        <SpatialOperator<Traits, ModelTypes::Groundwater, Discretization::CellCenteredFiniteVolume, DirectionType> >,
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
        enum {doLambdaSkeleton = false};
        enum {doLambdaBoundary = true};

        private:

        const ModelParameters<Traits,ModelTypes::Groundwater>&              parameters;
        const Boundary       <Traits,ModelTypes::Groundwater,DirectionType> boundary;
        const SourceTerm     <Traits,ModelTypes::Groundwater,DirectionType> sourceTerm;

        RF time;

        public:

        /**
         * @brief Constructor
         */
        SpatialOperator(
            const Traits& traits,
            const ModelParameters<Traits,ModelTypes::Groundwater>& parameters_
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
          }

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
            //const Domain& cellCenterInside  = referenceElement(is.inside() .geometry()).position(0,0);
            //const Domain& cellCenterOutside = referenceElement(is.outside().geometry()).position(0,0);

            // conductivity
            //const RF K_inside  = parameters.cond(is.inside(), cellCenterInside, time);
            //const RF K_outside = parameters.cond(is.outside(),cellCenterOutside,time);
            //const RF K = havg(K_inside,K_outside);
            const RF K = parameters.cond_intersection(is,time);

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

#endif // DUNE_FTG_GROUNDWATER_HH
