#ifndef DUNE_MODELLING_EQUATION_HH
#define DUNE_MODELLING_EQUATION_HH

#include<dune/common/timer.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>

#include<dune/modelling/solutionstorage.hh>

namespace Dune {
  namespace Modelling {

    /**
     * @brief Base class for Equations
     */
    template<typename Traits>
      class EquationBase
      {
        using ParameterList   = typename Traits::ParameterList;
        using MeasurementList = typename Traits::MeasurementList;
        using Measurements    = typename MeasurementList::SubMeasurements;

        public:

        /**
         * @brief Initialize equation to be able to solve again
         */
        virtual void initialize(
            const std::shared_ptr<const ParameterList>& parameterList,
            const std::shared_ptr<const MeasurementList>& measurementList,
            const std::shared_ptr<Measurements>& measurements
            ) = 0;

        /**
         * @brief Compute a time step
         */
        virtual void step(typename Traits::GridTraits::RangeField timestep) = 0;

        /**
         * @brief Suggest a time step to take
         */
        virtual typename Traits::GridTraits::RangeField suggestTimestep() const = 0;

        /**
         * @brief Return whether equation is solved
         */
        virtual bool finished() const = 0;
      };

    template<typename Traits>
      class DifferentialEquationBase
      : public EquationBase<Traits>
      {
        using ParameterList   = typename Traits::ParameterList;
        using MeasurementList = typename Traits::MeasurementList;
        using Measurements    = typename MeasurementList::SubMeasurements;

        public:

        /**
         *  @brief Initialize all objects and variables to compute a new transient solution
         */
        virtual void initialize(
            const std::shared_ptr<const ParameterList>& parameterList,
            const std::shared_ptr<const MeasurementList>& measurementList,
            const std::shared_ptr<Measurements>& measurements
            ) = 0;

        /**
         * @brief Return the time used in the last step
         */
        virtual typename Traits::GridTraits::RangeField getTime() const = 0;

        /**
         * @brief Return the timestep used in the last step
         */
        virtual typename Traits::GridTraits::RangeField getTimestep() const = 0;
      };

    /**
     * @brief Class for the solution of a transient PDE
     */
    template<typename Traits, typename ModelType, typename FormulationType, typename DirectionType>
      class DifferentialEquation
      : public DifferentialEquationBase<Traits>
      {
        public:

          using GridTraits = typename Traits::GridTraits;

          using Grid         = typename GridTraits::Grid;
          using RF           = typename GridTraits::RangeField;
          using Scalar       = typename GridTraits::Scalar;
          using Vector       = typename GridTraits::Vector;
          using Domain       = typename GridTraits::Domain;
          using IDomain      = typename GridTraits::IDomain;
          using Element      = typename GridTraits::Element;
          using Intersection = typename GridTraits::Intersection;

          using ParameterList   = typename Traits::ParameterList;
          using MeasurementList = typename Traits::MeasurementList;
          using Measurements    = typename MeasurementList::SubMeasurements;

          using GFS        = typename EquationTraits<Traits,ModelType,DirectionType>::GridFunctionSpace;
          using GridVector = typename EquationTraits<Traits,ModelType,DirectionType>::GridVector;
          using C          = typename GFS::template ConstraintsContainer<RF>::Type;

          // extract stationary solver from equation traits
          template<typename... T>
            using StationarySolver = typename EquationTraits<Traits,ModelType,DirectionType>
            ::template StationarySolver<T...>;
          
          // extract transient solver from equation traits
          template<typename... T>
            using TransientSolver  = typename EquationTraits<Traits,ModelType,DirectionType>
            ::template TransientSolver<T...>;
          
          // choose between stationary and transient solver based on formulation
          template<typename... T>
            using Solver = typename std::conditional
            <std::is_same<FormulationType, typename Formulation::Stationary>::value,
            StationarySolver<T...>,TransientSolver<T...> >::type;

        private:

          const Traits&                                        traits;
          const EquationTraits<Traits,ModelType,DirectionType> equationTraits;
          ModelParameters<Traits,ModelType>&                   parameters;

          RF startTime, endTime, maxTimestep, time, minTimestep, oldTime, usedTimestep;

          GridVector oldSolution, newSolution;
          C cg;

          std::shared_ptr<Solver<Traits,ModelType,DirectionType> >          solver;
          std::shared_ptr<SolutionStorage<Traits,ModelType,DirectionType> > storage;
          std::shared_ptr<Measurements>                                     measurements;

          mutable Dune::Timer stepTimer, printTimer;

        public:

          /**
           * @brief Constructor
           */
          DifferentialEquation(const Traits& traits_, ModelParameters<Traits,ModelType>& parameters_)
            : traits(traits_), equationTraits(traits), parameters(parameters_),
            oldSolution(equationTraits.gfs(),0.), newSolution(oldSolution),
            storage(new SolutionStorage<Traits,ModelType,DirectionType>(traits,parameters,equationTraits)),
            stepTimer(false), printTimer(false)
        {
          parameters.setStorage(storage);

          // force initialization when used
          if (DirectionType::isAdjoint())
          {
            time = -1.;
            endTime = 0.;
          }
          else
          {
            time = 1.;
            endTime = 0.;
          }
        }

          ~DifferentialEquation()
          {
            printTimers();
          }

          /**
           *  @brief Initialize all objects and variables to compute a new transient solution
           */
          void initialize(
              const std::shared_ptr<const ParameterList>& parameterList,
              const std::shared_ptr<const MeasurementList>& measurementList,
              const std::shared_ptr<Measurements>& measurements_
              )
          {
            // set up parameter class
            parameters.setParameterList(parameterList);
            parameters.setMeasurementList(measurementList);
            measurements = measurements_;

            // set up solution vectors
            oldSolution = GridVector(equationTraits.gfs(),0.);
            newSolution = GridVector(equationTraits.gfs(),0.);

            // define time domain
            if (DirectionType::isAdjoint())
            {
              endTime   = traits.config().template get<RF>("time.start");
              startTime = traits.config().template get<RF>("time.end");
            }
            else
            {
              startTime = traits.config().template get<RF>("time.start");
              endTime   = traits.config().template get<RF>("time.end");
            }

            minTimestep  = parameters.minTimestep();
            maxTimestep  = parameters.maxTimestep();
            time         = startTime;
            oldTime      = startTime;
            usedTimestep = 0.;

            // empty storage
            (*storage).clear();

            // interpolate initial conditions and constraints
            InitialValue<Traits,ModelType,DirectionType>
              initial(traits,equationTraits.gfs().gridView(),parameters);
            cg.clear();
            Dune::PDELab::constraints(equationTraits.gfs(),cg);
            Dune::PDELab::interpolate(initial,equationTraits.gfs(),oldSolution);

            newSolution = oldSolution;

            // set up solver object
            solver = std::make_shared<Solver<Traits,ModelType,DirectionType> >
              (traits,equationTraits,cg,parameters,endTime);

            printTimers();
          }

          /**
           * @brief Suggest timestep that would be optimal for equation
           */
          RF suggestTimestep() const
          {
            RF stepSuggested;
            RF timeSuggested;
            if (DirectionType::isAdjoint())
            {
              // ensure we don't walk past the end
              timeSuggested = std::max(endTime, time - maxTimestep);
              // ensure we don't skip measurements
              timeSuggested = (*measurements).suggestPreviousTime(timeSuggested, time);
              stepSuggested = timeSuggested - time;
              // ensure timesteps don't get too small
              stepSuggested = std::min(- minTimestep,stepSuggested);
            }
            else
            {
              // ensure we don't walk past the end
              timeSuggested = std::min(endTime, time + maxTimestep);
              // ensure we don't skip measurements
              timeSuggested = (*measurements).suggestNextTime(time, timeSuggested);
              stepSuggested = timeSuggested - time;
              // ensure timesteps don't get too small
              stepSuggested = std::max(minTimestep,stepSuggested);
            }

            return stepSuggested;
          }

          /**
           * @brief Compute a time step
           */
          void step(RF timestep)
          {
            if (finished())
              DUNE_THROW(Dune::Exception,"equation cannot step when it has finished or wasn't initialized");

            // perform time step
            stepTimer.start();
            (*measurements).setTimes(time,time+timestep);
            unsigned int subSteps = (*solver).step(time,timestep,oldSolution,newSolution);
            stepTimer.stop();

            if (traits.rank() == 0)
              std::cout << "from: " << time << " to: " << time + timestep
                << " stepwidth: " << timestep << " substeps: " << subSteps << std::endl;

            // accept time step
            oldTime = time;
            usedTimestep = timestep;

            time += timestep;

            (*storage).storeSolution(time,newSolution);

            // store solution for time == startTime
            if (oldTime == startTime)
            {
              if (FormulationType::isTransient())
                (*storage).storeSolution(oldTime,oldSolution);
              else
                (*storage).storeSolution(oldTime,newSolution);
            }

            oldSolution = newSolution;

            // print solution if selected
            if (traits.config().template get<bool>("output.writeVTK",false))
            {
              printTimer.start();
              
              std::stringstream ss;
              ss << time;
              std::string timeString(ss.str());
              if (DirectionType::isAdjoint())
              {
                (*storage).printValue   (time,parameters.name()+".adjointValue."
                    +timeString,parameters.name()+".adjointValue");
                (*storage).printFlux    (time,parameters.name()+".adjointFlux."
                    +timeString,parameters.name()+".adjointFlux");
                (*storage).printGradient(time,parameters.name()+".adjointGradient."
                    +timeString,parameters.name()+".adjointGradient");
              }
              else
              {
                (*storage).printValue   (time,parameters.name()+".forwardValue."
                    +timeString,parameters.name()+".forwardValue");
                //(*storage).printFlux    (time,parameters.name()+".forwardFlux."
                //    +timeString,parameters.name()+".forwardFlux");
                //(*storage).printGradient(time,parameters.name()+".forwardGradient."
                //    +timeString,parameters.name()+".forwardGradient");
              }
              
              printTimer.stop();
            }
          }

          /**
           * @brief Return time used in the last step
           */
          RF getTime() const
          {
            return time;
          }

          /**
           * @brief Return timestep used in the last step
           */
          RF getTimestep() const
          {
            return usedTimestep;
          }

          /**
           * @brief Whether the OSM has reached the other side of the time interval
           */
          bool finished() const
          {
            if (DirectionType::isAdjoint())
              return time <= endTime + 1e-6;
            else
              return time >= endTime - 1e-6;
          }

          /**
           * @brief Extract measurements from solution storage
           */
          void extractMeasurements()
          {
            (*measurements).extract(storage,time - usedTimestep,time);
          }

          /**
           * @brief Print timer information
           */
          void printTimers()
          {
            if (traits.rank() == 0)
            {
              if (stepTimer.elapsed()    > 1e-3)
                std::cout << "Time for "
                  << (DirectionType::isAdjoint()?"adjoint ":"")
                  << "equation " << parameters.name()
                  << " solve " << stepTimer.elapsed()  << std::endl;
              
              if (printTimer.elapsed()   > 1e-3)
                std::cout << "Time for "
                  << (DirectionType::isAdjoint()?"adjoint ":"")
                  << "equation " << parameters.name()
                  << " print " << printTimer.elapsed() << std::endl;
            }

            stepTimer.reset();
            printTimer.reset();
          }

      };
  }
}

#endif // DUNE_MODELLING_EQUATION_HH