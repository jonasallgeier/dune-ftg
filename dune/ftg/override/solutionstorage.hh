// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MODELLING_SOLUTIONSTORAGE_HH
#define DUNE_MODELLING_SOLUTIONSTORAGE_HH

#include<list>

#include<dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include<dune/modelling/fluxreconstruction.hh>

namespace Dune {
  namespace Modelling {

    /**
     * @brief Container for solution storage that keeps full history
     */
    template<typename Time, typename Solution>
      class FullContainer
      {
        std::list<std::pair<Time,std::shared_ptr<Solution> > > list;
        Time firstTime, lastTime;

        public:

        FullContainer()
        {
          clear();
        }

        /**
         * @brief Insert solution into container
         */
        void insertSolution(Time time, const Solution& solution)
        {
          std::shared_ptr<Solution> solutionPtr(new Solution(solution));

          if (time > lastTime)
          {
            list.push_back(std::pair<Time,std::shared_ptr<Solution> >(time, solutionPtr));
            lastTime = time;
            if (time < firstTime)
              firstTime = time;
          }
          else if (time < firstTime)
          {
            list.push_front(std::pair<Time,std::shared_ptr<Solution> >(time, solutionPtr));
            firstTime = time;
          }
          else
          {
            typename std::list<std::pair<Time, std::shared_ptr<Solution> > >::iterator it = list.begin();
            while (it != list.end() && it->first < time)
              ++it;
            list.insert(it, std::pair<Time,std::shared_ptr<Solution> >(time, solutionPtr));
          }

          if (time < firstTime)
            firstTime = time;
          if (time > lastTime)
            lastTime = time;
        }

        /**
         * @brief Retrieve the two solutions that are closest to a given time
         */
        void setSolutions(Time time, Time& timeLower, Time& timeUpper,
            std::shared_ptr<Solution>& solutionLower, std::shared_ptr<Solution>& solutionUpper) const
        {
          typename std::list<std::pair<Time,std::shared_ptr<Solution> > >::const_iterator it  = list.begin();
          typename std::list<std::pair<Time,std::shared_ptr<Solution> > >::const_iterator it2 = list.begin();

          if (it == list.end())
            DUNE_THROW(Dune::Exception,"SolutionContainer was empty");

          it2 = it;

          while (it != list.end() && time > it->first)
            it2 = it++;

          if (it == list.end())
            it = it2;

          timeLower = it2->first;
          timeUpper = it->first;

          solutionLower = it2->second;
          solutionUpper = it->second;

          if (std::abs(timeUpper - timeLower) < 1e-6)
          {
            timeLower -= 1e-6;
            timeUpper += 1e-6;
          }

          if (timeLower - time > 2e-6 || time - timeUpper > 2e-6)
            DUNE_THROW(Dune::Exception,"SolutionStorage didn't contain time");
        }

        /**
         * @brief Empty storage container
         */
        void clear()
        {
          list.clear();

          firstTime =  std::numeric_limits<Time>::max();
          lastTime  = -std::numeric_limits<Time>::max();
        }

      };

    /**
     * @brief Container for solution storage that keeps last two steps
     */
    template<typename Time, typename Solution>
      class LastTwoContainer
      {
        std::pair<Time,std::shared_ptr<Solution> > previous;
        std::pair<Time,std::shared_ptr<Solution> > current;

        public:

        LastTwoContainer()
          : previous(0.,std::shared_ptr<Solution>()), current(0.,std::shared_ptr<Solution>())
        {}

        /**
         * @brief Insert solution into container
         */
        void insertSolution(Time time, const Solution& solution)
        {
          std::shared_ptr<Solution> solutionPtr(new Solution(solution));

          previous = current;
          current  = std::pair<Time,std::shared_ptr<Solution> >(time,solutionPtr);
        }

        void setSolutions(Time time, Time& timeLower, Time& timeUpper,
            std::shared_ptr<Solution>& solutionLower, std::shared_ptr<Solution>& solutionUpper) const
        {
          if (current.second)
          {
            if (previous.second)
            {
              if (previous.first < current.first)
              {
                timeLower = previous.first;
                timeUpper = current.first;

                solutionLower = previous.second;
                solutionUpper = current.second;
              }
              else
              {
                timeLower = current.first;
                timeUpper = previous.first;

                solutionLower = current.second;
                solutionUpper = previous.second;
              }
            }
            else
            {
              timeLower = current.first;
              timeUpper = current.first;

              solutionLower = current.second;
              solutionUpper = current.second;
            }

            if (std::abs(timeUpper - timeLower) < 1e-6)
            {
              timeLower -= 1e-6;
              timeUpper += 1e-6;
            }

            if (timeLower - time > 2e-6 || time - timeUpper > 2e-6)
              DUNE_THROW(Dune::Exception,"SolutionStorage didn't contain time");
          }
          else
            DUNE_THROW(Dune::Exception,"SolutionContainer was empty");
        }

        /**
         * @brief Empty storage container
         */
        void clear()
        {
          previous = std::pair<Time,std::shared_ptr<Solution> >(0.,std::shared_ptr<Solution>());
          current  = std::pair<Time,std::shared_ptr<Solution> >(0.,std::shared_ptr<Solution>());
        }

      };

    /**
     * @brief Evaluation policy that picks previous time (explicit methods)
     */
    template<typename Traits>
      class PreviousTimestep
      {
        private:

          using Domain  = typename Traits::GridTraits::Domain;
          using Element = typename Traits::GridTraits::Element;

        public:

          template<typename Time, typename DGF, typename Range>
            static void evaluate(
                Time prevTime, const std::shared_ptr<DGF>& previous,
                Time nextTime, const std::shared_ptr<DGF>& next,
                Time time, const Element& elem, const Domain& x, Range& output
                )
            {
              (*previous).evaluate(elem,x,output);
            }
      };

    /**
     * @brief Evaluation policy that picks next time (implicit methods)
     */
    template<typename Traits>
      class NextTimestep
      {
        private:

          using Domain  = typename Traits::GridTraits::Domain;
          using Element = typename Traits::GridTraits::Element;

        public:

          template<typename Time, typename DGF, typename Range>
            static void evaluate(
                Time prevTime, const std::shared_ptr<DGF>& previous,
                Time nextTime, const std::shared_ptr<DGF>& next,
                Time time, const Element& elem, const Domain& x, Range& output
                )
            {
              (*next).evaluate(elem,x,output);
            }
      };

    /**
     * @brief Evaluation policy that interpolates linearly
     */
    template<typename Traits>
      class LinearInterpolation
      {
        private:

          using Domain  = typename Traits::GridTraits::Domain;
          using Element = typename Traits::GridTraits::Element;

        public:

          template<typename Time, typename DGF, typename Range>
            static void evaluate(
                Time prevTime, const std::shared_ptr<DGF>& previous,
                Time nextTime, const std::shared_ptr<DGF>& next,
                Time time, const Element& elem, const Domain& x, Range& output
                )
            {
              Range prevVal, nextVal;
              (*previous).evaluate(elem,x,prevVal);
              (*next    ).evaluate(elem,x,nextVal);

              for (unsigned int i = 0; i < output.size(); i++)
                output[i] = (nextVal[i] - prevVal[i]) / (nextTime - prevTime);
            }
      };

    /**
     * @brief Storage facility for transient solutions of PDEs
     */
    template<typename Traits, typename ModelType, typename DirectionType>
      class SolutionStorage
      {
        using GridTraits = typename Traits::GridTraits;

        enum {dim = GridTraits::dim};

        using GV           = typename GridTraits::GridView;
        using RF           = typename GridTraits::RangeField;
        using Scalar       = typename GridTraits::Scalar;
        using Vector       = typename GridTraits::Vector;
        using DF           = typename GridTraits::DomainField;
        using Domain       = typename GridTraits::Domain;
        using IDomain      = typename GridTraits::IDomain;
        using Element      = typename GridTraits::Element;
        using Intersection = typename GridTraits::Intersection;

        using GFS = typename EquationTraits<Traits,ModelType,DirectionType>::GridFunctionSpace;
        using GridVector = typename EquationTraits<Traits,ModelType,DirectionType>::GridVector;

        using ScalarDGF    = PDELab::DiscreteGridFunction<GFS,GridVector>;
        using LocalGradDGF = PDELab::DiscreteGridFunctionGradient<GFS,GridVector>;

        // extract flux reconstruction from equation traits
        template<typename... T>
          using FluxReconstruction = typename EquationTraits<Traits,ModelType,DirectionType>
          ::template FluxReconstruction<T...>;

        using GradientDGF = FluxReconstruction<Traits, ModelType, DirectionType, Reconstruction::Gradient>;
        using FluxDGF     = FluxReconstruction<Traits, ModelType, DirectionType, Reconstruction::Flux>;
        using LinFluxDGF  = FluxReconstruction<Traits, ModelType, DirectionType, Reconstruction::LinFlux>;
        
        const Traits&                                         traits;
        const ModelParameters<Traits,ModelType>&              parameters;
        const EquationTraits<Traits,ModelType,DirectionType>& equationTraits;
        const Boundary<Traits,ModelType,DirectionType>        boundary;

        // extract storage container policy from equation traits
        template<typename... T>
          using StorageContainer = typename EquationTraits<Traits,ModelType,DirectionType>
          ::template StorageContainer<T...>;

        StorageContainer<RF,GridVector> solutionContainer;

        // extract interpolation policy from equation traits
        template<typename... T>
          using TemporalInterpolation = typename EquationTraits<Traits,ModelType,DirectionType>
          ::template TemporalInterpolation<T...>;

        mutable RF timeLower, timeUpper, scalarTimeLower, scalarTimeUpper,
                gradientTimeLower, gradientTimeUpper, fluxTimeLower, fluxTimeUpper,
                linFluxTimeLower, linFluxTimeUpper;

        mutable std::shared_ptr<GridVector>   solutionLower, solutionUpper;
        mutable std::shared_ptr<ScalarDGF>    lowerScalarDGF, upperScalarDGF;
        mutable std::shared_ptr<LocalGradDGF> lowerLocalGradDGF, upperLocalGradDGF;
        mutable std::shared_ptr<GradientDGF>  lowerGradientDGF, upperGradientDGF;
        mutable std::shared_ptr<FluxDGF>      lowerFluxDGF, upperFluxDGF;
        mutable std::shared_ptr<LinFluxDGF>   lowerLinFluxDGF, upperLinFluxDGF;

        public:

        SolutionStorage(
            const Traits& traits_,
            const ModelParameters<Traits,ModelType>& parameters_,
            const EquationTraits<Traits,ModelType,DirectionType>& equationTraits_
            )
          : traits(traits_), parameters(parameters_), equationTraits(equationTraits_),
          boundary(traits,parameters.name())
        {
          clear();
        }

        /**
         * @brief Empty storage to facilitate next computation
         */
        void clear()
        {
          solutionContainer.clear();

          timeLower         =  std::numeric_limits<RF>::max();
          timeUpper         = -std::numeric_limits<RF>::max();
          scalarTimeLower   =  std::numeric_limits<RF>::max();
          scalarTimeUpper   = -std::numeric_limits<RF>::max();
          gradientTimeLower =  std::numeric_limits<RF>::max();
          gradientTimeUpper = -std::numeric_limits<RF>::max();
          fluxTimeLower     =  std::numeric_limits<RF>::max();
          fluxTimeUpper     = -std::numeric_limits<RF>::max();
          linFluxTimeLower  =  std::numeric_limits<RF>::max();
          linFluxTimeUpper  = -std::numeric_limits<RF>::max();
        }

        /**
         * @brief Store a solution for later use
         */
        void storeSolution(const RF time, const GridVector& solution)
        {
          solutionContainer.insertSolution(time,solution);
        }

        /**
         * @brief Evaluate solution at given coordinate and time
         */
        void value(const RF time, const Element& elem, const Domain& x, Scalar& output) const
        {
          updateScalarDGFs(time);

          if (DirectionType::isAdjoint())
            TemporalInterpolation<Traits>::evaluate(timeUpper,upperScalarDGF,
                timeLower,lowerScalarDGF,time,elem,x,output);
          else
            TemporalInterpolation<Traits>::evaluate(timeLower,lowerScalarDGF,
                timeUpper,upperScalarDGF,time,elem,x,output);
        }

        /**
         * @brief Evaluate temporal derivative of solution at given coordinate and time
         */
        void tempDeriv(const RF time, const Element& elem, const Domain& x, Scalar& output) const
        {
          updateScalarDGFs(time);

          if (std::abs(timeLower - timeUpper) < 1e-5)
            output = 0.;
          else
          {
            Scalar lower,upper;
            evaluateStorage(lowerScalarDGF,elem,x,lower);
            evaluateStorage(upperScalarDGF,elem,x,upper);

            for (unsigned int i = 0; i < output.size(); i++)
              output[i] = (upper[i] - lower[i]) / (timeUpper - timeLower);
          }
        }

        /**
         * @brief Evaluate gradient of solution at given coordinate and time
         */
        void gradient(const RF time, const Element& elem, const Domain& x, Vector& output) const
        {
          updateGradientDGFs(time);

          if (DirectionType::isAdjoint())
            TemporalInterpolation<Traits>::evaluate(timeUpper,upperGradientDGF,
                timeLower,lowerGradientDGF,time,elem,x,output);
          else
            TemporalInterpolation<Traits>::evaluate(timeLower,lowerGradientDGF,
                timeUpper,upperGradientDGF,time,elem,x,output);
        }

        /**
         * @brief Evaluate flux of solution at given coordinate and time
         */
        void flux(const RF time, const Element& elem, const Domain& x, Vector& output) const
        {
          updateFluxDGFs(time);

          if (DirectionType::isAdjoint())
            TemporalInterpolation<Traits>::evaluate(timeUpper,upperFluxDGF,
                timeLower,lowerFluxDGF,time,elem,x,output);
          else
            TemporalInterpolation<Traits>::evaluate(timeLower,lowerFluxDGF,
                timeUpper,upperFluxDGF,time,elem,x,output);
        }

        /**
         * @brief Evaluate derivative of flux of solution w.r.t. solution
         */
        void linFlux(const RF time, const Element& elem, const Domain& x, Vector& output) const
        {
          updateLinFluxDGFs(time);

          if (DirectionType::isAdjoint())
            TemporalInterpolation<Traits>::evaluate(timeUpper,upperLinFluxDGF,
                timeLower,lowerLinFluxDGF,time,elem,x,output);
          else
            TemporalInterpolation<Traits>::evaluate(timeLower,lowerLinFluxDGF,
                timeUpper,upperLinFluxDGF,time,elem,x,output);
        }

        /**
         * @brief Create VTK output of solution
         */
        void printValue(const RF time, const std::string& fileName, const std::string& dataName) const
        {
          updateScalarDGFs(time);

          /* this does not work... maybe has to be managed with some kind of pointers...
          std::shared_ptr<Dune::SubsamplingVTKWriter<GV> > ssvtk = std::make_shared<Dune::SubsamplingVTKWriter<GV> > (equationTraits.gfs().gridView(),1);
          Dune::VTKSequenceWriter<GV> vtk_writer(ssvtk,fileName,"vtk","");
 
          std::shared_ptr<Dune::PDELab::VTKGridFunctionAdapter<ScalarDGF> >
            dgfPtr(new Dune::PDELab::VTKGridFunctionAdapter<ScalarDGF>(lowerScalarDGF,dataName));
          vtk_writer.addVertexData(dgfPtr);
          vtk_writer.write(time);
          */
          Dune::SubsamplingVTKWriter<GV> vtkwriter(equationTraits.gfs().gridView(),1);
          std::shared_ptr<Dune::PDELab::VTKGridFunctionAdapter<ScalarDGF> >
            dgfPtr(new Dune::PDELab::VTKGridFunctionAdapter<ScalarDGF>(lowerScalarDGF,dataName));
          vtkwriter.addVertexData(dgfPtr);
          vtkwriter.pwrite(fileName,"vtk","",Dune::VTK::appendedraw);
        }

        /**
         * @brief Create VTK output of gradient of solution
         */
        void printGradient(const RF time, const std::string& fileName, const std::string& dataName) const
        {
          updateGradientDGFs(time);

          Dune::SubsamplingVTKWriter<GV> vtkwriter(equationTraits.gfs().gridView(),2);
          std::shared_ptr<Dune::PDELab::VTKGridFunctionAdapter<GradientDGF> >
            dgfPtr(new Dune::PDELab::VTKGridFunctionAdapter<GradientDGF>(lowerGradientDGF,dataName));
          vtkwriter.addVertexData(dgfPtr);
          vtkwriter.pwrite(fileName,"vtk","",Dune::VTK::appendedraw);
        }

        /**
         * @brief Create VTK output of flux of solution
         */
        void printFlux(const RF time, const std::string& fileName, const std::string& dataName) const
        {
          updateFluxDGFs(time);

          Dune::SubsamplingVTKWriter<GV> vtkwriter(equationTraits.gfs().gridView(),2);
          std::shared_ptr<Dune::PDELab::VTKGridFunctionAdapter<FluxDGF> >
            dgfPtr(new Dune::PDELab::VTKGridFunctionAdapter<FluxDGF>(lowerFluxDGF,dataName));
          vtkwriter.addVertexData(dgfPtr);
          vtkwriter.pwrite(fileName,"vtk","",Dune::VTK::appendedraw);
        }

        /**
         * @brief Create VTK output of derivative of flux of solution
         */
        void printLinFlux(const RF time, const std::string& fileName, const std::string& dataName) const
        {
          updateLinFluxDGFs(time);

          Dune::SubsamplingVTKWriter<GV> vtkwriter(equationTraits.gfs().gridView(),2);
          std::shared_ptr<Dune::PDELab::VTKGridFunctionAdapter<LinFluxDGF> >
            dgfPtr(new Dune::PDELab::VTKGridFunctionAdapter<LinFluxDGF>(lowerLinFluxDGF,dataName));
          vtkwriter.addVertexData(dgfPtr);
          vtkwriter.pwrite(fileName,"vtk","",Dune::VTK::appendedraw);
        }

        private:

        /**
         * @brief Update pair of solutions at interval boundary
         */
        void setSolutions(const RF time) const
        {
          solutionContainer.setSolutions(time, timeLower, timeUpper, solutionLower, solutionUpper);
        }

        /**
         * @brief Update pair of scalar DGFs
         */
        void updateScalarDGFs(const RF time) const
        {
          if (time < scalarTimeLower || time > scalarTimeUpper)
          {
            setSolutions(time);
            scalarTimeLower = timeLower;
            scalarTimeUpper = timeUpper;

            lowerScalarDGF = std::make_shared<ScalarDGF>(equationTraits.gfs(),*solutionLower);
            upperScalarDGF = std::make_shared<ScalarDGF>(equationTraits.gfs(),*solutionUpper);

            lowerLocalGradDGF = std::make_shared<LocalGradDGF>(equationTraits.gfs(),*solutionLower);
            upperLocalGradDGF = std::make_shared<LocalGradDGF>(equationTraits.gfs(),*solutionUpper);
          }
        }

        /**
         * @brief Update pair of gradient vector DGFs
         */
        void updateGradientDGFs(const RF time) const
        {
          if (time < scalarTimeLower || time > scalarTimeUpper
              || time < gradientTimeLower || time > gradientTimeUpper)
          {
            updateScalarDGFs(time);
            gradientTimeLower = timeLower;
            gradientTimeUpper = timeUpper;

            lowerGradientDGF = std::make_shared<GradientDGF>(traits,equationTraits,parameters,
                boundary,lowerScalarDGF,lowerLocalGradDGF,timeLower);
            upperGradientDGF = std::make_shared<GradientDGF>(traits,equationTraits,parameters,
                boundary,upperScalarDGF,upperLocalGradDGF,timeUpper);
          }
        }

        /**
         * @brief Update pair of flux vector DGFs
         */
        void updateFluxDGFs(const RF time) const
        {
          if (time < scalarTimeLower || time > scalarTimeUpper
              || time < fluxTimeLower || time > fluxTimeUpper)
          {
            updateScalarDGFs(time);
            fluxTimeLower = timeLower;
            fluxTimeUpper = timeUpper;

            lowerFluxDGF = std::make_shared<FluxDGF>(traits,equationTraits,parameters,
                boundary,lowerScalarDGF,lowerLocalGradDGF,timeLower);
            upperFluxDGF = std::make_shared<FluxDGF>(traits,equationTraits,parameters,
                boundary,upperScalarDGF,upperLocalGradDGF,timeUpper);
          }
        }

        /**
         * @brief Update pair of flux derivative vector DGFs
         */
        void updateLinFluxDGFs(const RF time) const
        {
          if (time < scalarTimeLower || time > scalarTimeUpper
              || time < linFluxTimeLower || time > linFluxTimeUpper)
          {
            updateScalarDGFs(time);
            linFluxTimeLower = timeLower;
            linFluxTimeUpper = timeUpper;

            lowerLinFluxDGF = std::make_shared<LinFluxDGF>(traits,equationTraits,parameters,
                boundary,lowerScalarDGF,lowerLocalGradDGF,timeLower);
            upperLinFluxDGF = std::make_shared<LinFluxDGF>(traits,equationTraits,parameters,
                boundary,upperScalarDGF,upperLocalGradDGF,timeUpper);
          }
        }

      };

  }
}

#endif // DUNE_MODELLING_SOLUTIONSTORAGE_HH
