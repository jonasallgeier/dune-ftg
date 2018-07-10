// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MODELLING_FORWARD_HH
#define DUNE_MODELLING_FORWARD_HH

#include<list>

#include<dune/modelling/declarations.hh>
#include<dune/modelling/equation.hh>

namespace Dune {
  namespace Modelling {

    /**
     * @brief Helper that allows reversed direction in range-based for loops
     */
    template<typename Iterable>
      class ReversionWrapper
      {
        Iterable& iterable;

        public:

        ReversionWrapper(Iterable& iterable_)
          : iterable(iterable_)
        {}

        auto begin()
        {
          return iterable.rbegin();
        }

        auto end()
        {
          return iterable.rend();
        }
      };

    /**
     * @brief Free function that reverses direction in range-based for loop
     */
    template<typename Iterable>
      auto reverse(Iterable& iterable)
      {
        return ReversionWrapper<Iterable>(iterable);
      }

    /**
     * @brief Abstract class for models
     */
    template<typename Traits>
      class ForwardModelBase
      {

        public:

          virtual ~ForwardModelBase() {};

          virtual void initForward(
              const std::shared_ptr<const typename Traits::ParameterList>& parameterList,
              //const std::shared_ptr<typename Traits::MeasurementList>& measurementList,
              //const std::shared_ptr<typename Traits::MeasurementList::SubMeasurements>& measurements
              const std::shared_ptr<typename Traits::MeasurementList>& measurements,
              const std::shared_ptr<typename Traits::ERTMatrixContainer>& ertMatrixContainer
              ) = 0;

          virtual void solveForward() = 0;
          virtual void stepForward(typename Traits::GridTraits::RangeField timestep) = 0;
          virtual void clearStorage() = 0;
          virtual typename Traits::GridTraits::RangeField suggestTimestepForward() const = 0;
          virtual bool forwardFinished() const = 0;
          virtual typename Traits::GridTraits::RangeField getTime() const = 0;
      };

    /**
     * @brief Concrete model based on PDE
     */
    template<typename Traits, typename ModelType, typename FormulationType,
      template <typename> class ModelBase = ForwardModelBase>
      class ForwardModel
      : public ModelBase<Traits>
      {
        protected:

          const std::string                                             name;
          std::shared_ptr<ModelParameters<Traits,ModelType> >           parameters;
          Equation<Traits,ModelType,FormulationType,Direction::Forward> equation;

        public:

          /**
           * @brief Constructor
           */
          ForwardModel(
              const Traits& traits,
              const std::string& name_
              )
            : name(name_), parameters(new ModelParameters<Traits,ModelType>(traits,name)),
            equation(traits,*parameters)
        {}

          /**
           * @brief Models are fully defined by their name and should not be copied
           */
          ForwardModel(const ForwardModel& other) = delete;

          /**
           * @brief Access name of model
           */
          const std::string& getName() const
          {
            return name;
          }

          /**
           * @brief Return parameters of the model
           */
          const std::shared_ptr<ModelParameters<Traits,ModelType> >& getParameters() const
          {
            return parameters;
          }

          /**
           * @brief Solve the equation
           */
          virtual void solveForward()
          {
            while(!forwardFinished())
            {
              using RF = typename Traits::GridTraits::RangeField;
              const RF timestep = equation.suggestTimestep();
              equation.step(timestep);
              equation.extractMeasurements();
            }
          }

          /**
           * @brief Initialize the equation to compute again
           */
          virtual void initForward(
              const std::shared_ptr<const typename Traits::ParameterList>& parameterList,
              //const std::shared_ptr<typename Traits::MeasurementList>& measurementList,
              //const std::shared_ptr<typename Traits::MeasurementList::SubMeasurements>& measurements
              const std::shared_ptr<typename Traits::MeasurementList>& measurements,
              const std::shared_ptr<typename Traits::ERTMatrixContainer>& ertMatrixContainer
              )
          {
            equation.initialize(parameterList,measurements,ertMatrixContainer);//,measurementList,measurements);
          }

          /**
           * @brief Compute next step of the equation
           */
          virtual void stepForward(typename Traits::GridTraits::RangeField timestep)
          {
            equation.step(timestep);
            equation.extractMeasurements();
          }

          /**
           * @brief Suggest optimal time stepsize for next step
           */
          virtual typename Traits::GridTraits::RangeField suggestTimestepForward() const
          {
            return equation.suggestTimestep();
          }

          virtual void clearStorage()
          {
            return equation.clearStorage();
          }

          virtual typename Traits::GridTraits::RangeField getTime() const
          {
            return equation.getTime();
          }

          /**
           * @brief Check if equation is solved
           */
          virtual bool forwardFinished() const
          {
            return equation.finished();
          }
      };

    /**
     * @brief Collection of forward models with given names
     */
    template<typename Traits, template<typename> class ModelBase = ForwardModelBase>
      class ForwardModelList
      {
        protected:

          using ParameterList   = typename Traits::ParameterList;
          using MeasurementList = typename Traits::MeasurementList;
          using ERTMatrixContainer = typename Traits::ERTMatrixContainer;

          using ModelPair   = std::pair<std::string, std::shared_ptr<ModelBase<Traits> > >;
          using ModelVector = std::vector<std::pair<std::string, std::shared_ptr<ModelBase<Traits> > > >;

          const Traits& traits;
          unsigned int highest_moment;
          ModelVector list;

          std::vector<bool> isInitialized;
        public:

          /**
           * @brief Constructor
           */
          ForwardModelList(const Traits& traits_,unsigned int highest_moment_)
            : traits(traits_), highest_moment(highest_moment_)
          {}

          /**
           * @brief Add a model to the list
           */
          template<typename Model, typename ...Args>
            void add(
                const std::string& name,
                const std::list<std::string>& otherNames = std::list<std::string>()
                )
            {
              // check for name clashes
              for (const ModelPair& pair : list)
                if (pair.first == name)
                  DUNE_THROW(Dune::Exception,"Model of name " + name + " was already in the list");

              // create model and insert it in list
              const std::shared_ptr<Model> model(new Model(traits,name));
              list.push_back({name,model});
              isInitialized.push_back(false);

              // register parameters of given models in each other
              std::list<std::string> otherNamesCopy = otherNames;
              registerModel<Model,Args...>(name,model,otherNamesCopy);
            }

          /**
           * @brief Iterate until all models have been computed
           */
          //virtual void solve(
          void solve(
              const std::shared_ptr<const ParameterList>& parameterList,
              const std::shared_ptr<MeasurementList>& measurementList,
              const std::shared_ptr<ERTMatrixContainer>& ertMatrixContainer
              )
          {
            solveForward(parameterList,measurementList,ertMatrixContainer);
          }

          /**
           * @brief Iterate forward until all equations are computed, in bulk
           */
          void solveForward(
              const std::shared_ptr<const ParameterList>& parameterList,
              const std::shared_ptr<MeasurementList>& measurementList,
              const std::shared_ptr<ERTMatrixContainer>& ertMatrixContainer
              )
          {
            //initForward(parameterList,measurementList,ertMatrixContainer);

            //while (!forwardFinished())
            //{
              using RF = typename Traits::GridTraits::RangeField;

              for (unsigned int i = 0; i!=list.size();i++)
              {
                if (isInitialized[i] == false)
                {
                  list[i].second->initForward(parameterList,measurementList,ertMatrixContainer);
                  isInitialized[i] = true;
                }
                RF timestep = list[i].second->suggestTimestepForward();
                list[i].second->stepForward(timestep);
                
                // if the last moment model for an electrode was evaluated -> empty the storages!
                if (list[i].first.find("momentsERT_") == 0)
                {
                  std::string str = list[i].first;
                  std::string common_base = "momentsERT_";
                  auto start_position_to_erase = str.find(common_base);
                  str.erase(start_position_to_erase, common_base.size());      
                  std::replace(str.begin(), str.end(), '_', ' ');

                  std::stringstream ss; 
                  ss << str;
                  std::string temp;

                  ss >> temp; // get the ERT model number
                  unsigned int model_number;
                  std::stringstream(temp) >> model_number;

                  if (model_number > 0) // delete all models after usage; except for the ones of electrode 1
                  {
                    ss >> temp; // get the moment number k
                    unsigned int k;
                    std::stringstream(temp) >> k; 

                    if (k==highest_moment)
                    {
                      for (unsigned int j = 0; j<=highest_moment;j++)
                      {
                        list[i-j].second->clearStorage();
                        (list[i-j].second).reset();                    
                      }
                      list[i-highest_moment-1].second->clearStorage();
                      (list[i-highest_moment-1].second).reset();
                    }
                  }
                }
              //}
            }
          }

          /**
           * @brief Print the models that are members of this list
           */
          void report(std::ostream& out) const
          {
            out << "Models in list:" << std::endl;
            for(const ModelPair& pair : list)
              out << pair.first << "; ";
            out << std::endl;
          }

        private:

          /**
           * @brief Allow mutual registration of parameters of different models
           */
          template<typename Model, typename OtherModel, typename ...Args>
            void registerModel(
                const std::string& name,
                const std::shared_ptr<Model>& model,
                std::list<std::string> otherNames
                )
            {
              if (otherNames.empty())
                DUNE_THROW(Dune::Exception,"Expected additional names for models to register in list");

              // try to find model with given name
              const std::string otherName(otherNames.front());
              const std::function<bool(const ModelPair&)> hasOtherName
                = [&](const ModelPair& elem){return elem.first == otherName;};
              const auto iter = std::find_if(list.begin(),list.end(),hasOtherName);
              if (iter == list.end())
                DUNE_THROW(Dune::Exception,"Model of name " + otherName + " not found in the list");

              // try to extract model with correct type from given name
              std::shared_ptr<OtherModel> otherModel = std::dynamic_pointer_cast<OtherModel>(iter->second);
              if (!otherModel)
                DUNE_THROW(Dune::Exception,"Model of name " + otherName + " didn't have expected type");

              // register pointers to model parameters in each other
              otherModel->getParameters()->template registerModel(name,     model     ->getParameters());
              model     ->getParameters()->template registerModel(otherName,otherModel->getParameters());
              otherNames.pop_front();

              // handle additional model clients (if any)
              registerModel<Model,Args...>(name,model,otherNames);
            }

          /**
           * @brief Specialization for termination of recursion
           */
          template<typename Model>
            void registerModel(
                const std::string& name,
                const std::shared_ptr<Model>& model,
                std::list<std::string> otherNames
                )
            {
              if (!otherNames.empty())
                DUNE_THROW(Dune::Exception,"Found additional names for models to register in list");
            }

        protected:

          /**
           * @brief Reset all equations to be able to compute again
           */
          void initForward(
              const std::shared_ptr<const ParameterList>& parameterList,
              const std::shared_ptr<MeasurementList>& measurementList,
              const std::shared_ptr<ERTMatrixContainer>& ertMatrixContainer
              )
          {
            for (ModelPair& pair : list)
              pair.second->initForward(parameterList,measurementList,ertMatrixContainer);//,
                  //measurementList->get(pair.first));
          }

          /**
           * @brief Check if all equations have been computed
           */
          bool forwardFinished() const
          {
            bool finished = true;
            for(const ModelPair& pair : list)
              finished = finished && pair.second->forwardFinished();

            return finished;
          }

      };
  }
}

#endif // DUNE_MODELLING_FORWARD_HH
