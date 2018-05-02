// -*- tab-width: 2; indent-tabs-mode: nil -*-
#ifndef DUNE_RANDOMFIELD_RANDOMFIELD_HH
#define	DUNE_RANDOMFIELD_RANDOMFIELD_HH

#include<dune/common/parametertree.hh>

#if HAVE_DUNE_FUNCTIONS
// for VTK output functionality
#include<dune/grid/yaspgrid.hh>
#include<dune/grid/io/file/vtk.hh>
#include<dune/functions/gridfunctions/analyticgridviewfunction.hh>
#endif // HAVE_DUNE_FUNCTIONS

#include<dune/randomfield/io.hh>
#include<dune/ftg/override/fieldtraits.hh> //use modified version!
#include<dune/randomfield/trend.hh>
#include<dune/randomfield/stochastic.hh>
#include<dune/randomfield/matrix.hh>
#include<dune/randomfield/mutators.hh>
#include<dune/randomfield/legacyvtk.hh>

namespace Dune {
  namespace RandomField {

    /**
     * @brief Gaussian random field in 1D, 2D or 3D
     */
    template<typename GridTraits, bool storeInvMat = true, bool storeInvRoot = false>
      class RandomField
      {

        public:

          typedef RandomFieldTraits<GridTraits,storeInvMat,storeInvRoot> Traits;
          typedef typename Traits::RF                                    RF;

        protected:

          // to allow reading in constructor
          class ParamTreeHelper
          {
            Dune::ParameterTree config;

            public:

            ParamTreeHelper(const std::string& fileName = "")
            {
              if (fileName != "")
              {
                Dune::ParameterTreeParser parser;
                parser.readINITree(fileName+".field",config);
              }
            }

            const Dune::ParameterTree& get() const
            {
              return config;
            }
          };

          const ParamTreeHelper                            treeHelper;
          const Dune::ParameterTree                        config;
          const ValueTransform<RF>                         valueTransform;
          std::shared_ptr<Traits>                          traits;
          std::shared_ptr<RandomFieldMatrix<Traits> >      matrix;
          TrendPart<Traits>                                trendPart;
          StochasticPart<Traits>                           stochasticPart;
          mutable std::shared_ptr<StochasticPart<Traits> > invMatPart;
          mutable bool                                     invMatValid;
          mutable std::shared_ptr<StochasticPart<Traits> > invRootPart;
          mutable bool                                     invRootValid;

        public:

          /**
           * @brief Constructor reading from file or creating homogeneous field
           */
          template<typename LoadBalance = DefaultLoadBalance<GridTraits::dim> >
            explicit RandomField(const Dune::ParameterTree& config_, const std::string& fileName = "", const LoadBalance& loadBalance = LoadBalance(), const MPI_Comm comm = MPI_COMM_WORLD)
            : config(config_), valueTransform(config), traits(new Traits(config,loadBalance,comm)), matrix(new RandomFieldMatrix<Traits>(traits)),
            trendPart(config,traits,fileName), stochasticPart(traits,fileName),
            invMatValid(false), invRootValid(false)
            {
              if (storeInvMat)
                invMatPart = std::shared_ptr<StochasticPart<Traits> >(new StochasticPart<Traits>(stochasticPart));

              if (storeInvRoot)
                invRootPart = std::shared_ptr<StochasticPart<Traits> >(new StochasticPart<Traits>(stochasticPart));
            }

          /**
           * @brief Constructor reading field and config from file
           */
          template<typename LoadBalance = DefaultLoadBalance<GridTraits::dim> >
            RandomField(const std::string& fileName, const LoadBalance& loadBalance = LoadBalance(), const MPI_Comm comm = MPI_COMM_WORLD)
            : treeHelper(fileName), config(treeHelper.get()), valueTransform(config), traits(new Traits(config,loadBalance,comm)), matrix(new RandomFieldMatrix<Traits>(traits)),
            trendPart(config,traits,fileName), stochasticPart(traits,fileName),
            invMatValid(false), invRootValid(false)
            {
              if (storeInvMat)
                invMatPart = std::shared_ptr<StochasticPart<Traits> >(new StochasticPart<Traits>(stochasticPart));

              if (storeInvRoot)
                invRootPart = std::shared_ptr<StochasticPart<Traits> >(new StochasticPart<Traits>(stochasticPart));
            }

          /**
           * @brief Constructor copying traits and covariance matrix
           */
          RandomField(const RandomField& other, const std::string& fileName)
            : config(other.config), valueTransform(other.valueTransform), traits(other.traits), matrix(other.matrix),
            trendPart(config,traits,fileName), stochasticPart(traits,fileName),
            invMatValid(false), invRootValid(false)
        {
          if (storeInvMat)
            invMatPart = std::shared_ptr<StochasticPart<Traits> >(new StochasticPart<Traits>(stochasticPart));

          if (storeInvRoot)
            invRootPart = std::shared_ptr<StochasticPart<Traits> >(new StochasticPart<Traits>(stochasticPart));
        }

#if HAVE_DUNE_PDELAB
           /**
           * @brief Constructor converting from GridFunctionSpace and GridVector
           */
          template<typename GFS, typename Field>
            RandomField(const RandomField& other, const GFS& gfs, const Field& field)
            : config(other.config), valueTransform(other.valueTransform), traits(other.traits), matrix(other.matrix),
            trendPart(other.trendPart,gfs,field), stochasticPart(other.stochasticPart,gfs,field),
            invMatValid(false), invRootValid(false)
        {
          if (storeInvMat)
            invMatPart = std::shared_ptr<StochasticPart<Traits> >(new StochasticPart<Traits>(stochasticPart));

          if (storeInvRoot)
            invRootPart = std::shared_ptr<StochasticPart<Traits> >(new StochasticPart<Traits>(stochasticPart));
        }
#endif // HAVE_DUNE_PDELAB

          /**
           * @brief Copy constructor
           */
          RandomField(const RandomField& other)
            : config(other.config), valueTransform(other.valueTransform), traits(other.traits), matrix(other.matrix),
            trendPart(other.trendPart), stochasticPart(other.stochasticPart),
            invMatValid(other.invMatValid), invRootValid(other.invRootValid)
        {
          if (storeInvMat)
            invMatPart = std::shared_ptr<StochasticPart<Traits> >(new StochasticPart<Traits>(*(other.invMatPart)));

          if (storeInvRoot)
            invRootPart = std::shared_ptr<StochasticPart<Traits> >(new StochasticPart<Traits>(*(other.invRootPart)));
        }

          /**
           * @brief Assignment operator
           */
          RandomField& operator=(const RandomField& other)
          {
            config         = other.config;
            valueTransform = other.valueTransform;
            traits         = other.traits;
            matrix         = other.matrix;
            trendPart      = other.trendPart;
            stochasticPart = other.stochasticPart;
            invMatValid    = other.invMatValid;
            invRootValid   = other.invRootValid;

            if (storeInvMat && invMatValid)
              invMatPart = std::shared_ptr<StochasticPart<Traits> >(new StochasticPart<Traits>(*(other.invMatPart)));

            if (storeInvRoot && invRootValid)
              invRootPart = std::shared_ptr<StochasticPart<Traits> >(new StochasticPart<Traits>(*(other.invRootPart)));

            return *this;
          }

          /**
           * @brief Cell volume of the random field discretization
           */
          RF cellVolume() const
          {
            return (*traits).cellVolume;
          }

          /**
           * @brief Generate a field with the desired correlation structure
           */
          void generate(bool allowNonWorldComm = false)
          {
            // create seed out of current time
            unsigned int seed = (unsigned int) clock();
            // different seeds for different fields
            seed += static_cast<int>(reinterpret_cast<uintptr_t>(&stochasticPart));

            generate(seed,allowNonWorldComm);
          }

          /**
           * @brief Generate a field with desired correlation structure using seed
           */
          void generate(unsigned int seed, bool allowNonWorldComm = false)
          {
            if (((*traits).comm != MPI_COMM_WORLD) && !allowNonWorldComm)
              DUNE_THROW(Dune::Exception,
                  "generation of inconsistent fields prevented, set allowNonWorldComm = true if you really want this");

            std::cout << "generate with seed: " << seed << std::endl;

            (*matrix).generateField(seed,stochasticPart);
            trendPart.generate(seed);
          }

          /**
           * @brief Generate a field without correlation structure (i.e. noise)
           */
          void generateUncorrelated(bool allowNonWorldComm = false)
          {
            // create seed out of current time
            unsigned int seed = (unsigned int) clock();
            // different seeds for different fields
            seed += static_cast<int>(reinterpret_cast<uintptr_t>(&stochasticPart));

            generate(seed,allowNonWorldComm);
          }

          /**
           * @brief Generate a field containing noise using seed
           */
          void generateUncorrelated(unsigned int seed, bool allowNonWorldComm = false)
          {
            if (((*traits).comm != MPI_COMM_WORLD) && !allowNonWorldComm)
              DUNE_THROW(Dune::Exception,
                  "generation of inconsistent fields prevented, set allowNonWorldComm = true if you really want this");

            (*matrix).generateUncorrelatedField(stochasticPart);
            trendPart.generateUncorrelated();
          }

#if HAVE_DUNE_GRID
          /**
           * @brief Evaluate the random field in the coordinates of an element
           */
          template<typename Element>
            void evaluate(
                const Element& elem,
                const typename Traits::DomainType& xElem,
                typename Traits::RangeType& output
                ) const
            {
              const typename Traits::DomainType location = elem.geometry().global(xElem);
              evaluate(location,output);
            }
#endif // HAVE_DUNE_GRID

          /**
           * @brief Evaluate the random field at given coordinates
           */
          void evaluate(const typename Traits::DomainType& location, typename Traits::RangeType& output) const
          {
            typename Traits::RangeType stochastic = 0., trend = 0.;

            stochasticPart.evaluate(location,stochastic);
            trendPart     .evaluate(location,trend);

            output = stochastic + trend;
            valueTransform.apply(output);
          }

          /**
           * @brief Export random field to files on disk
           */
          void writeToFile(const std::string& fileName) const
          {
            stochasticPart.writeToFile(fileName);
            trendPart     .writeToFile(fileName);

            std::ofstream file(fileName+".field",std::ofstream::trunc);
            config.report(file);
          }

          /**
           * @brief Export random field as flat unstructured VTK file, requires dune-grid and dune-functions
           */
          template<typename GridView>
            void writeToVTK(const std::string& fileName, const GridView& gv) const
            {
#if HAVE_DUNE_FUNCTIONS
              Dune::VTKWriter<GridView> vtkWriter(gv,Dune::VTK::conforming);
              auto f = Dune::Functions::makeAnalyticGridViewFunction([&](auto x)
                  {typename Traits::RangeType output; this->evaluate(x,output); return output;},gv);
              vtkWriter.addCellData(f,VTK::FieldInfo(fileName,VTK::FieldInfo::Type::scalar,1));
              vtkWriter.pwrite(fileName,"","",Dune::VTK::appendedraw);
#else //HAVE_DUNE_FUNCTIONS
              DUNE_THROW(Dune::NotImplemented,"Unstructured VTK output requires dune-grid and dune-functions");
#endif //HAVE_DUNE_FUNCTIONS
            }

          /**
           * @brief Export random field as unstructured VTK file, requires dune-grid and dune-functions
           */
          template<typename GridView>
            void writeToVTKSeparate(const std::string& fileName, const GridView& gv) const
            {
#if HAVE_DUNE_FUNCTIONS
              Dune::VTKWriter<GridView> vtkWriter(gv,Dune::VTK::conforming);
              {
                auto f = Dune::Functions::makeAnalyticGridViewFunction([&](auto x)
                    {typename Traits::RangeType output; stochasticPart.evaluate(x,output); return output;},gv);
                vtkWriter.addCellData(f,VTK::FieldInfo("stochastic",VTK::FieldInfo::Type::scalar,1));
              }
              for (unsigned int i = 0; i < trendPart.size(); i++)
              {
                const TrendComponent<Traits>& component = trendPart.getComponent(i);
                auto f = Dune::Functions::makeAnalyticGridViewFunction([&](auto x)
                    {typename Traits::RangeType output; component.evaluate(x,output); return output;},gv);
                vtkWriter.addCellData(f,VTK::FieldInfo(component.name(),VTK::FieldInfo::Type::scalar,1));
              }
              vtkWriter.pwrite(fileName,"","",Dune::VTK::appendedraw);
#else //HAVE_DUNE_FUNCTIONS
              DUNE_THROW(Dune::NotImplemented,"Unstructured VTK output requires dune-grid and dune-functions");
#endif //HAVE_DUNE_FUNCTIONS
            }

          /**
           * @brief Export random field as flat Legacy VTK file
           */
          void writeToLegacyVTK(const std::string& fileName) const
          {
            if ((*traits).commSize > 1)
              DUNE_THROW(Dune::NotImplemented,"legacy VTK output doesn't work for parallel runs");

            LegacyVTKWriter<Traits> legacyWriter(config,fileName);
            legacyWriter.writeScalarData("field",*this);
          }

          /**
           * @brief Export random field as separate Legacy VTK entries
           */
          void writeToLegacyVTKSeparate(const std::string& fileName) const
          {
            if ((*traits).commSize > 1)
              DUNE_THROW(Dune::NotImplemented,"legacy VTK output doesn't work for parallel runs");

            LegacyVTKWriter<Traits> legacyWriter(config,fileName);
            legacyWriter.writeScalarData("stochastic",stochasticPart);
            for (unsigned int i = 0; i < trendPart.size(); i++)
            {
              const TrendComponent<Traits>& component = trendPart.getComponent(i);
              legacyWriter.writeScalarData(component.name(),component);
            }
          }

          /**
           * @brief Make random field homogeneous
           */
          void zero()
          {
            trendPart     .zero();
            stochasticPart.zero();

            if (storeInvMat)
            {
              (*invMatPart).zero();
              invMatValid = true;
            }

            if (storeInvRoot)
            {
              (*invRootPart).zero();
              invRootValid = true;
            }
          }

          /**
           * @brief Double spatial resolution of covariance matrix
           */
          void refineMatrix()
          {
            (*traits).refine();
            (*matrix).update();
          }

          /**
           * @brief Double spatial resolution of random field
           */
          void refine()
          {
            if (storeInvMat && invMatValid)
            {
              (*invMatPart).refine();
              stochasticPart = (*matrix) * (*invMatPart);

              if ((*traits).dim == 3)
              {
                stochasticPart *= 1./8.;
                *invMatPart    *= 1./8.;
              }
              else if ((*traits).dim == 2)
              {
                stochasticPart *= 1./4.;
                *invMatPart    *= 1./4.;
              }
              else if ((*traits).dim == 1)
              {
                stochasticPart *= 1./2.;
                *invMatPart    *= 1./2.;
              }
              else
                DUNE_THROW(Dune::Exception,"dimension of field has to be 1, 2 or 3");

              if (storeInvRoot)
              {
                *invRootPart = (*matrix).multiplyRoot(*invMatPart);

                if ((*traits).dim == 3)
                  *invRootPart *= 1./8.;
                else if ((*traits).dim == 2)
                  *invRootPart *= 1./4.;
                else if ((*traits).dim == 1)
                  *invRootPart *= 1./2.;
                else
                  DUNE_THROW(Dune::Exception,"dimension of field has to be 1, 2 or 3");

                invRootValid = true;
              }

            }
            else if (storeInvRoot && invRootValid)
            {
              (*invRootPart).refine();
              stochasticPart = (*matrix).multiplyRoot(*invRootPart);

              if ((*traits).dim == 3)
              {
                stochasticPart *= 1./8.;
                *invRootPart   *= 1./8.;
              }
              else if ((*traits).dim == 2)
              {
                stochasticPart *= 1./4.;
                *invRootPart   *= 1./4.;
              }
              else if ((*traits).dim == 1)
              {
                stochasticPart *= 1./2.;
                *invRootPart   *= 1./2.;
              }
              else
                DUNE_THROW(Dune::Exception,"dimension of field has to be 1, 2 or 3");

              if (storeInvMat)
              {
                *invMatPart = stochasticPart;
                invMatValid = false;
              }
            }
            else
            {
              stochasticPart.refine();

              if (storeInvMat)
                (*invMatPart).refine();

              if (storeInvRoot)
                (*invRootPart).refine();
            }
          }

          /**
           * @brief Addition assignment operator
           */
          RandomField& operator+=(const RandomField& other)
          {
            trendPart      += other.trendPart;
            stochasticPart += other.stochasticPart;

            if (storeInvMat)
            {
              *invMatPart += *(other.invMatPart);
              invMatValid = invMatValid && other.invMatValid;
            }

            if (storeInvRoot)
            {
              *invRootPart += *(other.invRootPart);
              invRootValid = invRootValid && other.invRootValid;
            }

            return *this;
          }

          /**
           * @brief Subtraction assignment operator
           */
          RandomField& operator-=(const RandomField& other)
          {
            trendPart      -= other.trendPart;
            stochasticPart -= other.stochasticPart;

            if (storeInvMat)
            {
              *invMatPart -= *(other.invMatPart);
              invMatValid = invMatValid && other.invMatValid;
            }

            if (storeInvRoot)
            {
              *invRootPart -= *(other.invRootPart);
              invRootValid = invRootValid && other.invRootValid;
            }

            return *this;
          }

          /**
           * @brief Multiplication with scalar
           */
          RandomField& operator*=(const RF alpha)
          {
            trendPart      *= alpha;
            stochasticPart *= alpha;

            if (storeInvMat)
              *invMatPart *= alpha;

            if (storeInvRoot)
              *invRootPart *= alpha;

            return *this;
          }

          /**
           * @brief AXPY scaled addition
           */
          RandomField& axpy(const RandomField& other, const RF alpha)
          {
            trendPart     .axpy(other.trendPart     ,alpha);
            stochasticPart.axpy(other.stochasticPart,alpha);

            if (storeInvMat)
            {
              (*invMatPart).axpy(*(other.invMatPart),alpha);
              invMatValid = invMatValid && other.invMatValid;
            }

            if (storeInvRoot)
            {
              (*invRootPart).axpy(*(other.invRootPart),alpha);
              invRootValid = invRootValid && other.invRootValid;
            }

            return *this;
          }

          /**
           * @brief Scalar product
           */
          RF operator*(const RandomField& other) const
          {
            RF output = 0.;

            output += (*this).stochasticPart * other.stochasticPart;
            output += (*this).trendPart * other.trendPart;

            return output;
          }

          /**
           * @brief Multiply random field with covariance matrix
           */
          void timesMatrix()
          {
            if (storeInvMat)
            {
              *invMatPart = stochasticPart;
              invMatValid  = true;
            }

            if (storeInvRoot)
            {
              *invRootPart = (*matrix).multiplyRoot(stochasticPart);
              invRootValid = true;
            }

            stochasticPart = (*matrix) * stochasticPart;

            trendPart.timesMatrix();
          }

          /**
           * @brief Multiply random field with inverse of covariance matrix
           */
          void timesInverseMatrix()
          {
            if (storeInvMat && invMatValid)
            {
              if (storeInvRoot)
              {
                *invRootPart = (*matrix).multiplyRoot(*invMatPart);
                invRootValid = true;
              }

              stochasticPart = *invMatPart;
              invMatValid = false;
            }
            else
            {
              stochasticPart = (*matrix).multiplyInverse(stochasticPart);

              if (storeInvMat)
                invMatValid = false;

              if (storeInvRoot)
                invRootValid = false;
            }

            trendPart.timesInverseMatrix();
          }

          /**
           * @brief Multiply random field with approximate root of cov. matrix
           */
          void timesMatrixRoot()
          {
            if (storeInvMat && storeInvRoot)
            {
              *invMatPart = *invRootPart;
              invMatValid = invRootValid;
            }

            if (storeInvRoot)
            {
              *invRootPart = stochasticPart;
              invRootValid = true;
            }

            stochasticPart = (*matrix).multiplyRoot(stochasticPart);

            trendPart.timesMatrixRoot();
          }

          /**
           * @brief Multiply random field with approximate inverse root of cov. matrix
           */
          void timesInvMatRoot()
          {
            if (storeInvRoot && invRootValid)
            {
              stochasticPart = *invRootPart;
              invRootValid = false;

              if (storeInvMat)
              {
                *invRootPart = *invMatPart;
                invRootValid = invMatValid;
                invMatValid  = false;
              }
            }
            else
            {
              stochasticPart = (*matrix).multiplyInverse(stochasticPart);

              if (storeInvRoot)
              {
                *invRootPart = stochasticPart;
                invRootValid = true;
              }

              stochasticPart = (*matrix).multiplyRoot(stochasticPart);

              if (storeInvMat)
                invMatValid = false;
            }

            trendPart.timesInvMatRoot();
          }

          void localize(const typename Traits::DomainType& center, const RF radius)
          {
            stochasticPart.localize(center,radius);

            if (storeInvMat)
              invMatValid = false;

            if (storeInvRoot)
              invRootValid = false;
          }
      };

    /**
     * @brief List of Gaussian random fields in 1D, 2D or 3D
     */
    template<typename GridTraits, bool storeInvMat = true, bool storeInvRoot = false,
      template<typename, bool, bool> class RandomField = Dune::RandomField::RandomField>
      class RandomFieldList
      {
        public:

          typedef RandomField<GridTraits, storeInvMat, storeInvRoot> SubRandomField;

        protected:

          // to allow reading in constructor
          class ParamTreeHelper
          {
            Dune::ParameterTree config;

            public:

            ParamTreeHelper(const std::string& fileName = "")
            {
              if (fileName != "")
              {
                Dune::ParameterTreeParser parser;
                parser.readINITree(fileName+".fieldList",config);
              }
            }

            const Dune::ParameterTree& get() const
            {
              return config;
            }
          };

          const ParamTreeHelper     treeHelper;
          const Dune::ParameterTree config;
          std::vector<std::string>  fieldNames;
          std::vector<std::string>  activeTypes;

          std::map<std::string, std::shared_ptr<SubRandomField> > list;
          std::shared_ptr<SubRandomField> emptyPointer;

          typedef typename GridTraits::RangeField RF;

        public:

          /**
           * @brief Constructor reading random fields from file
           */
          template<typename LoadBalance = DefaultLoadBalance<GridTraits::dim> >
            RandomFieldList(
                const Dune::ParameterTree& config_,
                const std::string& fileName = "",
                const LoadBalance loadBalance = LoadBalance(),
                const MPI_Comm comm = MPI_COMM_WORLD
                )
            : config(config_)
            {
              std::stringstream typeStream(config.get<std::string>("randomField.types"));
              std::string type;
              while(std::getline(typeStream, type, ' '))
              {
                fieldNames.push_back(type);

                Dune::ParameterTree subConfig;
                Dune::ParameterTreeParser parser;
                parser.readINITree(type+".field",subConfig);

                // copy general keys to subConfig if necessary
                if (!subConfig.hasKey("grid.extensions") && config.hasKey("grid.extensions"))
                  subConfig["grid.extensions"] = config["grid.extensions"];
                if (!subConfig.hasKey("grid.cells") && config.hasKey("grid.cells"))
                  subConfig["grid.cells"] = config["grid.cells"];
                if (!subConfig.hasKey("randomField.cgIterations")
                    && config.hasKey("randomField.cgIterations"))
                  subConfig["randomField.cgIterations"] = config["randomField.cgIterations"];

                std::string subFileName = fileName;
                if (subFileName != "")
                  subFileName += "." + type;

                list.insert({type,
                    std::make_shared<SubRandomField>(subConfig,subFileName,loadBalance,comm)});
              }

              if (fieldNames.empty())
                DUNE_THROW(Dune::Exception,"List of randomField types is empty");

              activateFields(config.get<int>("randomField.active",fieldNames.size()));
            }

          /**
           * @brief Constructor reading random fields and config from file
           */
          template<typename LoadBalance = DefaultLoadBalance<GridTraits::dim> >
            RandomFieldList(
                const std::string& fileName,
                const LoadBalance loadBalance = LoadBalance(),
                const MPI_Comm comm = MPI_COMM_WORLD
                )
            : treeHelper(fileName), config(treeHelper.get())
            {
              std::stringstream typeStream(config.get<std::string>("randomField.types"));
              std::string type;
              while(std::getline(typeStream, type, ' '))
              {
                fieldNames.push_back(type);

                std::string subFileName = fileName + "." + type;

                list.insert({type, std::make_shared<SubRandomField>(subFileName,loadBalance,comm)});
              }

              if (fieldNames.empty())
                DUNE_THROW(Dune::Exception,"List of randomField types is empty");

              activateFields(config.get<int>("randomField.active",fieldNames.size()));
            }

          /**
           * @brief Constructor reading random fields from file, but reusing covariance matrices
           */
          RandomFieldList(const RandomFieldList& other, const std::string& fileName)
            : fieldNames(other.fieldNames), activeTypes(other.activeTypes)
          {
            for(const std::pair<std::string,std::shared_ptr<SubRandomField> >& pair : other.list)
              list.insert({pair.first,
                  std::make_shared<SubRandomField>(*(pair.second),fileName + "." + pair.first)});
          }

#if HAVE_DUNE_PDELAB
           /**
           * @brief Constructor converting from GridFunctionSpace and GridVector
           */
          template<typename GFS, typename FieldList>
            RandomFieldList(const RandomFieldList& other, const GFS& gfs, const FieldList& fieldList)
            : fieldNames(other.fieldNames), activeTypes(other.activeTypes)
            {
              for (const std::string& type : activeTypes)
              {
                if (fieldList.find(type) == fieldList.end())
                  DUNE_THROW(Dune::Exception,"Field name " + type + " not found in grid function list");

                std::shared_ptr<SubRandomField> otherField = other.list.find(type)->second;
                list.insert({type,
                    std::make_shared<SubRandomField>(*otherField,gfs,*(fieldList.find(type)->second))});
              }

              for (const std::string& type : fieldNames)
                if (fieldList.find(type) == fieldList.end())
                  list.insert({type,
                      std::make_shared<SubRandomField>(*(other.list.find(type)->second))});
            }
#endif // HAVE_DUNE_PDELAB

          /**
           * @brief Copy constructor
           */
          RandomFieldList(const RandomFieldList& other)
            : fieldNames(other.fieldNames), activeTypes(other.activeTypes)
          {
            for(const std::pair<std::string,std::shared_ptr<SubRandomField> >& pair : other.list)
              list.insert({pair.first, std::make_shared<SubRandomField>(*(pair.second))});
          }

          /**
           * @brief Assignment operator
           */
          RandomFieldList& operator=(const RandomFieldList& other)
          {
            fieldNames  = other.fieldNames;
            activeTypes = other.activeTypes;

            list.clear();
            for(const std::pair<std::string,std::shared_ptr<SubRandomField> >& pair : other.list)
              list.insert({pair.first, std::make_shared<SubRandomField>(*(pair.second))});

            return *this;
          }

          /**
           * @brief Define subset of fields kept constant (i.e. not changed by calculus operators)
           */
          void activateFields(const unsigned int number)
          {
            if (number > fieldNames.size())
              DUNE_THROW(Dune::Exception,"Too many randomFields activated");

            activeTypes.clear();
            for (unsigned int i = 0; i < number; i++)
              activeTypes.push_back(fieldNames[i]);
          }

          /**
           * @brief Generate fields with the desired correlation structure
           */
          void generate(bool allowNonWorldComm = false)
          {
            for(const std::string& type : fieldNames)
              list.find(type)->second->generate(allowNonWorldComm);
          }

          /**
           * @brief Generate fields without correlation structure (i.e. noise)
           */
          void generateUncorrelated(bool allowNonWorldComm = false)
          {
            for(const std::string& type : fieldNames)
              list.find(type)->second->generateUncorrelated(allowNonWorldComm);
          }

          /**
           * @brief Vector of random field types currently active
           */
          const std::vector<std::string> types() const
          {
            return activeTypes;
          }

          /**
           * @brief Access to individual random field
           */
          const std::shared_ptr<SubRandomField>& get(const std::string& type) const
          {
            if (list.find(type) != list.end())
              return (list.find(type))->second;

            return emptyPointer;
          }

          /**
           * @brief Export random fields to files on disk
           */
          void writeToFile(const std::string& fileName) const
          {
            for(const std::string& type : fieldNames)
              list.find(type)->second->writeToFile(fileName + "." + type);

            std::ofstream file(fileName+".fieldList",std::ofstream::trunc);
            config.report(file);
          }

          /**
           * @brief Export random fields as flat unstructured VTK files, requires dune-grid and dune-functions
           */
          template<typename GridView>
            void writeToVTK(const std::string& fileName, const GridView& gv) const
            {
#if HAVE_DUNE_FUNCTIONS
              for (const std::string& type : fieldNames)
                list.find(type)->second->writeToVTK(fileName + "." + type,gv);
#else //HAVE_DUNE_FUNCTIONS
              DUNE_THROW(Dune::NotImplemented,"Unstructured VTK output requires dune-grid and dune-functions");
#endif //HAVE_DUNE_FUNCTIONS
            }

          /**
           * @brief Export random fields as unstructured VTK files, requires dune-grid and dune-functions
           */
          template<typename GridView>
            void writeToVTKSeparate(const std::string& fileName, const GridView& gv) const
            {
#if HAVE_DUNE_FUNCTIONS
              for (const std::string& type : fieldNames)
                list.find(type)->second->writeToVTKSeparate(fileName + "." + type,gv);
#else //HAVE_DUNE_FUNCTIONS
              DUNE_THROW(Dune::NotImplemented,"Unstructured VTK output requires dune-grid and dune-functions");
#endif //HAVE_DUNE_FUNCTIONS
            }

          /**
           * @brief Export random fields as flat Legacy VTK files
           */
          void writeToLegacyVTK(const std::string& fileName) const
          {
            for (const std::string& type : fieldNames)
              list.find(type)->second->writeToLegacyVTK(fileName + "." + type);
          }

          /**
           * @brief Export random fields as separate Legacy VTK entries
           */
          void writeToLegacyVTKSeparate(const std::string& fileName) const
          {
            for (const std::string& type : fieldNames)
              list.find(type)->second->writeToLegacyVTKSeparate(fileName + "." + type);
          }

          /**
           * @brief Set the random fields to zero
           */
          void zero()
          {
            for(const std::string& type : activeTypes)
              list.find(type)->second->zero();
          }

          /**
           * @brief Double spatial resolution of covariance matrix
           */
          void refineMatrix()
          {
            for(const std::string& type : activeTypes)
              list.find(type)->second->refineMatrix();
          }

          /**
           * @brief Double spatial resolution of random fields
           */
          void refine()
          {
            for(const std::string& type : activeTypes)
              list.find(type)->second->refine();
          }

          /**
           * @brief Addition assignment operator
           */
          RandomFieldList& operator+=(const RandomFieldList& other)
          {
            for(const std::string& type : activeTypes)
            {
              if (other.list.find(type) == other.list.end())
                DUNE_THROW(Dune::Exception,"RandomFieldLists don't match in operator+=");

              list.find(type)->second->operator+=(*(other.list.find(type)->second));
            }

            return *this;
          }

          /**
           * @brief Subtraction assignment operator
           */
          RandomFieldList& operator-=(const RandomFieldList& other)
          {
            for(const std::string& type : activeTypes)
            {
              if (other.list.find(type) == other.list.end())
                DUNE_THROW(Dune::Exception,"RandomFieldLists don't match in operator+=");

              list.find(type)->second->operator-=(*(other.list.find(type)->second));
            }

            return *this;
          }

          /**
           * @brief Multiplication with scalar
           */
          RandomFieldList& operator*=(const RF alpha)
          {
            for(const std::string& type : activeTypes)
              list.find(type)->second->operator*=(alpha);

            return *this;
          }

          /**
           * @brief AXPY scaled addition
           */
          RandomFieldList& axpy(const RandomFieldList& other, const RF alpha)
          {
            for(const std::string& type : activeTypes)
            {
              if (other.list.find(type) == other.list.end())
                DUNE_THROW(Dune::Exception,"RandomFieldLists don't match in axpy");

              list.find(type)->second->axpy(*(other.list.find(type)->second),alpha);
            }

            return *this;
          }

          /**
           * @brief Scalar product
           */
          RF operator*(const RandomFieldList& other) const
          {
            RF output = 0.;

            for(const std::string& type : activeTypes)
            {
              if (other.list.find(type) == other.list.end())
                DUNE_THROW(Dune::Exception,"RandomFieldLists don't match in operator*");

              output += list.find(type)->second->operator*(*(other.list.find(type)->second));
            }

            return output;
          }

          /**
           * @brief Multiply random fields with covariance matrix
           */
          void timesMatrix()
          {
            for(const std::string& type : activeTypes)
              list.find(type)->second->timesMatrix();
          }

          /**
           * @brief Multiply random fields with inverse of covariance matrix
           */
          void timesInverseMatrix()
          {
            for(const std::string& type : activeTypes)
              list.find(type)->second->timesInverseMatrix();
          }

          /**
           * @brief Multiply random fields with approximate root of cov. matrix
           */
          void timesMatrixRoot()
          {
            for(const std::string& type : activeTypes)
              list.find(type)->second->timesMatrixRoot();
          }

          /**
           * @brief Multiply random fields with approximate inverse root of cov. matrix
           */
          void timesInvMatRoot()
          {
            for(const std::string& type : activeTypes)
              list.find(type)->second->timesInvMatRoot();
          }

          void localize(const typename GridTraits::Domain& center, const RF radius)
          {
            for(const std::string& type : activeTypes)
              list.find(type)->second->localize(center,radius);
          }

      };
  }
}

#endif // DUNE_RANDOMFIELD_RANDOMFIELD_HH
