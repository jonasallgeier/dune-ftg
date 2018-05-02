#ifndef DUNE_MODELLING_BOUNDARY_HH
#define	DUNE_MODELLING_BOUNDARY_HH

#include<dune/common/parametertreeparser.hh>

#include<dune/modelling/declarations.hh>

namespace Dune {
  namespace Modelling {

    template<typename Traits, typename ModelType, typename DirectionType>
      class Boundary
      {};

    /**
     * @brief Default implementation for boundary conditions
     */
    template<typename Traits, typename ModelType>
      class Boundary<Traits,ModelType,Direction::Forward>
      {
        using GridTraits = typename Traits::GridTraits;

        using RF           = typename GridTraits::RangeField;
        using Domain       = typename GridTraits::Domain;
        using IDomain      = typename GridTraits::IDomain;
        using Intersection = typename GridTraits::Intersection;

        enum {dim = GridTraits::dim};

        //const std::vector<RF> extensions;
        const std::vector<RF> corner_ll;
        const std::vector<RF> corner_ur;
        const RF eps;

        /**
         * @brief Internal class for boundary segments
         */
        class BoundarySegment
        {
          private:

            typedef std::list<std::pair<RF,std::pair<unsigned int,RF> > > ConditionList;
            typedef std::pair<RF,std::pair<unsigned int,RF> >             ConditionPair;

            std::vector<RF> limitsVector;
            RF firstTime, lastTime;
            mutable typename ConditionList::const_iterator lower, upper;
            ConditionList conditionList;

          public:

            BoundarySegment()
              : firstTime(std::numeric_limits<RF>::max()), lastTime(-firstTime)
            {}

            std::vector<RF>& limits()
            {
              return limitsVector;
            }

            const std::vector<RF>& limits() const
            {
              return limitsVector;
            }

            void storeCondition(const std::vector<RF>& conditionVector)
            {
              if (conditionVector.size() != 3)
                DUNE_THROW(Dune::Exception,"syntax error in boundary condition");

              const RF time = conditionVector[0];
              const std::pair<unsigned int,RF> condPair(std::round(conditionVector[1]),conditionVector[2]);

              if (time > lastTime)
              {
                conditionList.push_back(ConditionPair(time,condPair));
                lastTime = time;
              }
              else if (time < firstTime)
              {
                conditionList.push_front(ConditionPair(time,condPair));
                firstTime = time;
              }
              else
              {
                typename ConditionList::iterator it = conditionList.begin();
                while (it != conditionList.end() && it->first < time)
                  ++it;
                conditionList.insert(it, ConditionPair(time,condPair));
              }
            }

            BoundaryCondition::Type bc(RF time) const
            {
              setCondition(time);

              if (upper->second.first == 0)
                return BoundaryCondition::Dirichlet;
              if (upper->second.first == 1)
                return BoundaryCondition::Neumann;
              if (upper->second.first == 2)
                return BoundaryCondition::LimitedFlux;

              return BoundaryCondition::Other; // unknown
            }

            RF g(RF time) const
            {
              setCondition(time);

              return upper->second.second;
            }

            RF j(RF time) const
            {
              setCondition(time);

              return upper->second.second;
            }

          private:

            void setCondition(RF time) const
            {
              typename ConditionList::const_iterator it  = conditionList.begin();
              typename ConditionList::const_iterator it2 = conditionList.begin();

              if (it == conditionList.end())
                DUNE_THROW(Dune::Exception,"BoundarySegment was empty");

              while (it != conditionList.end() && it->first < time - 1e-6)
                it2 = it++;

              if (it == conditionList.end())
                it = it2;

              lower = it2;
              upper = it;
            }

        };

        std::vector<std::vector<BoundarySegment> > segmentVector;

        mutable RF storedTime;

        public:

        Boundary(const Traits& traits, const std::string& name)
          : corner_ll(traits.corners_ll()), corner_ur(traits.corners_ur()),
          eps(1e-5), storedTime(std::numeric_limits<RF>::max())
        {
          Dune::ParameterTree boundaryConfig;
          Dune::ParameterTreeParser parser;
          
	  std::string temporaryname = name;
          temporaryname.erase (temporaryname.begin()+12, temporaryname.end());
          if (temporaryname=="geoelectrics")
          {
            parser.readINITree("boundary.geoelectrics",boundaryConfig);
          }
          else
          {		
            parser.readINITree("boundary."+name,boundaryConfig);
          }
          for (unsigned int i = 0; i < 2*dim; i++)
          {
            segmentVector.push_back(std::vector<BoundarySegment>());

            std::string side;
            if (i == 0)
              side = "left";
            else if (i == 1)
              side = "right";
            else if (dim == 3)
            {
              if (i == 2)
                side = "front";
              else if (i == 3)
                side = "back";
              else if (i == 4)
                side = "bottom";
              else
                side = "top";
            }
            else
            {
              if (i == 2)
                side = "bottom";
              else
                side = "top";
            }

            const unsigned int numSegments = boundaryConfig.get<unsigned int>("segments."+side);

            for (unsigned int j = 0; j < numSegments; j++)
            {
              std::string segmentName;
              std::stringstream s;
              s << side << j;
              s >> segmentName;

              std::vector<RF> limits = boundaryConfig.get<std::vector<RF> >(segmentName+".limits");
              if (limits.size() != 2*(dim-1))
                DUNE_THROW(Dune::Exception,"unable to read limits of boundary segment "+segmentName);

              BoundarySegment boundarySegment;
              readSegment(boundaryConfig,segmentName,boundarySegment);
              segmentVector[i].push_back(boundarySegment);
            }
          }
        }

        /**
         * @brief Access to corresponding forward boundary conditions (here: identity op)
         */
        const Boundary<Traits,ModelType,Direction::Forward>& forward() const
        {
          return *this;
        }

        /**
         * @brief Boundary condition type
         */
        BoundaryCondition::Type bc (const Intersection& is, const IDomain& x, RF time) const
        {
          const Domain& global = is.geometry().global(x);
          unsigned int sideId, segmentId;

          findSegment(global,sideId,segmentId);
          return segmentVector[sideId][segmentId].bc(time);
        }

        /**
         * @brief Dirichlet boundary condition value
         */
        RF g (const Intersection& is, const IDomain& x, RF time) const
        {
          storedTime = time;

          const Domain& global = is.geometry().global(x);
          unsigned int sideId, segmentId;

          findSegment(global,sideId,segmentId);
          return segmentVector[sideId][segmentId].g(time);
        }

        /**
         * @brief Neumann boundary condition value
         */
        RF j (const Intersection& is, const IDomain& x, RF time) const
        {
          storedTime = time;

          const Domain& global = is.geometry().global(x);
          unsigned int sideId, segmentId;

          findSegment(global,sideId,segmentId);
          return segmentVector[sideId][segmentId].j(time);
        }

        private:

        /**
         * @brief Parse boundary segment definition
         */
        void readSegment(
            const Dune::ParameterTree& boundaryConfig,
            const std::string& segmentName,
            BoundarySegment& boundarySegment
            )
        {
          boundarySegment.limits() = boundaryConfig.get<std::vector<RF> >(segmentName+".limits");

          int count = 0;
          std::stringstream s;
          std::vector<RF> emptyVector;
          bool endReached = false;

          RF firstTime =  std::numeric_limits<RF>::max();
          RF lastTime  = -std::numeric_limits<RF>::max();

          while(!endReached)
          {
            s.clear();
            s.str(std::string());
            s << count;
            std::vector<RF> conditionVector
              = boundaryConfig.get<std::vector<RF> >(segmentName+".time"+s.str(),emptyVector);

            if (conditionVector.empty())
            {
              endReached = true;
            }
            else
            {
              boundarySegment.storeCondition(conditionVector);

              count++;
            }
          }
        }

        /**
         * @brief Find boundary segment for coordinate
         */
        void findSegment(const Domain& global, unsigned int& sideId, unsigned int& segmentId) const
        {
          // left boundary
          if (global[0] < eps+corner_ll[0])
            sideId = 0;
          // right boundary
          else if (global[0] > corner_ur[0] - eps)
            sideId = 1;
          else if (dim == 3)
          {
            // front boundary
            if (global[1] < eps+corner_ll[1])
              sideId = 2;
            // back boundary
            else if (global[1] > corner_ur[1] - eps)
              sideId = 3;
            // bottom boundary
            else if (global[2] < eps+corner_ll[2])
              sideId = 4;
            // top boundary
            else
              sideId = 5;
          }
          else
          {
            // bottom boundary
            if (global[1] < eps+corner_ll[1])
              sideId = 2;
            // top boundary
            else
              sideId = 3;
          }

          std::vector<unsigned int> sideDims(2,0);
          for (unsigned int i = 0; i < dim - 1; i++)
            sideDims[i] = (sideId/2 + i + 1)%dim;

          bool found = false;
          for (unsigned int i = 0; i < segmentVector[sideId].size(); i++)
          {
            bool inside = true;
            for (unsigned int j = 0; j < dim - 1; j++)
            {
              const std::vector<RF>& limits = segmentVector[sideId][i].limits();
              if ( global[sideDims[j]] < limits[j] - eps
                  || global[sideDims[j]] > limits[(dim - 1) + j] + eps)
                inside = false;
            }

            if (inside)
            {
              found = true;
              segmentId = i;
              break;
            }
          }

          if (!found)
            DUNE_THROW(Dune::Exception,"boundary segment not found, check limits");
        }

      };

    /**
     * @brief Homogeneous boundary conditions for adjoint equation
     */
    template<typename Traits, typename ModelType>
      class Boundary<Traits,ModelType,Direction::Adjoint>
      {
        using GridTraits = typename Traits::GridTraits;

        using RF           = typename GridTraits::RangeField;
        using IDomain      = typename GridTraits::IDomain;
        using Intersection = typename GridTraits::Intersection;

        const Boundary<Traits,ModelType,Direction::Forward> forwardBoundary;

        public:

        Boundary(const Traits& traits, const std::string& name)
          : forwardBoundary(traits,name)
        {}

        BoundaryCondition::Type bc (const Intersection& is, const IDomain& x, const RF time) const
        {
          return forwardBoundary.bc(is,x,time);
        }

        /**
         * @brief Access to corresponding forward boundary conditions
         */
        const Boundary<Traits,ModelType,Direction::Forward>& forward() const
        {
          return forwardBoundary;
        }

        /**
         * @brief Dirichlet boundary condition
         */
        RF g (const Intersection& is, const IDomain& x, const RF time) const
        {
          return 0.;
        }

        /**
         * @brief Neumann boundary condition
         */
        RF j (const Intersection& is, const IDomain& x, const RF time) const
        {
          return 0.;
        }
      };

  }
}

#endif // DUNE_MODELLING_BOUNDARY_HH
