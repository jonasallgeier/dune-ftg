#ifndef DUNE_INVERSION_REORDEREDGRIDVIEW_HH
#define DUNE_INVERSION_REORDEREDGRIDVIEW_HH

#include <list>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/shared_ptr.hh>

#include <dune/grid/common/indexidset.hh>
#include <dune/grid/common/grid.hh> 

#include <dune/pdelab/common/elementmapper.hh>

namespace Dune{
  namespace Modelling{

    struct PressureLikeOrdering
    {
      PressureLikeOrdering() {};

      template <typename RF>
        bool operator()(const RF& value1, const RF& value2) const
        {
          if( value1 < value2 )
            return true;
          else
            return false;
        }
    };

    struct ReversePressureLikeOrdering
    {
      ReversePressureLikeOrdering() {};

      template <typename RF>
        bool operator()(const RF& value1, const RF& value2) const
        {
          if( value1 > value2 )
            return true;
          else
            return false;
        }
    };


    template<typename GridView, typename Parameter>
      class ReorderedIndexSet
      : public Dune::IndexSet<GridView,ReorderedIndexSet<GridView,Parameter>,typename GridView::IndexSet::IndexType>
      {
        typedef typename GridView::Grid     Grid;
        typedef typename GridView::IndexSet Backend;

        typedef typename Parameter::RF_public RF;

        public:

        // Type of base class index type
        typedef typename Backend::IndexType IndexType;

        private:

        // Extract dimension
        enum {dim = Grid::dimension };   

        const GridView&  gridView;
        const Parameter& param;
        const Backend&   backend;

        std::vector<IndexType> reorderedIndicesVertices;
        std::vector<IndexType> reorderedIndicesElements;

        public:

        // called from constructor of class ReorderedGridView:
        ReorderedIndexSet(const GridView& gridView_, const Parameter& param_)
          : gridView(gridView_), param(param_), backend(gridView.indexSet())
        {}

        template<typename Compare>
          void reorder(const Compare& compare)
          {
            generateReorderedIndices<0>  (compare,reorderedIndicesElements);
            generateReorderedIndices<dim>(compare,reorderedIndicesVertices);
          }

        template<typename ScalarDGF, typename Compare>
          void reorder(const ScalarDGF& potential, const Compare& compare)
          {
            generateReorderedIndices<0>  (potential,compare,reorderedIndicesElements);
            generateReorderedIndices<dim>(potential,compare,reorderedIndicesVertices);
          }

        template<int cc>
          IndexType index (const typename std::remove_const<Grid>::type::Traits::template Codim<cc>::Entity& e) const
          {
            return this->index(e);
          }

        // Return index of entity
        template<typename Entity>
          IndexType index (const Entity& e) const
          {
            int cc = Entity::codimension;
            assert( cc == 0 || cc == Backend::dimension );
            if( cc==0 )
            {
              return reorderedIndicesElements[backend.index(e)]; // new index [old index]
            }
            else if( cc==dim )
            {
              return reorderedIndicesVertices[backend.index(e)]; 
            }
            else
            {
              DUNE_THROW(Dune::RangeError,"rgv.index() for codim " << Entity::codimension << " should never be called.");
              return 0;
            }
          }

        template<int cc>
          IndexType subIndex (const typename std::remove_const<Grid>::type::Traits::template Codim<cc>::Entity &e, int i, unsigned int codim ) const
          {
            return this->subIndex(e, i, codim);
          }

        // Return index of subentity
        template<typename Entity>
          IndexType subIndex (const Entity &e, int i, unsigned int codim ) const
          {
            int cc = Entity::codimension;
            assert( cc == 0 || cc == Backend::dimension );
            int cc_codim = cc+codim;
            if ( cc_codim == 0 )
              return reorderedIndicesElements[backend.subIndex(e,i,codim)]; 
            else if (cc_codim == dim) 
              return reorderedIndicesVertices[backend.subIndex(e,i,codim)];
            else {
              DUNE_THROW(Dune::RangeError,"subIndex for entity with codim " << Entity::codimension << " should never be called.");
              return 0;
            }
          }

        const std::vector<Dune::GeometryType>& geomTypes (int codim) const
        {
          return backend.geomTypes(codim);
        }

        IndexType size (Dune::GeometryType type) const
        {
          return backend.size(type);
        }

        IndexType size (int codim) const
        {
          return backend.size(codim);
        }

        template<typename EntityType>
          bool contains (const EntityType& e) const
          {
            return backend.contains(e);
          }

        private:

        template<typename Compare>
          struct PositionIndexCompare
          {
            PositionIndexCompare(const Compare& compare_)
              : compare(compare_)
            {}

            bool operator()(const typename std::pair<RF,IndexType>& v1i, const typename std::pair<RF,IndexType>& v2i)
            {
              return compare(v1i.first,v2i.first);
            }

            Compare compare;
          };

        // Generate the mapping table for a given codim.
        template<int codim, typename Compare>
          void generateReorderedIndices (const Compare& compare, std::vector<IndexType>& reorderedIndices)
          {
            typedef typename std::pair<RF,IndexType> PositionIndex;
            typedef typename std::list<PositionIndex> IndexList;

            // Loop over all entities it and generate a list of pairs 
            // {v(it),index(it)} where v(it) is the center of the entity
            // and index(it) its index 
            IndexList indexList;

            // Map each cell to unique id
            Dune::PDELab::ElementMapper<GridView> cell_mapper(gridView);

            if( codim==dim )
            {
              typedef typename GridView::template Codim<dim>::template Partition<Dune::All_Partition>::Iterator Iterator;
              typedef typename GridView::template Codim<dim>::Geometry::GlobalCoordinate CoordType;

              for (Iterator it=gridView.template begin<dim,Dune::All_Partition>(); it!=gridView.template end<dim,Dune::All_Partition>();++it)
              {
                const CoordType xglobal = it->geometry().corner(0);

                const IndexType idx = gridView.indexSet().index(*it);
                const RF potential  = (RF) idx; // vertices order remain unchanged!
                indexList.push_back( std::make_pair(potential,idx));
              }
            }
            else if( codim==0 )
            {
              typedef typename GridView::Traits::template Codim<0>::template Partition<Dune::All_Partition>::Iterator Iterator;
              typedef typename GridView::Traits::template Codim<0>::Geometry::GlobalCoordinate CoordType;

              for (Iterator it=gridView.template begin<0,Dune::All_Partition>(); it!=gridView.template end<0,Dune::All_Partition>();++it)
              {
                const CoordType xlocal(0.5); // evaluate pressure on the center of the element

                const IndexType idx = gridView.indexSet().index(*it);
                const RF pot = param.potential(*it,xlocal);
                indexList.push_back(std::make_pair(pot,idx));
              }
            }

            PositionIndexCompare<Compare> comparePairs(compare);

            // Sort list according to first entry in each pair
            indexList.sort(comparePairs);
            reorderedIndices.resize(indexList.size());

            int i=0;
            // Generate the mapping table using the index information in the
            // reordered list.
            for(typename IndexList::const_iterator it=indexList.begin(); it!=indexList.end();it++)
            {
              reorderedIndices[it->second] = i;  // new index [old index]
              i++;
            }
          }

        // Generate the mapping table for a given codim.
        template<int codim, typename ScalarDGF, typename Compare>
          void generateReorderedIndices (const ScalarDGF& potential, const Compare& compare, std::vector<IndexType>& reorderedIndices)
          {
            typedef typename std::pair<RF,IndexType> PositionIndex;
            typedef typename std::list<PositionIndex> IndexList;

            // Loop over all entities it and generate a list of pairs 
            // {v(it),index(it)} where v(it) is the center of the entity
            // and index(it) its index 
            IndexList indexList;

            // Map each cell to unique id
            Dune::PDELab::ElementMapper<GridView> cell_mapper(gridView);

            if( codim==dim )
            {
              typedef typename GridView::template Codim<dim>::template Partition<Dune::All_Partition>::Iterator Iterator;
              typedef typename GridView::template Codim<dim>::Geometry::GlobalCoordinate CoordType;

              for (Iterator it=gridView.template begin<dim,Dune::All_Partition>(); it!=gridView.template end<dim,Dune::All_Partition>();++it)
              {
                const CoordType xglobal = it->geometry().corner(0);

                const IndexType idx = gridView.indexSet().index(*it);
                const RF potential  = (RF) idx; // vertices order remain unchanged!
                indexList.push_back( std::make_pair(potential,idx));
              }
            }
            else if( codim==0 )
            {
              typedef typename GridView::Traits::template Codim<0>::template Partition<Dune::All_Partition>::Iterator Iterator;
              typedef typename GridView::Traits::template Codim<0>::Geometry::GlobalCoordinate CoordType;

              for (Iterator it=gridView.template begin<0,Dune::All_Partition>(); it!=gridView.template end<0,Dune::All_Partition>();++it)
              {
                const CoordType xlocal(0.5); // evaluate pressure on the center of the element

                const IndexType idx = gridView.indexSet().index(*it);
                Dune::FieldVector<RF,1> pot;
                potential.evaluate(*it,xlocal,pot);
                indexList.push_back(std::make_pair(pot[0],idx));
              }
            }

            PositionIndexCompare<Compare> comparePairs(compare);

            // Sort list according to first entry in each pair
            indexList.sort(comparePairs);
            reorderedIndices.resize(indexList.size());

            int i=0;
            // Generate the mapping table using the index information in the
            // reordered list.
            for(typename IndexList::const_iterator it=indexList.begin(); it!=indexList.end();it++)
            {
              reorderedIndices[it->second] = i;  // new index [old index]
              i++;
            }
          }

      };


    /* ************************************************************ *
     * Grid view with reordered index set.
     *
     * The ordering of the index set is defined by the Compare class
     * ************************************************************ */
    template<typename GridView, typename Parameter>
      class ReorderedGridView
      : public GridView
      {
        public:

          /* ************************************************************ *
           *  subclass for reordered index set for elements and vertices.
           *  Reordering is realised with the (std::vector) mapping tables 
           *
           *    reorderedindicesElements[] and 
           *    reorderedindicesVertices[]
           *
           *  which return the reordered index for an index returned
           *  by the base class index function. 
           * ************************************************************ */
          // Internal index set with vertical index running fastest
          typedef ReorderedIndexSet<GridView, Parameter> IndexSet;
          typedef typename IndexSet::IndexType IndexType;

          const Parameter& param;
          IndexSet reorderedIndexSet; 

          // Copy (from base class) constructor
          ReorderedGridView(const ReorderedGridView& other)
            : GridView(other.getGridView()), param(other.param), reorderedIndexSet(this->getGridView(), param)
          {}

          // constructor of class ReorderedGridView 
          ReorderedGridView(const GridView& other, const Parameter& param_)
            : GridView(other), param(param_), reorderedIndexSet(this->getGridView(), param)
          {}

          template<typename Compare>
            void reorder(const Compare& compare)
            {
              reorderedIndexSet.reorder(compare);
            }

          template<typename ScalarDGF, typename Compare>
            void reorder(const ScalarDGF& potential, const Compare& compare)
            {
              reorderedIndexSet.reorder(potential,compare);
            }

          const IndexSet& indexSet() const
          {
            return reorderedIndexSet;
          }

          const GridView& getGridView() const
          {
            return *static_cast<const GridView*>(this);
          }

      };

  }
}
#endif
