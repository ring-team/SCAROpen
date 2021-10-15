// Copyright (c) 2021 ASGA and Universite de Lorraine. All Rights Reserved.

// This program (SCAR) was developed in the frame of the RING project managed by
// ASGA and Universite de Lorraine.

// It is distributed under a dual licensing scheme:

// 1. RING Consortium License
// Members of the RING-COCAD Consortium may only use this file in
// accordance with the terms of described in the GOCAD Advancement Agreement,
// without the prior written authorization of the ASGA.
// Licencee agrees to attach or embed this Notice on all copies
// of the program, including partial copies or modified versions thereof.
// Please use: contact at ring dash team dot org, for more information.

// 2. GNU General Public License Usage
// Alternatively, this file may be used under the terms of the
// GNU General Public license version 3. The licenses are as published by
// the Free Software Foundation and appearing in the file
// https://www.gnu.org/licenses/gpl-3.0.html

#pragma once

#include <scar/basic/common.h>

#include <deque>
#include <queue>
#include <numeric>
#include <vector>

#include <geogram/basic/algorithm.h>

#include <ringpcl/utils.h>

// no forward declaration

namespace SCAR
{
    /*! @class CorrelationMap
     * 	@brief Define an adjacency matrix for undirected graph.
     *
     * 	The template parameters can be:
     * 	 - bool: 			defines if an edge link the two vertices
     * 	 - int/char/string: categoric edge definition
     * 	 - double: 			weighted edges
     *
     */
    template <typename T>
    class ringpcl_api CorrelationMap
    {
    public:
        typedef typename std::vector<T> vector_type;
        typedef typename vector_type::reference refT;
        typedef typename vector_type::const_reference const_refT;

    public:
        CorrelationMap(index_t nb_elem, const_refT init_elem)
            : nb_elem_(nb_elem),
              size_(nb_elem_ * (nb_elem_ - 1) / 2),
              map_(nb_elem_ * (nb_elem_ - 1) / 2, init_elem)
        {
        }

        ~CorrelationMap()
        {
        }

        void reset_and_resize(index_t nb_elem, const_refT init_elem)
        {
            nb_elem_ = nb_elem;
            size_ = nb_elem_ * (nb_elem_ - 1) / 2;
            map_ = vector_type(nb_elem_ * (nb_elem_ - 1) / 2, init_elem);
        }

        /*! @return a reference to the value linking the element i and j */
        refT operator()(index_t i, index_t j)
        {
            return map_[do_get_vector_index(i, j)];
        }
        /*! @return a const reference to the value linking the element i and j */
        const_refT operator()(index_t i, index_t j) const
        {
            return map_[do_get_vector_index(i, j)];
        }

        /*! @brief Print the index and the values */
        void print(std::ostream &file) const
        {
            //if( file.is_open() ) {
            file << " \n Correlation Map  \n";
            file << "---  size of correlation map =  " << nb_elem_ << "\n";

            for (index_t i = 0; i < nb_elem_; ++i)
            {
                for (index_t j = 0; j < nb_elem_; ++j)
                {
                    if (i == j)
                    {
                        file << " ( diag ) ";
                    }
                    else
                    {
                        index_t test = do_get_vector_index(i, j);
                        file << " ( " << test << " ; " << map_[test] << " ) ";
                    }
                }
                file << " endline \n \n";
            }
        }

        index_t nb_elem() const
        {
            return nb_elem_;
        }

    private:
        /*!
         * @param i line index
         * @param j column index
         * @return the (i, j) value index in the correlation map
         * @warning valid for for i < j
         */
        index_t get_vector_index(index_t i, index_t j) const
        {
            ringpcl_assert(i < j);
            return (j - (((i + 1) * (i + 2)) / 2) + i * nb_elem_);
        }

        /*!
         * @param i line index
         * @param j column index
         * @return the (i, j) by calling get_vector_index
         * @warning valid if i!=j
         */
        index_t do_get_vector_index(index_t i, index_t j) const
        {
            ringpcl_assert(i < nb_elem_);
            ringpcl_assert(j < nb_elem_);
            ringpcl_assert(i != j);
            index_t e_id = 0;

            if (i < j)
            {
                e_id = get_vector_index(i, j);
            }
            else
            {
                e_id = get_vector_index(j, i);
            }

            if (!(e_id < map_.size()))
            {
                std::cout << " ERROR " << e_id << "  " << map_.size() << std::endl;
                std::cout << " i j size_ " << i << " " << j << " " << size_
                          << std::endl;
            }

            ringpcl_assert(e_id < map_.size());
            return e_id;
        }

    private:
        CorrelationMap();

        /*! @brief Number of elements that should be correlated */
        index_t nb_elem_;
        /*! @brief Size of the map_ (size_*(size_-1)/2) */
        index_t size_;
        /*! @brief Vector of size_ that represent the
         * 	upper-triangular matrix defining the CorrelationMap	 */
        vector_type map_;

        friend class CorrelationMapAPI;
    };

#ifdef WIN32
    // The class has been speecialized to allow dll export ringpcl_api
    // http://stackoverflow.com/questions/17519879/
    template ringpcl_api class CorrelationMap<double>;
    template ringpcl_api class CorrelationMap<bool>;
    template ringpcl_api class CorrelationMap<int>;
    template ringpcl_api class CorrelationMap<unsigned int>;
#endif

    /*! @class CorrelationMapAPI
     *  @brief Implements algorithm using the CorrelationMap data structure
     *  Is friend of CorrelationMap
     */
    class ringpcl_api CorrelationMapAPI
    {
    public:
        /*! @brief Return the number of edges of the given graph
         * (there is an edge if value is not not_connected_value) */
        template <typename T>
        static index_t size(
            const CorrelationMap<T> &map,
            const T &not_connected_value)
        {

            index_t counter = 0;
            for (index_t cur_el = 0; cur_el < map.size_; ++cur_el)
            {
                if (map.map_[cur_el] != not_connected_value)
                {
                    counter++;
                }
            }
            return counter;
        }

        /*!	@brief All the elements connected to the given one
         * (linked by an edge whose value is not not_connected_value)
         * 	are returned */
        template <typename T>
        static std::vector<index_t> directly_linked(
            const CorrelationMap<T> &map,
            index_t elem_index,
            const T &not_connected_value,
            bool add_cur_elem_in_list = true)
        {
            ringpcl_assert(elem_index < map.nb_elem_);

            std::vector<index_t> connected;
            connected.reserve(map.nb_elem_);

            for (index_t cur_el = 0; cur_el < map.nb_elem_; ++cur_el)
            {

                if (cur_el == elem_index)
                {
                    continue;
                }
                else if (map(elem_index, cur_el) != not_connected_value)
                {
                    connected.push_back(cur_el);
                }
            }
            if (add_cur_elem_in_list)
            {
                connected.push_back(elem_index);
            }
            return connected;
        }

        static std::vector<index_t> directly_linked(
            const CorrelationMap<bool> &map,
            index_t elem_index,
            bool add_cur_elem_in_list = true)
        {
            return directly_linked<bool>(map, elem_index, false,
                                         add_cur_elem_in_list);
        }

        /*! @brief labels all connected components
         *  @param[in] adjacency_map: edge-based representation of a graph
         *  @param[out] connected: contains the label of the connected component
         *      for each element
         *  @return the number of connected component
         */
        template <typename T>
        static index_t compute_connected_components_label(
            const CorrelationMap<T> &map,
            const T &not_connected_value,
            std::vector<index_t> &marks)
        {
            // using  Breadth First Search as suggested by
            // http://stackoverflow.com/questions/8124626/

            // allocate mark<int> of length n, where n is the number of vertex in graph and fill it with zeros
            marks = std::vector<index_t>(map.nb_elem_, 0);
            index_t nb_components = 0; //Components = 0;

            for (index_t cur_el = 0; cur_el < map.nb_elem_; ++cur_el)
            { // Enumerate all vertices,
                if (marks[cur_el] == 0)
                {                    // if for vertex number i, marks[i] == 0
                    ++nb_components; // ++Components;

                    std::queue<index_t> cur_queue; // Put this vertex into queue,
                    cur_queue.push(cur_el);

                    while (!cur_queue.empty())
                    {                                             // while queue is not empty
                        index_t cur_neigh_el = cur_queue.front(); // pop vertex v from q
                        cur_queue.pop();
                        marks[cur_neigh_el] = nb_components; //  marks[v] = Components;

                        for (index_t other_neigh_el = 0;
                             other_neigh_el < map.nb_elem_; ++other_neigh_el)
                        {
                            if (cur_neigh_el == other_neigh_el)
                            { // can not be correlated with oneself
                                continue;
                            }
                            if (map(cur_neigh_el, other_neigh_el) == not_connected_value)
                            { // Put all adjacent vertices
                                continue;
                            }
                            if (marks[other_neigh_el] != 0)
                            { //with marks equal to zero into queue.
                                continue;
                            }
                            cur_queue.push(other_neigh_el);
                        }
                    }
                }
            }
            return nb_components;
        }

        /*! @brief labels all connected components
         *  @param[in] adjacency_map: edge-based representation of a graph
         *  @param[out] connected: contains the label of the connected component
         *  	for each element
         *	@return the number of connected component
         */
        static index_t compute_connected_components_label(
            const CorrelationMap<bool> &map,
            std::vector<index_t> &marks);

        /*! @brief identifies elements belonging to the same connected component
         *  @param[in] map: edge-based representation of a graph
         *  @param[out] labels: sets of nodes belonging to the same connected components
         *  @return the number of connected component
         * */
        template <typename T>
        static index_t get_connected_components(
            const CorrelationMap<T> &map,
            std::vector<std::vector<index_t>> &labels,
            const T &not_connected_value,
            bool show = false)
        {
            labels.clear();

            std::vector<index_t> marks;
            index_t nb_component = compute_connected_components_label(map,
                                                                      not_connected_value, marks);

            labels.resize(nb_component);
            for (index_t el = 0; el < map.nb_elem_; ++el)
            {
                ringpcl_assert(marks[el] - 1 < nb_component);
                labels[marks[el] - 1].push_back(el);
            }

            if (show)
            {
                std::cout << "Going to show" << std::endl;
                for (unsigned int cc = 0; cc < labels.size(); ++cc)
                {
                    for (unsigned int el = 0; el < labels[cc].size(); ++el)
                    {
                        std::cout << labels[cc][el] << "\t" << std::flush;
                    }
                    std::cout << std::endl;
                }
            }

            return nb_component;
        }

        /*! @brief identifies elements belonging to the same connected component
         *  @param[in] map: edge-based representation of a graph
         *  @param[out] marks: sets of nodes belonging to the same connected components
         * 	@return the number of connected component
         * */
        static index_t get_connected_components(
            const CorrelationMap<bool> &map,
            std::vector<std::vector<index_t>> &labels,
            bool show = false);
        /*! @brief Add edges so that all connected component is a complete
         * 	subgraph */
        static void set_all_subgraph_complete(CorrelationMap<bool> &map);

        // todo set all subgrap complete

        /*!
         * @brief Get all the indices of correlation value
         * @param[in] i
         * @param[in] j
         * @param[out] correlated_ids between a given point
         * @return the indices that represent a line in the correlation map.
         */
        template <typename T>
        static void point_get_correlated_indexes(
            const CorrelationMap<T> &map,
            index_t i,
            index_t j,
            std::vector<index_t> &correlated_ids)
        {

            index_t point_id = map.do_get_vector_index(i, j);
            ringpcl_assert(point_id < map.size_);

            correlated_ids.resize(map.size_ - 1);
            for (index_t n = 0; n < point_id; ++n)
            {
                correlated_ids[n] = map.do_get_vector_index(n, point_id);
            }

            index_t to_get = map.size_ - (point_id + 1);
            if (to_get > 0)
            {
                index_t cm_ind = map.do_get_vector_index(point_id, point_id + 1);
                for (index_t n = point_id; n < map.size_; ++n)
                {
                    correlated_ids[n] = cm_ind;
                    ++cm_ind;
                }
            }
        }
        template <typename T>
        static void point_get_sorted_correlated_indexes(
            const CorrelationMap<T> &map,
            index_t i,
            index_t j,
            std::vector<index_t> &correlated_ids)
        {
            point_get_correlated_indexes<T>(map, i, j, correlated_ids);
            GEO::sort(correlated_ids.begin(), correlated_ids.end());
            //ValueComparator< std::vector< index_t > >( correlated_ids ) ) ;
        }

        template <typename T>
        static void sorted_correlation_map_indexes(
            const CorrelationMap<T> &map,
            std::vector<index_t> &correlation_map_ids)
        {
            //fill in the correlation_map_ids vector with unsorted indexes
            correlation_map_ids.resize(map.size());
            std::iota(correlation_map_ids.begin(), correlation_map_ids.end(), 0);
            //reimpl_std_iota(correlation_map_ids);
            //sort correlation_map_ids according to correlation map values
            std::sort(correlation_map_ids.begin(), correlation_map_ids.end(),
                      ValueComparator<std::vector<T>>(map.map_));
        }
    };

    /*! @class CentralityMeasureAPI
     * 	@brief Compute centrality measures of CorrelationMap<bool>
     *
     * 	Centrality measures are used to rank nodes according to their topological
     * 	importance. Several centrality measures can be defined.
     *
     * 	CentralityMeasureAPI allows to compute centrality measure on undirected graph
     * 	defined using CorrelationMap<bool>, or weighted graph defined using a
     * 	CorrelationMap<double>
     *
     * 	Centrality measures are generally discussed for oriented graph, but can also
     * 	be extended to undericted graph [Perra et al, 2008]
     *
     * 	If the case of fault / fracture observations correlation, such centrality
     * 	measures can be used to determine the observations that can be correlated
     * 	to a lot of other observations, an thus identify area with a high
     * 	"correlation uncertainty"
     *
     * 	https://en.wikipedia.org/wiki/Centrality
     *
     */
    class ringpcl_api CentralityMeasureAPI : public CorrelationMapAPI
    {
    public:
        /*! @brief Compute the PageRank for each element of a CorrelationMap
         * 	@see Page el al, 1998.  Perra et al, 2008
         * 	@see https://en.wikipedia.org/wiki/PageRank
         *
         * 	@param[in] map:	represents the undericted graph on which we compute the PageRank
         * 	@param[in] damping_factor: avoid to sink in whole for directed graph...
         * 			Commonly set to 0.85
         * 	@param[in] nb_iter:	number of iteration for computation
         *
         * 	@warning: there is a patent on this algorithm. Can we use it?
         * 	@note: there are probably much more efficient implementation
         *
         * 	*/
        static std::vector<float> compute_page_rank(
            const CorrelationMap<bool> &map,
            double damping_factor = 0.85,
            unsigned int nb_iter = 100);

        /*! @brief Computes the degree for each node of an unweighted graph */
        static std::vector<index_t> compute_node_degree(
            const CorrelationMap<bool> &map);

        /*! @brief Computes the degree for each node of a weighted graph */
        static std::vector<double> compute_node_degree(
            const CorrelationMap<double> &map);

        /*! https://en.wikipedia.org/wiki/Floyd%E2%80%93Warshall_algorithm
         * 	By convention, if input_graph[i][j] == 0, there is no edge
         * 	between the nodes i and j */
        static CorrelationMap<double> compute_betweenness_centrality_FloydWashall(
            const CorrelationMap<double> &input_graph);
    };

    class ringpcl_api MaximalCliquesAPI : public CorrelationMapAPI
    {
    public:
        // TODO there is an implementation (probably better) in bost
        static std::vector<std::vector<unsigned int>> maximal_cliques_listing(
            const CorrelationMap<bool> &map);
    };

    /*! @class ClusterMap
     * 	@brief Data structure to group Elements into different Cluster
     *
     *  Clusters are stored in a std::vector<t>
     * 	The template parameters can be:
     * 	 - bool: 	defines if the element is in the cluster or not
     * 	 - double:	defines the probability that the element is in the cluster or not
     *  The size of the vector depend on both the number of elements and the number of cluster.
     *  |   |  nb elements   |
     *  | n |----------------|
     *  | b |----------------|
     *  |   |----------------|
     *  | c |----------------|
     *  | l |----------------|
     *  | u |----------------|
     *  | s |----------------|
     *  | t |----------------|
     *  | e |----------------|
     *  | r |----------------|
     *
     *  @note The number of cluster can be adjusted during the clustering.
     *  @note If you want to change the number of elements, You have to create a new ClusterMap.
     */
    template <typename T>
    class ringpcl_api ClusterMap
    {
    public:
        typedef typename std::vector<T> vector_type;
        typedef typename vector_type::reference refT;
        typedef typename vector_type::const_reference const_refT;

        ClusterMap(index_t nb_elem, index_t nb_cluster, const_refT initial_value)
            : nb_elem_(nb_elem),
              nb_cluster_(nb_cluster),
              map_size_(nb_elem * nb_cluster),
              map_(map_size_, initial_value)
        {
            ringpcl_assert(nb_elem_ > 0);
        }

        ~ClusterMap(){};

        /*!
         * @param el_ind: line index representing the element
         * @param cl_ind: column index for the clustering
         * @return a reference to the value linking the element el_ind and cl_ind
         */
        refT operator()(index_t el_ind, index_t cl_ind)
        {
            return map_[index(el_ind, cl_ind)];
        }
        /*!
         * @param el_ind: line index representing the element
         * @param cl_ind: column index for the clustering
         * @return a const reference to the value linking the element el_ind andcl_ind
         */
        const_refT operator()(index_t el_ind, index_t cl_ind) const
        {

            return map_[index(el_ind, cl_ind)];
        }

        index_t nb_cluster() const
        {
            return nb_cluster_;
        }

        index_t nb_elem() const
        {
            return nb_elem_;
        }

    private:
        // Forget about this
        ClusterMap();
        /*! @brief Number of Elements that should be clustered */
        const index_t nb_elem_;
        /*! @brief Pre-defined number of Clusters
         *	A cluster can be empty
         *	One of the Cluster can represents the unassigned Elements
         */
        const index_t nb_cluster_;

        /*! @brief  Size of the map_ : nb_cluster_* nb_elem_ */
        index_t map_size_;
        /*! @brief Vector of nb_cluster_* nb_elem_ defining if an
         * 	Element is in the Cluster of not */
        vector_type map_;

        /*!
         * @param el_ind: line index representing the element
         * @param cl_ind: column index for the clustering
         * @return the (i, j) value index in the cluster map
         */
        index_t index(index_t el_ind, index_t cl_ind) const
        {
            ringpcl_assert(el_ind < nb_elem_);
            ringpcl_assert(cl_ind < nb_cluster_);
            return cl_ind * nb_elem_ + el_ind;
        }

        friend class ClusterMapAPI;
    };

    /*! @class ClusterMapAPI
     *  @brief Implements algorithm using the CorrelationMap data structure
     *
     *  Is friend of ClusterMap
     */
    class ringpcl_api ClusterMapAPI
    {
    public:
        /*!@brief create ClusterMap from a list of cluster indexes
         */
        static ClusterMap<bool> create_cluster_map(
            std::deque<std::deque<index_t>> &lst_clusters,
            index_t total_nb_elements);

        /*! @return true if each element is assigned once and only once to a cluster */
        static bool assigned_only_once(const ClusterMap<bool> &);

        /*!@brief test the assignement of element in a cluster
         * @param[in] el_ind the index of the element to test
         * @return true if the element \param el_ind is assigned to at least one cluster return false if not.
         * @param[out] cl_ind the index of the cluster that hold the element \param el_ind
         * @note if the element is not assign to any cluster the \param cl_ind return the index of the last cluster of the map
         */
        static bool is_assigned_to_cluster(
            const ClusterMap<bool> &,
            index_t el_ind,
            index_t &cl_ind);

        /*! @brief Set a constant value in the cluster map */
        template <typename T>
        static void set_value(ClusterMap<T> &map, const T &val)
        {
            for (index_t ind = 0; ind < map.map_.size(); ++ind)
            {
                map.map_[ind] = val;
            }
        }

        /*! @return Unassigned the Element from old cluster and assigned it to the the Cluster
         *  @warning crashes if the Element is not set in the old_cluster
         */
        static void move_element(
            ClusterMap<bool> &,
            index_t el_ind,
            index_t old_cluster,
            index_t new_cluster);

        /*! @return the index of the first Cluster index the Element belongs to */
        static index_t get_cluster_from_element(
            const ClusterMap<bool> &,
            index_t el_ind);

        /*! @return all the Element indices that belong to the given cluster  */
        static std::vector<index_t> get_elements_from_cluster(
            const ClusterMap<bool> &,
            index_t cl_ind);

        /*! @return a vector of vector of elements index_t that constitutes clusters */
        static std::vector<std::vector<index_t>> get_all_clusters(
            const ClusterMap<bool> &);

        /*! \brief Brief Print the content */
        template <typename T>
        static void print(const ClusterMap<T> &clm, std::ostream &stream)
        {

            stream << "\n --- ClusterMap(i,j,index,value) ---" << std::endl;
            for (unsigned int cl_ind = 0; cl_ind < clm.nb_cluster_; ++cl_ind)
            {
                for (unsigned int el_ind = 0; el_ind < clm.nb_elem_; ++el_ind)
                {
                    stream << "(" << el_ind << "," << cl_ind << ","
                           << clm(el_ind, cl_ind) << ")\t" << std::flush;
                }
                stream << std::endl;
            }
        }
    };

}
//namespace
