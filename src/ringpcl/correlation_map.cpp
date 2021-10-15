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

/*[
 * correlation_map.cpp
 *
 *  Created on: Jun 9, 2016
 *      Author: bonneau
 *
 ]*/

#include <ringpcl/correlation_map.h>
#include <queue>
#include <stack>
#include <list>

namespace SCAR
{

	index_t CorrelationMapAPI::compute_connected_components_label(
		const CorrelationMap<bool> &map,
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
			{					 // if for vertex number i, marks[i] == 0
				++nb_components; // ++Components;

				std::queue<index_t> cur_queue; // Put this vertex into queue,
				cur_queue.push(cur_el);

				while (!cur_queue.empty())
				{											  // while queue is not empty
					index_t cur_neigh_el = cur_queue.front(); // pop vertex v from q
					cur_queue.pop();
					marks[cur_neigh_el] = nb_components; //  marks[v] = Components;

					for (index_t other_neigh_el = 0; other_neigh_el < map.nb_elem_;
						 ++other_neigh_el)
					{
						if (cur_neigh_el == other_neigh_el)
						{ // can not be correlated with oneself
							continue;
						}
						if (!map(cur_neigh_el, other_neigh_el))
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

	index_t CorrelationMapAPI::get_connected_components(
		const CorrelationMap<bool> &map,
		std::vector<std::vector<index_t>> &labels,
		bool show)
	{
		labels.clear();

		std::vector<index_t> marks;
		index_t nb_component = compute_connected_components_label(map, marks);

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

	void CorrelationMapAPI::set_all_subgraph_complete(CorrelationMap<bool> &map)
	{
		std::vector<std::vector<index_t>> labels;
		get_connected_components(map, labels);
		for (unsigned int cc = 0; cc < labels.size(); ++cc)
		{
			for (unsigned int el1 = 0; el1 < labels[cc].size(); ++el1)
			{
				for (unsigned int el2 = el1 + 1; el2 < labels[cc].size(); ++el2)
				{
					map(labels[cc][el1], labels[cc][el2]) = true;
				}
			}
		}
	}

	std::vector<float> CentralityMeasureAPI::compute_page_rank(
		const CorrelationMap<bool> &map,
		double damping_factor,
		unsigned int nb_iter)
	{
		//////////////////////////////////////////////////////////
		//	https://en.wikipedia.org/wiki/PageRank#Computation	//
		//////////////////////////////////////////////////////////

		const index_t N = map.nb_elem();

		// L(p_j) is the number of outbound links on page  p{j}
		std::vector<unsigned int> L(map.nb_elem());
		for (unsigned int i = 0; i < N; ++i)
		{
			for (unsigned int j = 0; j < N; ++j)
			{
				if (i == j)
				{
					continue;
				}

				if (map(i, j))
				{
					L[i] += 1;
				}
			}
		}
		// PageRank value for each elements of the correlation map
		// Initialize to 1./nb_elem for each element
		std::vector<float> page_rank(map.nb_elem(), (float)(1.0 / N));

		std::vector<float> page_rank_tmp(N, 0.0f);
		for (unsigned int iter = 0; iter < nb_iter; ++iter)
		{ // iterative computation

			for (unsigned int elem = 0; elem < N; ++elem)
			{ // re-compute the PageRank for each element

				page_rank_tmp[elem] =
					(float)((1. - damping_factor) / map.nb_elem());

				for (unsigned int neigh = 0; neigh < N; ++neigh)
				{
					if (neigh == elem)
					{
						continue;
					}
					// linked
					if (map(elem, neigh))
					{
						page_rank_tmp[elem] += (float)(damping_factor * (double)page_rank[neigh] / L[neigh]);
					}
				}
			}
			page_rank = page_rank_tmp;
		}

		return page_rank;
	}

	std::vector<index_t> CentralityMeasureAPI::compute_node_degree(
		const CorrelationMap<bool> &map)
	{
		const index_t N = map.nb_elem();
		std::vector<index_t> node_degree(N, 0);
		for (unsigned int i = 0; i < N; ++i)
		{
			for (unsigned int j = i + 1; j < N; ++j)
			{
				if (map(i, j))
				{
					node_degree[i] += 1;
					node_degree[j] += 1;
				}
			}
		}
		return node_degree;
	}

	std::vector<double> CentralityMeasureAPI::compute_node_degree(
		const CorrelationMap<double> &map)
	{
		const unsigned int N = map.nb_elem();
		std::vector<double> node_degree(N, 0.);
		for (unsigned int i = 0; i < N; ++i)
		{
			for (unsigned int j = i + 1; j < N; ++j)
			{
				node_degree[i] += map(i, j);
				node_degree[j] += map(i, j);
			}
		}
		return node_degree;
	}

	CorrelationMap<double> CentralityMeasureAPI::compute_betweenness_centrality_FloydWashall(
		const CorrelationMap<double> &input_graph)
	{
		const unsigned int number_of_elements = input_graph.nb_elem();
		const double no_edge_value = 0.;
		CorrelationMap<double> centrality_graph(number_of_elements, no_edge_value);

		for (unsigned int diag_elem = 0; diag_elem < number_of_elements; ++diag_elem)
		{
			centrality_graph(diag_elem, diag_elem) = 0.;
		}

		for (unsigned int i = 0; i < number_of_elements; ++i)
		{
			for (unsigned int j = i + 1; j < number_of_elements; ++j)
			{
				if (input_graph(i, j) != no_edge_value)
				{
					centrality_graph(i, j) = input_graph(i, j);
				}
			}
		}

		for (unsigned int k = 0; k < number_of_elements; ++k)
		{
			for (unsigned int j = 0; j < number_of_elements; ++j)
			{
				for (unsigned int i = 0; i < number_of_elements; ++i)
				{
					if (centrality_graph(i, j) > centrality_graph(i, k) + centrality_graph(k, j))
					{
						centrality_graph(i, j) = centrality_graph(i, k) + centrality_graph(k, j);
					}
				}
			}
		}
		return centrality_graph;
	}
	//https://en.wikipedia.org/wiki/Bron%E2%80%93Kerbosch_algorithm#Without_pivoting
	static void bron_kerbosch1(
		const CorrelationMap<bool> &map,
		std::set<unsigned int> &R,
		std::set<unsigned int> &P,
		std::set<unsigned int> &X,
		std::vector<std::vector<unsigned int>> &max_cliques)
	{

		if (P.empty() && X.empty())
		{																		  //  if P and X are both empty:
			max_cliques.push_back(std::vector<unsigned int>(R.begin(), R.end())); // report R as a maximal clique
		}

		while (!P.empty())
		{ // for each vertex v in P:
			unsigned int v = *P.begin();
			// BronKerbosch1(R ⋃ {v}, P ⋂ N(v), X ⋂ N(v))
			std::vector<unsigned int> neigh_v = CorrelationMapAPI::directly_linked(map, v, false); // N(v)

			std::set<unsigned int> R_union_v = R; // R ⋃ {v}
			R_union_v.insert(v);

			std::set<unsigned int> P_inter_neigh_v;
			std::set_intersection(P.begin(), P.end(), neigh_v.begin(), neigh_v.end(),
								  std::inserter(P_inter_neigh_v, P_inter_neigh_v.begin()));

			std::set<unsigned int> X_inter_neigh_v;
			std::set_intersection(X.begin(), X.end(), neigh_v.begin(), neigh_v.end(),
								  std::inserter(X_inter_neigh_v, X_inter_neigh_v.begin()));

			bron_kerbosch1(map, R_union_v, P_inter_neigh_v, X_inter_neigh_v, max_cliques);

			// end BronKerbosch1(R ⋃ {v}, P ⋂ N(v), X ⋂ N(v))
			P.erase(v);	 //    P := P \ {v}
			X.insert(v); // X := X ⋃ {v}
		}
	}

	std::vector<std::vector<unsigned int>> MaximalCliquesAPI::maximal_cliques_listing(
		const CorrelationMap<bool> &map)
	{
		std::vector<std::vector<unsigned int>> maximal_cliques;
		std::set<unsigned int> R, P, X;
		for (unsigned int i = 0; i < map.nb_elem(); ++i)
		{
			P.insert(i);
		}
		bron_kerbosch1(map, R, P, X, maximal_cliques);
		return maximal_cliques;
	}

	ClusterMap<bool> ClusterMapAPI::create_cluster_map(std::deque<std::deque<index_t>> &lst_clusters, index_t total_nb_elements)
	{
		ClusterMap<bool> map(total_nb_elements, static_cast<index_t>(lst_clusters.size()), false);
		for (index_t cl = 0; cl < lst_clusters.size(); ++cl)
		{
			for (index_t el = 0; el < lst_clusters[cl].size(); ++el)
			{
				map(lst_clusters[cl][el], cl) = true;
			}
		}
		return map;
	}

	bool ClusterMapAPI::assigned_only_once(const ClusterMap<bool> &map)
	{

		for (index_t el_ind = 0; el_ind < map.nb_elem_; ++el_ind)
		{
			index_t counter = 0;
			for (index_t cl_ind = 0; cl_ind < map.nb_cluster_; ++cl_ind)
			{
				if (map(el_ind, cl_ind) == true)
				{
					++counter;
				}
			}
			if (counter != 1)
			{
				return false;
			}
		}
		return true;
	}

	bool ClusterMapAPI::is_assigned_to_cluster(const ClusterMap<bool> &map, index_t el_ind, index_t &cl_ind)
	{
		cl_ind = 0;
		for (; cl_ind < map.nb_cluster_; ++cl_ind)
		{
			if (map(el_ind, cl_ind) == true)
			{
				return true;
			}
		}
		return false;
	}

	void ClusterMapAPI::move_element(
		ClusterMap<bool> &map,
		index_t el_ind,
		index_t old_cluster,
		index_t new_cluster)
	{

		ringpcl_assert(el_ind < map.nb_elem_);
		ringpcl_assert(old_cluster < map.nb_cluster_);
		ringpcl_assert(new_cluster < map.nb_cluster_);

		ringpcl_assert(map(el_ind, old_cluster) == true);
		ringpcl_assert(map(el_ind, new_cluster) == false);

		map(el_ind, old_cluster) = false;
		map(el_ind, new_cluster) = true;
	}

	index_t ClusterMapAPI::get_cluster_from_element(
		const ClusterMap<bool> &map,
		index_t el_ind)
	{

		for (index_t cl_ind = 0; cl_ind < map.nb_cluster_; ++cl_ind)
		{
			if (map(el_ind, cl_ind) == true)
			{
				return cl_ind;
			}
		}
		ringpcl_assert(0 == 1);
		std::cout << "Error, this element is not assigned to any cluster"
				  << std::endl;
		return 0;
	}

	std::vector<index_t> ClusterMapAPI::get_elements_from_cluster(
		const ClusterMap<bool> &map,
		index_t cl_ind)
	{

		std::vector<index_t> indices;
		indices.reserve(map.nb_elem_);
		for (index_t el_ind = 0; el_ind < map.nb_elem_; ++el_ind)
		{
			if (map(el_ind, cl_ind) == true)
			{
				indices.push_back(el_ind);
			}
		}
		indices.shrink_to_fit();
		return indices;
	}

	std::vector<std::vector<index_t>> ClusterMapAPI::get_all_clusters(const ClusterMap<bool> &map)
	{
		std::vector<std::vector<index_t>> clusters(map.nb_cluster_);
		for (index_t cl_ind = 0; cl_ind < map.nb_cluster_; ++cl_ind)
		{
			clusters[cl_ind] = get_elements_from_cluster(map, cl_ind);
		}
		return clusters;
	}

}
//namespace
