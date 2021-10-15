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
#include <ringpcl/export.h>

namespace SCAR
{
    /*!
     * ValueComparator is a class used by the std::sort() algorithm to sort a std::vector from values stored in templated container T
     */
    template <class T>
    class ringpcl_api ValueComparator
    {
    public:
        ValueComparator(const T &values)
            : values_(values)
        {
        }

        bool operator()(index_t id1, index_t id2) const
        {
            return (values_[id1] < values_[id2]);
        }

    private:
        const T &values_;
    };

    //    static double compute_log_sampling( double min, double max, index_t nb_sampling )
    //    {
    //        ringpcl_assert( min > max || min != 0.);
    //        return ( log10( max / min ) ) / static_cast<double> (nb_sampling ) ;
    //    }
    //
    //    static std::vector<std::pair<double, double> > build_log_sampling_plot( double min, double max, index_t nb_sampling )
    //    {
    //        ringpcl_assert( min > max || min != 0. );
    //        double delta = compute_log_sampling( min, max, nb_sampling );
    //
    //        std::vector<std::pair<double, double> > results( nb_sampling );
    //        for( RINGpcl::index_t gitr = 0; gitr < nb_sampling; ++gitr ) {
    //            results[gitr].first = pow( 10, log10( min ) + ( double ) gitr * delta );
    //        }
    //
    //        return results ;
    //    }
    //    static size_t compute_nb_paire_of_point( size_t nb_point)
    //    {
    //        ringpcl_assert( nb_point != 0 );
    //        return static_cast<size_t>( nb_point *  ( nb_point - 1 ) / 2. );
    //    }

} // namespace
