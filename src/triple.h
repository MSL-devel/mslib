/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries) 
 Copyright (C) 2008-2012 The MSL Developer Group (see README.TXT)
 MSL Libraries: http://msl-libraries.org

If used in a scientific publication, please cite: 
 Kulp DW, Subramaniam S, Donald JE, Hannigan BT, Mueller BK, Grigoryan G and 
 Senes A "Structural informatics, modeling and design with a open source 
 Molecular Software Library (MSL)" (2012) J. Comput. Chem, 33, 1645-61 
 DOI: 10.1002/jcc.22968

This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, 
 USA, or go to http://www.gnu.org/copyleft/lesser.txt.
----------------------------------------------------------------------------
*/

#ifndef __GLIBCPP_INTERNAL_TRIPLE_H
#define __GLIBCPP_INTERNAL_TRIPLE_H
// Code from http://www.bioinf.uni-leipzig.de/Software/bbq/doc/triple_8h-source.html
// by Axel Mosig, Bioinformatics, University of Leipzig.

namespace std {
    template <class _T1, class _T2, class _T3>
    struct triple {
        typedef _T1 first_type;
        typedef _T2 second_type;
        typedef _T3 third_type;

        _T1 first;
        _T2 second;
        _T3 third;

        triple() : first(_T1()), second(_T2()), third (_T3()) {
        }

        triple(const _T1& __a, const _T2& __b, const _T3 & __c) : first(__a), second(__b), third(__c) {
        }

        template <class _U1, class _U2, class _U3 >
        triple(const triple<_U1, _U2, _U3>& __p) : first(__p.first), second(__p.second), third(__p.third) {
        }
    };

    template <class _T1, class _T2, class _T3>
    inline bool operator==(const triple<_T1, _T2, _T3>& __x, const triple<_T1, _T2, _T3>& __y) {
        return __x.first == __y.first && __x.second == __y.second && __x.third == __y.third;
    }

    template <class _T1, class _T2, class _T3>
    inline bool operator<(const triple<_T1, _T2, _T3>& __x, const triple<_T1, _T2, _T3>& __y) {
        return __x.first < __y.first ||
                (!(__y.first < __x.first) && __x.second < __y.second) ||
                ((!(__y.first < __x.first) && !(__y.second < __x.second) && __x.third < __y.third));
    }

    template <class _T1, class _T2, class _T3>
    inline bool operator!=(const triple<_T1, _T2, _T3>& __x, const triple<_T1, _T2, _T3>& __y) {
        return !(__x == __y);
    }

    template <class _T1, class _T2, class _T3>
    inline bool operator>(const triple<_T1, _T2, _T3>& __x, const triple<_T1, _T2, _T3>& __y) {
        return __y < __x;
    }

    template <class _T1, class _T2, class _T3>
    inline bool operator<=(const triple<_T1, _T2, _T3>& __x, const triple<_T1, _T2, _T3>& __y) {
        return !(__y < __x);
    }

    template <class _T1, class _T2, class _T3>
    inline bool operator>=(const triple<_T1, _T2, _T3>& __x, const triple<_T1, _T2, _T3>& __y) {
        return !(__x < __y);
    }

    template <class _T1, class _T2, class _T3>
    inline triple<_T1, _T2, _T3> make_triple(const _T1& __x, const _T2& __y, const _T3& __z)
    {
        return triple<_T1, _T2, _T3 > (__x, __y, __z);
    }

} // namespace std

#endif /* __GLIBCPP_INTERNAL_TRIPLE_H */
