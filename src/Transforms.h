/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Simulation Library)n
 Copyright (C) 2009 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan

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

#ifndef TRANSFORMS_H
#define TRANSFORMS_H

#include "AtomVector.h"
#include "Residue.h"
#include "Quaternion.h"
#include "SphericalPoint.h"
#include <math.h>

using namespace std;

class Transforms {

	public:
		Transforms();
		Transforms(const Transforms & _transform);
		~Transforms();

		/***************************************************************************************
		 *
		 *  TRANSFORMATIONS APPLIED TO A SINGLE ATOM
		 *
		 ***************************************************************************************/
		// translation
		void translate(Atom & _atom, const CartesianPoint & _p);

		// rotation
		void Xrotate(Atom & _atom, double _degrees);
		void Yrotate(Atom & _atom, double _degrees);
		void Zrotate(Atom & _atom, double _degrees);
		void rotate(Atom & _atom, double _degrees, const CartesianPoint & _axisFromRotCenter, const CartesianPoint & _rotCenter=CartesianPoint(0.0, 0.0, 0.0));
		void rotate(Atom & _atom, const Matrix & _rotMatrix, const CartesianPoint & _rotCenter=CartesianPoint(0.0, 0.0, 0.0));

		/*****************************************************
		 *  Align rotates the atom so that it points in the
		 *  same direction of the target (i.e. the angle
		 *  _atom - _rotCenter - _target is = 0)
		 *
		 *  Orient operates a torsional rotation to orient the
		 *  _atom in the same direction of the _target with
		 *  respect to the axis of rotation defined bv _axis1 - _axis2 
		 *  (i.e. the dihedral(_atom, _axis1, _axis2, _target = 0)
		 *****************************************************/
		void align(Atom & _atom, const CartesianPoint & _target, const CartesianPoint & _rotCenter=CartesianPoint(0.0, 0.0, 0.0));
		void orient(Atom & _atom, const CartesianPoint & _target, const CartesianPoint & _axis1, const CartesianPoint & _axis2);

		// Internal functions of above, but also useful to be public.
		bool align(CartesianPoint & _object, const CartesianPoint & _target, const CartesianPoint & _rotCenter=CartesianPoint(0.0, 0.0, 0.0));
		bool orient(CartesianPoint & _object, const CartesianPoint & _target, const CartesianPoint & _axis1, const CartesianPoint & _axis2);


		/***************************************************************************************
		 *
		 *  TRANSFORMATIONS APPLIED TO AN ATOM VECTOR
		 *
		 ***************************************************************************************/
		void translate(AtomVector & _atoms, CartesianPoint _p);
		void Xrotate(AtomVector & _atoms, double _degrees);
		void Yrotate(AtomVector & _atoms, double _degrees);
		void Zrotate(AtomVector & _atoms, double _degrees);
		void rotate(AtomVector & _atoms, double _degrees, const CartesianPoint & _axisFromRotCenter, const CartesianPoint & _rotCenter=CartesianPoint(0.0, 0.0, 0.0));
		void rotate(AtomVector & _atoms, const Matrix & _rotMatrix, const CartesianPoint & _rotCenter=CartesianPoint(0.0, 0.0, 0.0));
		
		/*****************************************************
		 *  For the atom vector, the align and orient operations 
		 *  are applied to an external point _reference, and the
		 *  _atoms are moved according to the same transformation
		 *****************************************************/
		void align(AtomVector & _atoms, const CartesianPoint & _reference, const CartesianPoint & _target, const CartesianPoint & _rotCenter=CartesianPoint(0.0, 0.0, 0.0));
		void orient(AtomVector & _atoms, const CartesianPoint & _reference, const CartesianPoint & _target, const CartesianPoint & _axis1, const CartesianPoint & _axis2);


		/*******************************************************
		 *  For two atom vectors, align atoms via moving geometric centers on top, 
                 *  solving for a optimal rotation matrix.  Return value is a bool
                 *  which is true when alignment is proper.
		 *******************************************************/
		bool align(AtomVector &_align, AtomVector &_ref);
		bool align(AtomVector &_align, AtomVector &_ref, AtomVector &_moveable);
		/*
		double align(vector<Residue *> &_align, vector<Residue *> &_ref); // Align backbone atoms of each residue
		double align(vector<Residue *> &_align, vector<Residue *> &_ref,vector<Residue *> &_moveable); // Align backbone atoms of each residue
		double align(vector<Residue *> &_align, vector<Residue *> &_ref,AtomVector &_moveable); // Align backbone atoms of each residue
		*/

		Matrix createBasisTransformation(vector<vector<double> > &_basis1, vector<vector<double> > &_basis2);


		SphericalPoint transform(CartesianPoint &_p1);


		/*****************************************************
		 *  The object can store the history of all transofmations
		 *  as applied to a local frame.  If the history is saved,
		 *  the result of all recorded conformations can be applied
		 *  at once to any other atom or atom vector
		 *****************************************************/
		void setStoreTransformHistory(bool _flag);
		bool getStoreTransformHistory() const;
		void resetHistory();
		void applyHistory(Atom & _atom);
		void applyHistory(AtomVector & _atoms);

		Matrix getLastRotationMatrix() const;
		CartesianPoint getLastTranslation() const;

		/*******************************************************
                 * Functions to allow grid search (from Cinque Soto)
                 * RotatePdbAboutZYX rotates an AtomVector and the three
                 * localAxes (centered on the origin) based on the input
                 * degree rotations in local Z, Y, and X.  Center is
                 * not changed.
                 * TranslateRigidBody translates AtomVector and center
                 * my a specified distance.
                 *******************************************************/
		void RotatePdbAboutZYX(AtomVector & _theAtoms, CartesianPoint & _center, CartesianPoint & _localZ, CartesianPoint & _localY, CartesianPoint & _localX, double _RotationAlongZ, double _RotationAlongY, double _RotationAlongX);
		void TranslateRigidBodyPdbResidue(AtomVector & _theAtoms, CartesianPoint & _center, CartesianPoint & _TranslationVector, double _TranslationAmount);

	private:


		Quaternion q;
		
		Matrix lastRotMatrix;
		CartesianPoint lastTranslation;
		bool saveHistory_flag;

		map<string, CartesianPoint> frame;
};

// INLINE FUNCTIONS
inline void Transforms::setStoreTransformHistory(bool _flag) {saveHistory_flag = _flag;}
inline bool Transforms::getStoreTransformHistory() const {return saveHistory_flag;}
inline void Transforms::resetHistory() {frame["O"] = CartesianPoint(0.0, 0.0, 0.0); frame["X"] = CartesianPoint(1.0, 0.0, 0.0); frame["Y"] = CartesianPoint(0.0, 1.0, 0.0);}
inline Matrix Transforms::getLastRotationMatrix() const {return lastRotMatrix;}
inline CartesianPoint Transforms::getLastTranslation() const {return lastTranslation;}
#endif
